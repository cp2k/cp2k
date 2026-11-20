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
   USE fft_api,                         ONLY: fft_3d_local
   USE fft_tools,                       ONLY: BWFFT,&
                                              FFT_RADIX_NEXT,&
                                              FWFFT,&
                                              fft_allocate,&
                                              fft_deallocate,&
                                              fft_radix_operations
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_walltime
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
   PUBLIC :: k_grid_to_cell_panel_width
   PUBLIC :: lattice_fft_profile_type, lattice_fft_execute_pair_major, lattice_fft_prepare_geometry

   ! Optional local FFT diagnostics; enable explicitly for profiling builds.
   LOGICAL, PARAMETER, PUBLIC :: fft_profile_enabled = .FALSE.
   PUBLIC :: lattice_fft_real_sink_type

   ! Keep pair-major grid filling in cache-sized strips of AO elements and images.
   INTEGER, PARAMETER :: fft_pack_rows = 8, fft_pack_cells = 32
   ! Avoid OpenMP launch overhead unless a packing/extraction stage has enough values.
   INTEGER(KIND=int_8), PARAMETER :: fft_parallel_values = 32768_int_8
   ! Keep the default pair-major workspace bounded, as in the master lattice FFT.
   INTEGER(KIND=int_8), PARAMETER :: fft_workspace_bytes = 64_int_8*1024_int_8**2
   INTEGER(KIND=int_8), PARAMETER :: fft_panel_element_bytes = 48_int_8

   ! Optional local-panel diagnostics. Axis times include any plan creation.
   TYPE lattice_fft_profile_type
      REAL(KIND=dp) :: seconds(9) = 0.0_dp
      INTEGER(KIND=int_8) :: batches = 0_int_8, elements = 0_int_8
      INTEGER(KIND=int_8) :: transforms = 0_int_8
      INTEGER :: nfft(3) = 0, tile_min = HUGE(0), tile_max = 0
      LOGICAL :: pair_major = .FALSE., two_real = .FALSE.
   END TYPE lattice_fft_profile_type

   ! Reusable forward lattice geometry and bounded numerical FFT workspace.
   TYPE cell_to_k_grid_work_type
      LOGICAL :: regular = .FALSE.
      INTEGER :: nfft(3) = 0
      INTEGER, ALLOCATABLE, PRIVATE :: k_slot(:), cell_slot(:)
      COMPLEX(KIND=dp), ALLOCATABLE, PRIVATE :: cell_factor(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER, PRIVATE :: fft_in(:) => NULL(), fft_out(:) => NULL()
   END TYPE cell_to_k_grid_work_type

   ! Reusable inverse lattice maps, buffers and plans; the caller releases them after use.
   TYPE k_grid_to_cell_work_type
      LOGICAL, PRIVATE :: regular = .FALSE., uniform_weights = .FALSE.
      INTEGER :: nfft(3) = 0
      INTEGER, ALLOCATABLE :: k_slot(:), cell_slot(:), k_partner(:)
      INTEGER, ALLOCATABLE, PRIVATE :: cells(:, :)
      REAL(KIND=dp), ALLOCATABLE, PRIVATE :: xkp(:, :), weight(:)
      COMPLEX(KIND=dp), ALLOCATABLE, PRIVATE :: cell_phase(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER, PRIVATE :: fft_in(:) => NULL(), fft_out(:) => NULL()
   END TYPE k_grid_to_cell_work_type

   ! A real-output consumer may finish extraction directly in its next operator's layout.
   ! The FFT algorithm and its bounded panels remain shared with the array-output interface.
   TYPE, ABSTRACT :: lattice_fft_real_sink_type
      ! Whether each complex transform in the supplied panel packs two real fields.
      LOGICAL :: two_real = .FALSE.
   CONTAINS
      PROCEDURE(lattice_fft_consume_real), DEFERRED :: consume
   END TYPE lattice_fft_real_sink_type

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief Consume one inverse FFT panel for real output.
!> \param sink consumer storing the output or fusing a subsequent operation
!> \param first first flattened field index in the panel
!> \param nbatch number of physical fields in the panel
!> \param spectrum pair-major complex FFT output
!> \param cell_slot FFT slot of each requested cell
!> \param cell_phase inverse Bloch phase with any common weight
! **************************************************************************************************
      SUBROUTINE lattice_fft_consume_real(sink, first, nbatch, spectrum, cell_slot, cell_phase)
         IMPORT :: dp, int_8, lattice_fft_real_sink_type
         CLASS(lattice_fft_real_sink_type), INTENT(INOUT) :: sink
         INTEGER(KIND=int_8), INTENT(IN)                  :: first
         INTEGER, INTENT(IN)                              :: nbatch, cell_slot(:)
         COMPLEX(KIND=dp), INTENT(IN)                     :: spectrum(:, :), cell_phase(:)
      END SUBROUTINE lattice_fft_consume_real
   END INTERFACE

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
         IF (allow_arbitrary) THEN
            ! The PW grid radix list is a cost heuristic, not a capability limit.
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
!> \brief Prepare shared lattice FFT geometry in internal z,y,x storage order.
!>        An optional translated origin maps R to R-cell_origin. The output phase
!>        restores exp(2*pi*i*m*cell_origin/nkp_grid) for a forward transform;
!>        the physical k offset is returned separately for caller-specific factors.
!> \param xkp ...
!> \param nkp_grid ...
!> \param cells ...
!> \param grid_index ...
!> \param nfft ...
!> \param k_offset ...
!> \param cell_index ...
!> \param cell_linear ...
!> \param cell_origin ...
!> \param origin_phase ...
!> \param allow_incomplete ...
!> \param grid_linear ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION lattice_fft_prepare_geometry(xkp, nkp_grid, cells, grid_index, nfft, k_offset, &
                                                 cell_index, cell_linear, cell_origin, origin_phase, &
                                                 allow_incomplete, grid_linear) RESULT(compatible)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3), cells(:, :)
      INTEGER, INTENT(OUT)                               :: grid_index(:, :), nfft(3)
      REAL(KIND=dp), INTENT(OUT)                         :: k_offset(3)
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: cell_index(:, :), cell_linear(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: cell_origin(3)
      COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL            :: origin_phase(:)
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_incomplete
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL        :: grid_linear(:)

      INTEGER                                            :: d, icell, ik, INDEX(3), origin(3), &
                                                            physical_axis, SHAPE(3)
      REAL(KIND=dp)                                      :: arg

      compatible = regular_kpoint_grid(xkp, nkp_grid, grid_index, k_offset, allow_incomplete)
      IF (.NOT. compatible) RETURN
      compatible = lattice_fft_shape(nkp_grid, shape, allow_arbitrary=.TRUE.)
      IF (.NOT. compatible) RETURN

      nfft = SHAPE([3, 2, 1])
      grid_index = grid_index([3, 2, 1], :)
      IF (PRESENT(grid_linear)) THEN
         ALLOCATE (grid_linear(SIZE(xkp, 2)))
         DO ik = 1, SIZE(xkp, 2)
            grid_linear(ik) = grid_index(1, ik) + nfft(1)* &
                              (grid_index(2, ik) - 1 + nfft(2)*(grid_index(3, ik) - 1))
         END DO
      END IF
      origin = 0
      IF (PRESENT(cell_origin)) origin = cell_origin
      IF (PRESENT(cell_index)) ALLOCATE (cell_index(3, SIZE(cells, 2)))
      IF (PRESENT(cell_linear)) ALLOCATE (cell_linear(SIZE(cells, 2)))
      DO icell = 1, SIZE(cells, 2)
         DO d = 1, 3
            physical_axis = 4 - d
            INDEX(d) = 1 + MODULO(cells(physical_axis, icell) - origin(physical_axis), &
                                  nkp_grid(physical_axis))*(nfft(d)/nkp_grid(physical_axis))
         END DO
         IF (PRESENT(cell_index)) cell_index(:, icell) = index
         IF (PRESENT(cell_linear)) cell_linear(icell) = &
            INDEX(1) + nfft(1)*(INDEX(2) - 1 + nfft(2)*(INDEX(3) - 1))
      END DO
      IF (PRESENT(origin_phase)) THEN
         CPASSERT(SIZE(origin_phase) == SIZE(xkp, 2))
         DO ik = 1, SIZE(xkp, 2)
            arg = SUM(REAL(grid_index([3, 2, 1], ik) - 1, dp)*REAL(origin, dp)/REAL(nkp_grid, dp))
            origin_phase(ik) = EXP(gaussi*twopi*arg)
         END DO
      END IF

   END FUNCTION lattice_fft_prepare_geometry

! **************************************************************************************************
!> \brief Transform pair-major FFT input directly to pair-major output.
!>        Report the complete 3D transform under the first FFT stage in the profile.
!> \param dir ...
!> \param n ...
!> \param nbatch ...
!> \param buffer ...
!> \param scratch ...
!> \param profile ...
! **************************************************************************************************
   SUBROUTINE lattice_fft_execute_pair_major(dir, n, nbatch, buffer, scratch, profile)
      LOGICAL, INTENT(IN)                                :: dir
      INTEGER, INTENT(IN)                                :: n(3), nbatch
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: buffer(:), scratch(:)
      TYPE(lattice_fft_profile_type), INTENT(INOUT), &
         OPTIONAL                                        :: profile

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: swap(:)
      REAL(KIND=dp)                                      :: started

      IF (PRESENT(profile)) started = m_walltime()
      CALL fft_3d_local(dir, n, nbatch, .TRUE., .TRUE., buffer, scratch)
      IF (PRESENT(profile)) profile%seconds(5) = profile%seconds(5) + m_walltime() - started
      swap => buffer
      buffer => scratch
      scratch => swap
   END SUBROUTINE lattice_fft_execute_pair_major

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
!> \param max_workspace_bytes optional soft FFT workspace target; one element always fits the policy
!> \param accumulate add to values_k; abort on FFT execution failure to avoid partial accumulation
!> \param profile optional local-panel stage timings and batch counts
! **************************************************************************************************
   SUBROUTINE cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, values_k, &
                                 used_fft, deriv_direction, hmat, selected_kpoints, require_fft, work, cell_first, &
                                 max_workspace_bytes, accumulate, profile)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: values_k
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
      LOGICAL, INTENT(IN), OPTIONAL                      :: accumulate
      TYPE(lattice_fft_profile_type), INTENT(INOUT), &
         OPTIONAL                                        :: profile

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: input(:, :), output(:, :)
      INTEGER                                            :: handle, i, icell, ik, ik_out, image, &
                                                            image_first, matrix_column, &
                                                            matrix_row, nbatch, nout, row_first, &
                                                            slot, start, width
      INTEGER(KIND=int_8)                                :: budget, first, matrix_index, needed, &
                                                            ngrid, nvalues
      LOGICAL                                            :: add, success
      REAL(KIND=dp)                                      :: started
      TYPE(cell_to_k_grid_work_type), POINTER            :: plan
      TYPE(cell_to_k_grid_work_type), TARGET             :: local_work

      CALL timeset("cell_to_k_grid_fft", handle)
      add = .FALSE.
      IF (PRESENT(accumulate)) add = accumulate

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
         IF (PRESENT(profile)) started = m_walltime()
         CALL cell_to_k_grid_prepare(plan, index_to_cell, xkp, nkp_grid, deriv_direction, hmat)
         IF (PRESENT(profile)) profile%seconds(2) = profile%seconds(2) + m_walltime() - started
      END IF
      IF (.NOT. plan%regular .AND. PRESENT(require_fft)) THEN
         IF (require_fft) THEN
            CPABORT("Forced lattice FFT needs a supported regular grid or subset")
         END IF
      END IF
      success = .FALSE.
      IF (plan%regular) THEN
         start = 1
         IF (PRESENT(cell_first)) start = cell_first
         CPASSERT(start >= 1)
         IF (SIZE(values_rs, 3) > SIZE(plan%cell_factor) - start + 1) THEN
            CPABORT("Lattice FFT image range exceeds the prepared geometry")
         END IF
         ngrid = PRODUCT(INT(plan%nfft, int_8))
         ! A single pair may exceed this soft limit; otherwise keep workspaces bounded.
         budget = fft_workspace_bytes
         IF (PRESENT(max_workspace_bytes)) budget = MAX(48_int_8*ngrid, max_workspace_bytes)
         nvalues = SIZE(values_rs, 1, KIND=int_8)*SIZE(values_rs, 2, KIND=int_8)
         ! The FFT wrappers and allocated buffer sizes use default INTEGER.
         CPASSERT(ngrid <= HUGE(0))
         width = INT(MIN(nvalues, MAX(1_int_8, budget/(48_int_8*ngrid)), &
                         INT(HUGE(0), int_8)/ngrid))
         IF (width == 0) THEN
            success = .TRUE.
         ELSE
            needed = ngrid*width
            CPASSERT(needed <= HUGE(0))
            IF (PRESENT(profile)) THEN
               profile%nfft = plan%nfft
               started = m_walltime()
            END IF
            IF (ASSOCIATED(plan%fft_in)) THEN
               ! Reuse capacity while also allowing a smaller caller budget.
               IF (SIZE(plan%fft_in, KIND=int_8) < needed .OR. &
                   32_int_8*SIZE(plan%fft_in, KIND=int_8) + 16_int_8*needed > budget) THEN
                  CALL fft_deallocate(plan%fft_in)
                  CALL fft_deallocate(plan%fft_out)
               END IF
            END IF
            IF (.NOT. ASSOCIATED(plan%fft_in)) THEN
               CALL fft_allocate(plan%fft_in, [INT(needed)])
               CALL fft_allocate(plan%fft_out, [INT(needed)])
            END IF
            IF (PRESENT(profile)) THEN
               profile%pair_major = .TRUE.
               profile%seconds(3) = profile%seconds(3) + m_walltime() - started
            END IF
            ! Flatten both matrix dimensions so each FFT batch spans whole columns.
            DO first = 1_int_8, nvalues, INT(width, int_8)
               nbatch = INT(MIN(INT(width, int_8), nvalues - first + 1_int_8))
               IF (PRESENT(profile)) THEN
                  profile%batches = profile%batches + 1_int_8
                  profile%elements = profile%elements + INT(nbatch, int_8)
                  profile%tile_min = MIN(profile%tile_min, nbatch)
                  profile%tile_max = MAX(profile%tile_max, nbatch)
                  started = m_walltime()
               END IF
               input(1:nbatch, 1:INT(ngrid)) => plan%fft_in(1:INT(ngrid)*nbatch)
               input = z_zero
               ! Read adjacent AO elements together while retaining each element's image order.
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP IF(INT(nbatch, int_8)*SIZE(values_rs, 3, KIND=int_8) >= fft_parallel_values) &
!$OMP SHARED(input, plan, values_rs, first, nbatch, start) &
!$OMP PRIVATE(image_first, icell, image, slot, i, matrix_index, matrix_row, matrix_column)
               DO row_first = 1, nbatch, fft_pack_rows
                  DO image_first = 1, SIZE(values_rs, 3), fft_pack_cells
                     DO icell = image_first, MIN(image_first + fft_pack_cells - 1, SIZE(values_rs, 3))
                        image = start + icell - 1
                        slot = plan%cell_slot(image)
                        IF (SIZE(values_rs, 2) == 1) THEN
!$OMP SIMD
                           DO i = row_first, MIN(row_first + fft_pack_rows - 1, nbatch)
                              input(i, slot) = input(i, slot) + &
                                               plan%cell_factor(image)*values_rs(INT(first) + i - 1, 1, icell)
                           END DO
                        ELSE
!$OMP SIMD PRIVATE(matrix_index, matrix_row, matrix_column)
                           DO i = row_first, MIN(row_first + fft_pack_rows - 1, nbatch)
                              matrix_index = first + INT(i - 2, int_8)
                              matrix_row = INT(MODULO(matrix_index, INT(SIZE(values_rs, 1), int_8))) + 1
                              matrix_column = INT(matrix_index/INT(SIZE(values_rs, 1), int_8)) + 1
                              input(i, slot) = input(i, slot) + &
                                               plan%cell_factor(image)*values_rs(matrix_row, matrix_column, icell)
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
               IF (PRESENT(profile)) profile%seconds(4) = profile%seconds(4) + m_walltime() - started
               CALL lattice_fft_execute_pair_major(BWFFT, plan%nfft, nbatch, plan%fft_in, plan%fft_out, profile)
               IF (PRESENT(profile)) started = m_walltime()
               output(1:nbatch, 1:INT(ngrid)) => plan%fft_in(1:INT(ngrid)*nbatch)
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP IF(INT(nbatch, int_8)*SIZE(values_k, 3, KIND=int_8) >= fft_parallel_values) &
!$OMP SHARED(output, plan, values_k, selected_kpoints, first, nbatch, add) &
!$OMP PRIVATE(ik, i, matrix_index, matrix_row, matrix_column)
               DO ik_out = 1, SIZE(values_k, 3)
                  ik = ik_out
                  IF (PRESENT(selected_kpoints)) ik = selected_kpoints(ik_out)
                  IF (SIZE(values_k, 2) == 1) THEN
                     IF (add) THEN
                        values_k(INT(first):INT(first) + nbatch - 1, 1, ik_out) = &
                           values_k(INT(first):INT(first) + nbatch - 1, 1, ik_out) + &
                           output(:, plan%k_slot(ik))
                     ELSE
                        values_k(INT(first):INT(first) + nbatch - 1, 1, ik_out) = &
                           output(:, plan%k_slot(ik))
                     END IF
                  ELSE
                     DO i = 1, nbatch
                        matrix_index = first + INT(i - 2, int_8)
                        matrix_row = INT(MODULO(matrix_index, INT(SIZE(values_k, 1), int_8))) + 1
                        matrix_column = INT(matrix_index/INT(SIZE(values_k, 1), int_8)) + 1
                        IF (add) THEN
                           values_k(matrix_row, matrix_column, ik_out) = &
                              values_k(matrix_row, matrix_column, ik_out) + &
                              output(i, plan%k_slot(ik))
                        ELSE
                           values_k(matrix_row, matrix_column, ik_out) = &
                              output(i, plan%k_slot(ik))
                        END IF
                     END DO
                  END IF
               END DO
!$OMP END PARALLEL DO
               IF (PRESENT(profile)) profile%seconds(8) = profile%seconds(8) + m_walltime() - started
            END DO
            success = .TRUE.
         END IF
      END IF
      IF (.NOT. success) THEN
         IF (PRESENT(require_fft)) THEN
            IF (require_fft) THEN
               CPABORT("Forced lattice FFT execution failed")
            END IF
         END IF
         CALL direct_cell_to_k(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat, selected_kpoints, add)
      END IF
      IF (PRESENT(used_fft)) used_fft = success
      IF (.NOT. PRESENT(work)) THEN
         IF (PRESENT(profile)) started = m_walltime()
         CALL cell_to_k_grid_release(local_work)
         IF (PRESENT(profile)) profile%seconds(9) = profile%seconds(9) + m_walltime() - started
      END IF
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

      INTEGER                                            :: icell
      INTEGER, ALLOCATABLE                               :: grid_index(:, :)
      REAL(KIND=dp)                                      :: arg, cell_vector(3), k_offset(3)

      CALL cell_to_k_grid_release(work)
      ALLOCATE (grid_index(3, SIZE(xkp, 2)))
      work%regular = lattice_fft_prepare_geometry(xkp, nkp_grid, index_to_cell, grid_index, &
                                                  work%nfft, k_offset, cell_linear=work%cell_slot, &
                                                  grid_linear=work%k_slot, allow_incomplete=.TRUE.)
      IF (.NOT. work%regular) RETURN
      IF (PRESENT(deriv_direction)) THEN
         CPASSERT(PRESENT(hmat))
         CPASSERT(deriv_direction >= 1 .AND. deriv_direction <= 3)
      END IF
      ALLOCATE (work%cell_factor(SIZE(index_to_cell, 2)))
      DO icell = 1, SIZE(index_to_cell, 2)
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

      IF (ASSOCIATED(work%fft_in)) CALL fft_deallocate(work%fft_in)
      IF (ASSOCIATED(work%fft_out)) CALL fft_deallocate(work%fft_out)
      IF (ALLOCATED(work%k_slot)) DEALLOCATE (work%k_slot)
      IF (ALLOCATED(work%cell_slot)) DEALLOCATE (work%cell_slot, work%cell_factor)
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
!> \param work reusable maps, weights, buffers and backend plans
!> \param xkp reciprocal coordinates
!> \param nkp_grid candidate reciprocal-grid dimensions
!> \param index_to_cell requested lattice translations
!> \param weights final k weights, otherwise 1/number-of-points
!> \param allow_incomplete allow missing grid points, with zero values
!> \param cell_phase_forward unweighted forward Bloch phases, returned for a regular grid
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_prepare(work, xkp, nkp_grid, index_to_cell, weights, allow_incomplete, cell_phase_forward)

      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: weights
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_incomplete
      COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL            :: cell_phase_forward(:)

      INTEGER                                            :: half_shift(3), icell, ik, ncell, nkp, &
                                                            partner_index(3), partner_slot
      INTEGER, ALLOCATABLE                               :: grid_index(:, :), slot_to_k(:)
      REAL(KIND=dp), DIMENSION(3)                        :: k_offset

      CALL k_grid_to_cell_release(work)
      nkp = SIZE(xkp, 2)
      ncell = SIZE(index_to_cell, 2)
      IF (SIZE(xkp, 1) < 3 .OR. SIZE(index_to_cell, 1) < 3 .OR. nkp == 0) THEN
         CALL cp_abort(__LOCATION__, "Inconsistent reciprocal-grid mapping in inverse lattice FFT")
      END IF
      ALLOCATE (work%xkp(3, nkp), work%cells(3, ncell), work%weight(nkp), grid_index(3, nkp))
      work%xkp(:, :) = xkp(1:3, :)
      work%cells(:, :) = index_to_cell(1:3, :)
      work%weight(:) = 1.0_dp/REAL(nkp, dp)
      IF (PRESENT(weights)) THEN
         CPASSERT(SIZE(weights) == nkp)
         work%weight(:) = weights
      END IF
      work%regular = lattice_fft_prepare_geometry(xkp, nkp_grid, index_to_cell, grid_index, &
                                                  work%nfft, k_offset, cell_linear=work%cell_slot, &
                                                  grid_linear=work%k_slot, allow_incomplete=allow_incomplete)
      IF (.NOT. work%regular) RETURN
      ! For uniform k weights, apply the common scale once at extraction.
      work%uniform_weights = ALL(work%weight == work%weight(1))

      ! Real output depends only on S(k)=[w(k)A(k)+conjg(w(-k)A(-k))]/2.
      ! A complete inversion-closed mesh permits two such fields in one complex FFT.
      half_shift = NINT(2.0_dp*k_offset*REAL(nkp_grid, dp))
      IF (ALL(work%nfft == nkp_grid([3, 2, 1])) .AND. nkp == PRODUCT(nkp_grid) .AND. &
          ALL(ABS(2.0_dp*k_offset*REAL(nkp_grid, dp) - REAL(half_shift, dp)) <= 64.0_dp*EPSILON(1.0_dp))) THEN
         ALLOCATE (work%k_partner(nkp), slot_to_k(nkp))
         DO ik = 1, nkp
            slot_to_k(work%k_slot(ik)) = ik
         END DO
         DO ik = 1, nkp
            partner_index = 1 + MODULO(1 - grid_index(:, ik) - half_shift([3, 2, 1]), nkp_grid([3, 2, 1]))
            partner_slot = partner_index(1) + work%nfft(1)* &
                           (partner_index(2) - 1 + work%nfft(2)*(partner_index(3) - 1))
            work%k_partner(ik) = slot_to_k(partner_slot)
         END DO
      END IF

      ALLOCATE (work%cell_phase(ncell))
      IF (PRESENT(cell_phase_forward)) THEN
         CPASSERT(SIZE(cell_phase_forward) == ncell)
      END IF
      DO icell = 1, ncell
         work%cell_phase(icell) = EXP(-gaussi*twopi*SUM(k_offset*REAL(index_to_cell(1:3, icell), dp)))
         IF (PRESENT(cell_phase_forward)) cell_phase_forward(icell) = CONJG(work%cell_phase(icell))
         IF (work%uniform_weights) work%cell_phase(icell) = work%weight(1)*work%cell_phase(icell)
      END DO

   END SUBROUTINE k_grid_to_cell_prepare

! **************************************************************************************************
!> \brief Physical fields per inverse FFT panel, including two-real alignment.
!> \param work prepared inverse FFT geometry
!> \param nvalues number of physical fields
!> \param budget workspace limit in bytes
!> \param real_output whether the inverse uses the real-output path
!> \return panel width, zero for an empty input
! **************************************************************************************************
   FUNCTION k_grid_to_cell_panel_width(work, nvalues, budget, real_output) RESULT(width)
      TYPE(k_grid_to_cell_work_type), INTENT(IN)         :: work
      INTEGER(KIND=int_8), INTENT(IN)                    :: nvalues, budget
      LOGICAL, INTENT(IN)                                :: real_output
      INTEGER                                            :: width

      INTEGER(KIND=int_8)                                :: ngrid

      ngrid = PRODUCT(INT(work%nfft, int_8))
      CPASSERT(ngrid > 0 .AND. ngrid <= HUGE(0))
      width = INT(MIN(nvalues, MAX(1_int_8, budget/(fft_panel_element_bytes*ngrid)), INT(HUGE(0), int_8)/ngrid))
      IF (real_output .AND. ALLOCATED(work%k_partner) .AND. width > 1) width = 2*(width/2)
   END FUNCTION k_grid_to_cell_panel_width

! **************************************************************************************************
!> \brief Execute a matrix batch using the prepared inverse lattice transform.
!>        Call outside OpenMP worker regions: fft_tools owns its shared pool and planner.
!> \param work prepared transform; scratch is reused, never shared between concurrent callers
!> \param values_k matrices at the prepared reciprocal points
!> \param values_rs complex matrices at the prepared lattice translations; omit for values_real
!> \param used_fft whether this batch used FFT; false denotes a direct fallback or no result
!> \param require_fft abort on an unsupported mapping or FFT failure, without a direct fallback
!> \param allow_direct_fallback execute the generic direct sum when FFT is unavailable; defaults to true
!> \param values_real real part of the inverse transform, instead of values_rs
!> \param profile optional local-panel stage timings and batch counts
!> \param max_workspace_bytes optional soft FFT workspace target (default 64 MiB)
!> \param real_sink consume a real FFT panel directly, instead of materializing an output array
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_execute(work, values_k, values_rs, used_fft, require_fft, allow_direct_fallback, values_real, &
                                     profile, max_workspace_bytes, real_sink)

      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: values_k
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(OUT), OPTIONAL                           :: values_rs
      LOGICAL, INTENT(OUT), OPTIONAL                     :: used_fft
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_fft, allow_direct_fallback
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: values_real(:, :, :)
      TYPE(lattice_fft_profile_type), INTENT(INOUT), &
         OPTIONAL                                        :: profile
      INTEGER(KIND=int_8), INTENT(IN), OPTIONAL          :: max_workspace_bytes

      CLASS(lattice_fft_real_sink_type), INTENT(INOUT), &
         OPTIONAL                                        :: real_sink

      INTEGER                                            :: fft_count, fft_width, handle, nbatch, n_output, &
                                                            output_shape(3), width
      INTEGER(KIND=int_8)                                :: budget, first, ngrid, nvalues
      LOGICAL                                            :: pack_real, real_output
      REAL(KIND=dp)                                      :: started

      CALL timeset("k_grid_to_cell_fft", handle)
      CPASSERT(ALLOCATED(work%xkp))
      n_output = MERGE(1, 0, PRESENT(values_rs)) + MERGE(1, 0, PRESENT(values_real)) + MERGE(1, 0, PRESENT(real_sink))
      CPASSERT(n_output == 1)
      IF (PRESENT(values_rs)) output_shape = SHAPE(values_rs)
      IF (PRESENT(values_real)) output_shape = SHAPE(values_real)
      IF (PRESENT(real_sink)) THEN
         output_shape = [SIZE(values_k, 1), SIZE(values_k, 2), SIZE(work%cells, 2)]
      END IF
      IF (SIZE(values_k, 3) /= SIZE(work%xkp, 2) .OR. &
          ANY(output_shape /= [SIZE(values_k, 1), SIZE(values_k, 2), SIZE(work%cells, 2)])) THEN
         CALL cp_abort(__LOCATION__, "Inconsistent matrix batches in inverse lattice FFT")
      END IF

      IF (PRESENT(require_fft)) THEN
         IF (require_fft .AND. .NOT. work%regular) THEN
            CPABORT("Forced inverse lattice FFT needs a supported regular grid")
         END IF
      END IF
      IF (work%regular) THEN
         ngrid = PRODUCT(INT(work%nfft, int_8))
         nvalues = SIZE(values_k, 1, KIND=int_8)*SIZE(values_k, 2, KIND=int_8)
         budget = fft_workspace_bytes
         IF (PRESENT(max_workspace_bytes)) budget = max_workspace_bytes
         real_output = PRESENT(values_real) .OR. PRESENT(real_sink)
         width = k_grid_to_cell_panel_width(work, nvalues, budget, real_output)
         budget = MAX(fft_panel_element_bytes*ngrid, budget)
         IF (width == 0) THEN
            IF (PRESENT(used_fft)) used_fft = .TRUE.
            CALL timestop(handle)
            RETURN
         END IF
         pack_real = real_output .AND. ALLOCATED(work%k_partner) .AND. width > 1
         fft_width = width
         IF (pack_real) fft_width = width/2
         CPASSERT(ngrid*width <= HUGE(0))
         IF (PRESENT(profile)) profile%pair_major = .TRUE.
         IF (PRESENT(profile)) THEN
            profile%two_real = pack_real
            profile%nfft = work%nfft
            started = m_walltime()
         END IF
         IF (ASSOCIATED(work%fft_in)) THEN
            IF (SIZE(work%fft_in, KIND=int_8) < ngrid*fft_width .OR. &
                32_int_8*SIZE(work%fft_in, KIND=int_8) + 16_int_8*ngrid*fft_width > budget) THEN
               CALL fft_deallocate(work%fft_in)
               CALL fft_deallocate(work%fft_out)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(work%fft_in)) THEN
            CALL fft_allocate(work%fft_in, [INT(ngrid)*fft_width])
            CALL fft_allocate(work%fft_out, [INT(ngrid)*fft_width])
         END IF
         IF (PRESENT(profile)) profile%seconds(3) = profile%seconds(3) + m_walltime() - started
         ! Flatten both matrix dimensions so each FFT batch spans whole columns.
         DO first = 1_int_8, nvalues, INT(width, int_8)
            nbatch = INT(MIN(INT(width, int_8), nvalues - first + 1_int_8))
            fft_count = nbatch
            IF (pack_real) fft_count = (nbatch + 1)/2
            IF (PRESENT(profile)) THEN
               profile%batches = profile%batches + 1_int_8
               profile%elements = profile%elements + INT(nbatch, int_8)
               profile%transforms = profile%transforms + INT(fft_count, int_8)
               profile%tile_min = MIN(profile%tile_min, nbatch)
               profile%tile_max = MAX(profile%tile_max, nbatch)
               started = m_walltime()
            END IF
            CALL k_grid_to_cell_pair_major_batch(work, values_k, first, nbatch, values_rs, values_real, &
                                                 profile, pack_real, real_sink)
         END DO
      END IF
      IF (.NOT. work%regular) THEN
         IF (PRESENT(real_sink)) THEN
            CPABORT('A lattice FFT consumer requires successful FFT execution')
         END IF
         IF (PRESENT(require_fft)) THEN
            IF (require_fft) THEN
               CPABORT("Forced inverse lattice FFT execution failed")
            END IF
         END IF
         IF (PRESENT(allow_direct_fallback)) THEN
            IF (.NOT. allow_direct_fallback) THEN
               IF (PRESENT(values_rs)) values_rs = z_zero
               IF (PRESENT(values_real)) values_real = 0.0_dp
               IF (PRESENT(used_fft)) used_fft = .FALSE.
               CALL timestop(handle)
               RETURN
            END IF
         END IF
         CALL direct_k_to_cell(values_k, work%xkp, work%cells, values_rs, work%weight, values_real)
      END IF
      IF (PRESENT(used_fft)) used_fft = work%regular
      CALL timestop(handle)

   END SUBROUTINE k_grid_to_cell_execute

! **************************************************************************************************
!> \brief One K2R FFT panel with pair-major buffers at both ends.
!>        Mapping from k-point order to FFT slot is performed by contiguous pair strips.
!> \param work ...
!> \param values_k ...
!> \param first ...
!> \param nbatch ...
!> \param values_rs ...
!> \param values_real ...
!> \param profile ...
!> \param pack_real pack two independently symmetrized real-output transforms
!> \param real_sink real-output consumer
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_pair_major_batch(work, values_k, first, nbatch, values_rs, values_real, profile, &
                                              pack_real, real_sink)
      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work
      COMPLEX(KIND=dp), INTENT(IN)                       :: values_k(:, :, :)
      INTEGER(KIND=int_8), INTENT(IN)                    :: first
      INTEGER, INTENT(IN)                                :: nbatch
      COMPLEX(KIND=dp), INTENT(INOUT), OPTIONAL          :: values_rs(:, :, :)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: values_real(:, :, :)
      TYPE(lattice_fft_profile_type), INTENT(INOUT), &
         OPTIONAL                                        :: profile
      LOGICAL, INTENT(IN)                                :: pack_real

      CLASS(lattice_fft_real_sink_type), INTENT(INOUT), &
         OPTIONAL                                        :: real_sink

      COMPLEX(KIND=dp)                                   :: field_a, field_b, phase, spectrum
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: input(:, :), output(:, :)
      INTEGER                                            :: fft_count, i, icell, ik, ipartner, &
                                                            matrix_column, matrix_row, &
                                                            other_column, other_row, partner_slot, &
                                                            slot
      INTEGER(KIND=int_8)                                :: matrix_index, ngrid
      REAL(KIND=dp)                                      :: started

      ngrid = PRODUCT(INT(work%nfft, int_8))
      fft_count = nbatch
      IF (pack_real) fft_count = (nbatch + 1)/2
      input(1:fft_count, 1:INT(ngrid)) => work%fft_in(1:fft_count*INT(ngrid))
      IF (PRESENT(profile)) started = m_walltime()
      ! Missing k points and FFT padding must be zero-filled.
      IF (SIZE(values_k, 3, KIND=int_8) /= ngrid) input = z_zero
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP IF(INT(nbatch, int_8)*SIZE(values_k, 3, KIND=int_8) >= fft_parallel_values) &
!$OMP SHARED(input, work, values_k, first, nbatch, fft_count, pack_real) &
!$OMP PRIVATE(ik, slot, i, matrix_index, matrix_row, matrix_column, ipartner, partner_slot) &
!$OMP PRIVATE(other_row, other_column, field_a, field_b)
      DO ik = 1, SIZE(values_k, 3)
         slot = work%k_slot(ik)
         IF (pack_real) THEN
            ipartner = work%k_partner(ik)
            IF (ipartner < ik) CYCLE
            partner_slot = work%k_slot(ipartner)
            ! Visit each +/-k pair once, so every original input is read once.
            ! Keep the factor 1/2 in extraction, alongside the common Bloch phase.
            IF (SIZE(values_k, 2) == 1) THEN
               !$OMP SIMD PRIVATE(matrix_row, field_a, field_b)
               DO i = 1, fft_count
                  matrix_row = INT(first) + 2*i - 2
                  IF (work%uniform_weights) THEN
                     field_a = values_k(matrix_row, 1, ik) + CONJG(values_k(matrix_row, 1, ipartner))
                  ELSE
                     field_a = work%weight(ik)*values_k(matrix_row, 1, ik) + &
                               CONJG(work%weight(ipartner)*values_k(matrix_row, 1, ipartner))
                  END IF
                  field_b = z_zero
                  IF (2*i <= nbatch) THEN
                     IF (work%uniform_weights) THEN
                        field_b = values_k(matrix_row + 1, 1, ik) + CONJG(values_k(matrix_row + 1, 1, ipartner))
                     ELSE
                        field_b = work%weight(ik)*values_k(matrix_row + 1, 1, ik) + &
                                  CONJG(work%weight(ipartner)*values_k(matrix_row + 1, 1, ipartner))
                     END IF
                  END IF
                  input(i, slot) = field_a + gaussi*field_b
                  IF (ipartner /= ik) input(i, partner_slot) = CONJG(field_a) + gaussi*CONJG(field_b)
               END DO
            ELSE
               DO i = 1, fft_count
                  matrix_index = first + INT(2*i - 3, int_8)
                  matrix_row = INT(MODULO(matrix_index, INT(SIZE(values_k, 1), int_8))) + 1
                  matrix_column = INT(matrix_index/INT(SIZE(values_k, 1), int_8)) + 1
                  IF (work%uniform_weights) THEN
                     field_a = values_k(matrix_row, matrix_column, ik) + &
                               CONJG(values_k(matrix_row, matrix_column, ipartner))
                  ELSE
                     field_a = work%weight(ik)*values_k(matrix_row, matrix_column, ik) + &
                               CONJG(work%weight(ipartner)*values_k(matrix_row, matrix_column, ipartner))
                  END IF
                  field_b = z_zero
                  IF (2*i <= nbatch) THEN
                     matrix_index = matrix_index + 1_int_8
                     other_row = INT(MODULO(matrix_index, INT(SIZE(values_k, 1), int_8))) + 1
                     other_column = INT(matrix_index/INT(SIZE(values_k, 1), int_8)) + 1
                     IF (work%uniform_weights) THEN
                        field_b = values_k(other_row, other_column, ik) + &
                                  CONJG(values_k(other_row, other_column, ipartner))
                     ELSE
                        field_b = work%weight(ik)*values_k(other_row, other_column, ik) + &
                                  CONJG(work%weight(ipartner)*values_k(other_row, other_column, ipartner))
                     END IF
                  END IF
                  input(i, slot) = field_a + gaussi*field_b
                  IF (ipartner /= ik) input(i, partner_slot) = CONJG(field_a) + gaussi*CONJG(field_b)
               END DO
            END IF
         ELSE IF (SIZE(values_k, 2) == 1) THEN
            IF (work%uniform_weights) THEN
               input(:, slot) = values_k(INT(first):INT(first) + nbatch - 1, 1, ik)
            ELSE
               input(:, slot) = work%weight(ik)*values_k(INT(first):INT(first) + nbatch - 1, 1, ik)
            END IF
         ELSE
            DO i = 1, nbatch
               matrix_index = first + INT(i - 2, int_8)
               matrix_row = INT(MODULO(matrix_index, INT(SIZE(values_k, 1), int_8))) + 1
               matrix_column = INT(matrix_index/INT(SIZE(values_k, 1), int_8)) + 1
               IF (work%uniform_weights) THEN
                  input(i, slot) = values_k(matrix_row, matrix_column, ik)
               ELSE
                  input(i, slot) = work%weight(ik)*values_k(matrix_row, matrix_column, ik)
               END IF
            END DO
         END IF
      END DO
!$OMP END PARALLEL DO
      IF (PRESENT(profile)) profile%seconds(4) = profile%seconds(4) + m_walltime() - started
      CALL lattice_fft_execute_pair_major(FWFFT, work%nfft, fft_count, work%fft_in, work%fft_out, profile)
      output(1:fft_count, 1:INT(ngrid)) => work%fft_in(1:fft_count*INT(ngrid))
      IF (PRESENT(profile)) started = m_walltime()
      IF (PRESENT(real_sink)) THEN
         real_sink%two_real = pack_real
         CALL real_sink%consume(first, nbatch, output, work%cell_slot, work%cell_phase)
         IF (PRESENT(profile)) profile%seconds(8) = profile%seconds(8) + m_walltime() - started
         RETURN
      END IF
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP IF(INT(nbatch, int_8)*SIZE(work%cells, 2, KIND=int_8) >= fft_parallel_values) &
!$OMP SHARED(output, work, values_rs, values_real, values_k, first, nbatch, fft_count, pack_real) &
!$OMP PRIVATE(icell, slot, i, matrix_index, matrix_row, matrix_column, phase, spectrum)
      DO icell = 1, SIZE(work%cells, 2)
         slot = work%cell_slot(icell)
         IF (pack_real) THEN
            phase = 0.5_dp*work%cell_phase(icell)
            IF (SIZE(values_k, 2) == 1) THEN
               !$OMP SIMD PRIVATE(matrix_row, spectrum)
               DO i = 1, fft_count
                  matrix_row = INT(first) + 2*i - 2
                  spectrum = phase*output(i, slot)
                  values_real(matrix_row, 1, icell) = REAL(spectrum, dp)
                  IF (2*i <= nbatch) values_real(matrix_row + 1, 1, icell) = AIMAG(spectrum)
               END DO
            ELSE
               DO i = 1, nbatch
                  matrix_index = first + INT(i - 2, int_8)
                  matrix_row = INT(MODULO(matrix_index, INT(SIZE(values_k, 1), int_8))) + 1
                  matrix_column = INT(matrix_index/INT(SIZE(values_k, 1), int_8)) + 1
                  spectrum = phase*output((i + 1)/2, slot)
                  values_real(matrix_row, matrix_column, icell) = &
                     MERGE(REAL(spectrum, dp), AIMAG(spectrum), MODULO(i, 2) == 1)
               END DO
            END IF
            CYCLE
         END IF
         IF (SIZE(values_k, 2) == 1) THEN
            IF (PRESENT(values_rs)) THEN
               values_rs(INT(first):INT(first) + nbatch - 1, 1, icell) = &
                  work%cell_phase(icell)*output(:, slot)
            ELSE
               values_real(INT(first):INT(first) + nbatch - 1, 1, icell) = &
                  REAL(work%cell_phase(icell)*output(:, slot), dp)
            END IF
         ELSE
            DO i = 1, nbatch
               matrix_index = first + INT(i - 2, int_8)
               matrix_row = INT(MODULO(matrix_index, INT(SIZE(values_k, 1), int_8))) + 1
               matrix_column = INT(matrix_index/INT(SIZE(values_k, 1), int_8)) + 1
               IF (PRESENT(values_rs)) THEN
                  values_rs(matrix_row, matrix_column, icell) = work%cell_phase(icell)*output(i, slot)
               ELSE
                  values_real(matrix_row, matrix_column, icell) = REAL(work%cell_phase(icell)*output(i, slot), dp)
               END IF
            END DO
         END IF
      END DO
!$OMP END PARALLEL DO
      IF (PRESENT(profile)) profile%seconds(8) = profile%seconds(8) + m_walltime() - started
   END SUBROUTINE k_grid_to_cell_pair_major_batch

! **************************************************************************************************
!> \brief Release inverse lattice plans before their buffers and shared geometry.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_release(work)

      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work

      IF (ASSOCIATED(work%fft_in)) CALL fft_deallocate(work%fft_in)
      IF (ASSOCIATED(work%fft_out)) CALL fft_deallocate(work%fft_out)
      IF (ALLOCATED(work%xkp)) DEALLOCATE (work%xkp, work%cells, work%weight)
      IF (ALLOCATED(work%k_slot)) DEALLOCATE (work%k_slot)
      IF (ALLOCATED(work%k_partner)) DEALLOCATE (work%k_partner)
      IF (ALLOCATED(work%cell_slot)) DEALLOCATE (work%cell_slot, work%cell_phase)
      work%regular = .FALSE.
      work%uniform_weights = .FALSE.

   END SUBROUTINE k_grid_to_cell_release

! **************************************************************************************************
!> \brief Direct inverse phase sum for arbitrary reciprocal-point lists.
!> \param values_k ...
!> \param xkp ...
!> \param index_to_cell ...
!> \param values_rs ...
!> \param weights optional k weights, otherwise 1/number-of-points
!> \param values_real real part of the inverse transform, instead of values_rs
! **************************************************************************************************
   SUBROUTINE direct_k_to_cell(values_k, xkp, index_to_cell, values_rs, weights, values_real)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: values_k
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(OUT), OPTIONAL                           :: values_rs
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: weights
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: values_real(:, :, :)

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik
      REAL(KIND=dp)                                      :: arg, normalization

      IF (PRESENT(values_rs)) values_rs = z_zero
      IF (PRESENT(values_real)) values_real = 0.0_dp
      normalization = 1.0_dp/REAL(SIZE(xkp, 2), KIND=dp)
      DO icell = 1, SIZE(index_to_cell, 2)
         DO ik = 1, SIZE(xkp, 2)
            IF (PRESENT(weights)) normalization = weights(ik)
            arg = SUM(xkp(1:3, ik)*REAL(index_to_cell(1:3, icell), KIND=dp))
            factor = normalization*EXP(-gaussi*twopi*arg)
            IF (PRESENT(values_rs)) THEN
               values_rs(:, :, icell) = values_rs(:, :, icell) + factor*values_k(:, :, ik)
            ELSE
               values_real(:, :, icell) = values_real(:, :, icell) + REAL(factor*values_k(:, :, ik), dp)
            END IF
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
!> \param accumulate add to the existing output
! **************************************************************************************************
   SUBROUTINE direct_cell_to_k(values_rs, index_to_cell, xkp, values_k, &
                               deriv_direction, hmat, selected_kpoints, accumulate)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: values_k
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: selected_kpoints
      LOGICAL, INTENT(IN)                                :: accumulate

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik, ik_out
      REAL(KIND=dp)                                      :: arg
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector

      IF (.NOT. accumulate) values_k = z_zero
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

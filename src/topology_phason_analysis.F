!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native streaming mesh analysis with physical endpoint sewing.
! **************************************************************************************************
MODULE topology_phason_analysis
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE iso_fortran_env,                 ONLY: int64
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE topology_curvature,              ONLY: curvature_density,&
                                              link_plaquette,&
                                              polar_link
   USE topology_snapshot,               ONLY: check_snapshot,&
                                              check_snapshot_seam,&
                                              read_snapshot,&
                                              snapshot_overlap,&
                                              snapshot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_phason_analysis'
   PRIVATE
   TYPE :: phason_result_type
      INTEGER :: dimension = 0, rank = 0, SHAPE(4) = 0, nat = 0, point = 0
      REAL(KIND=dp) :: value = 0.0_dp, origin(3) = 0.0_dp, directions(3, 4) = 0.0_dp
      CHARACTER(LEN=default_path_length) :: first_file = ''
   END TYPE phason_result_type
   PUBLIC :: analyze_mesh, phason_result_type, compare_refinements

CONTAINS

! **************************************************************************************************
!> \brief Flatten a zero-based index with the first coordinate varying fastest.
!> \param index Zero-based coordinate along each axis
!> \param shape Number of entries along each axis
!> \return One-based linear index
! **************************************************************************************************
   INTEGER FUNCTION flat(index, shape) RESULT(value)
      INTEGER, INTENT(IN)                                :: index(:), SHAPE(:)

      INTEGER                                            :: d, stride

      value = 1
      stride = 1
      DO d = 1, SIZE(shape)
         value = value + INDEX(d)*stride
         stride = stride*SHAPE(d)
      END DO
   END FUNCTION flat

! **************************************************************************************************
!> \brief Inverse of flat; first coordinate varies fastest.
!> \param value One-based linear index
!> \param shape Number of entries along each axis
!> \param index Zero-based coordinate along each axis
! **************************************************************************************************
   SUBROUTINE unflat(value, shape, index)
      INTEGER, INTENT(IN)                                :: value, SHAPE(:)
      INTEGER, INTENT(OUT)                               :: index(:)

      INTEGER                                            :: d, remainder

      remainder = value - 1
      DO d = 1, SIZE(shape)
         INDEX(d) = MOD(remainder, SHAPE(d))
         remainder = remainder/SHAPE(d)
      END DO
   END SUBROUTINE unflat

! **************************************************************************************************
!> \brief Build physical links or read reference links; stream plaquettes through an anonymous file.
!> \param path Manifest or reference-link filename
!> \param physical True for physical STATE_EXPORT frames, false for reference links
!> \param report Mesh metadata and unrounded Chern estimate for refinement comparison
!> \param io_unit Unit receiving mesh diagnostics
! **************************************************************************************************
   SUBROUTINE analyze_mesh(path, physical, report, io_unit)
      CHARACTER(LEN=*), INTENT(IN) :: path
      LOGICAL, INTENT(IN) :: physical
      TYPE(phason_result_type), INTENT(OUT) :: report
      INTEGER, INTENT(IN) :: io_unit
      CHARACTER(LEN=256) :: header
      CHARACTER(LEN=default_path_length), ALLOCATABLE :: files(:)
      INTEGER :: input, store, stat, dim, rank, nvertex, nframe, nat, vertex, axis, mu, nu, p, i, j, na, nb, nbase, idx, record_size
      INTEGER, ALLOCATABLE :: shape(:), extended(:), INDEX(:), next(:), zero(:), points(:), permutation(:, :)
      LOGICAL, ALLOCATABLE :: checked(:)
      REAL(KIND=dp), ALLOCATABLE :: offsets(:), directions(:, :)
      REAL(KIND=dp) :: origin(3), ka(3), kb(3), tol_metric, tol_gap, tol_sv, phase_limit, minimum, local_min, &
                       max_metric, min_gap, err, gap, value, contribution, max_phase, local_phase, seam_error, &
                       valence_max, conduction_min, energy_shift, compensation, y, t
      COMPLEX(KIND=dp), ALLOCATABLE :: raw(:, :), link(:, :), sewing(:, :), u(:, :, :), pqs(:, :, :)
      TYPE(snapshot_type), ALLOCATABLE :: a, b, base
      LOGICAL :: common_reference, prefix
      INTEGER(KIND=int64) :: total

      CALL open_file(path, unit_number=input, file_status='OLD', file_action='READ')
      READ (input, '(A)', IOSTAT=stat) header
      CPASSERT(stat == 0)
      IF (physical) THEN
         IF (TRIM(header) /= 'CP2K_PHASON_MESH 1') THEN
            CPABORT('Expected CP2K_PHASON_MESH 1')
         END IF
      ELSE
         IF (TRIM(header) /= 'CP2K_TOPOLOGY_LINKS 1') THEN
            CPABORT('Expected CP2K_TOPOLOGY_LINKS 1')
         END IF
      END IF
      READ (input, *, IOSTAT=stat) dim, rank
      CPASSERT(stat == 0)
      IF ((dim /= 2 .AND. dim /= 4) .OR. rank < 1) THEN
         CPABORT('Invalid mesh dimension or rank')
      END IF
      ALLOCATE (SHAPE(dim), extended(dim), INDEX(dim), next(dim), zero(dim))
      READ (input, *, IOSTAT=stat) shape
      CPASSERT(stat == 0)
      IF (ANY(shape < 3)) THEN
         CPABORT('At least three samples per mesh axis are required')
      END IF
      total = 1_int64
      DO axis = 1, dim
         IF (total > INT(HUGE(1), int64)/INT(SHAPE(axis) + 1, int64)) THEN
            CPABORT('Mesh index overflow')
         END IF
         total = total*INT(SHAPE(axis) + 1, int64)
      END DO
      extended(:) = shape + 1
      nvertex = PRODUCT(shape)
      nframe = INT(total)
      report%dimension = dim
      report%rank = rank
      report%shape(1:dim) = shape
      READ (input, *, IOSTAT=stat) tol_metric, tol_gap, tol_sv, phase_limit
      CPASSERT(stat == 0)
      IF (.NOT. ALL(ieee_is_finite([tol_metric, tol_gap, tol_sv, phase_limit]))) THEN
         CPABORT('Nonfinite tolerances')
      END IF
      IF (MIN(tol_metric, tol_gap, tol_sv, phase_limit) <= 0.0_dp .OR. MAX(tol_metric, tol_sv) >= 1.0_dp .OR. &
          phase_limit >= ACOS(-1.0_dp)) THEN
         CPABORT('Invalid mesh tolerances')
      END IF
      ALLOCATE (raw(rank, rank), link(rank, rank), sewing(rank, rank), u(rank, rank, 4), pqs(rank, rank, dim*(dim - 1)/2))
      ! Anonymous scratch storage avoids persistent, stale link caches.
      CALL open_file('', unit_number=store, file_status='SCRATCH', file_access='STREAM', &
                     file_form='UNFORMATTED', file_action='READWRITE')
      INQUIRE (IOLENGTH=record_size) link
      minimum = HUGE(1.0_dp)
      max_metric = 0.0_dp
      min_gap = HUGE(1.0_dp)
      seam_error = 0.0_dp
      valence_max = -HUGE(1.0_dp)
      conduction_min = HUGE(1.0_dp)
      common_reference = .FALSE.
      prefix = .TRUE.
      IF (physical) THEN
         ALLOCATE (a, b, base)
         ALLOCATE (files(nframe), points(nframe), offsets(nframe), checked(nframe), directions(3, dim))
         checked(:) = .FALSE.
         READ (input, *, IOSTAT=stat) origin
         CPASSERT(stat == 0)
         IF (.NOT. ALL(ieee_is_finite(origin))) THEN
            CPABORT('Nonfinite mesh origin')
         END IF
         DO axis = 1, dim
            READ (input, *, IOSTAT=stat) directions(:, axis)
            CPASSERT(stat == 0)
         END DO
         IF (.NOT. ALL(ieee_is_finite(directions))) THEN
            CPABORT('Nonfinite reciprocal periods')
         END IF
         IF (MAXVAL(ABS(directions - ANINT(directions))) > 1.e-10_dp) THEN
            CPABORT('Noninteger reciprocal periods')
         END IF
         READ (input, *, IOSTAT=stat) nat, common_reference
         CPASSERT(stat == 0)
         IF (nat < 1) THEN
            CPABORT('Invalid mesh atom count')
         END IF
         ALLOCATE (permutation(nat, dim))
         DO axis = 1, dim
            READ (input, *, IOSTAT=stat) permutation(:, axis)
            CPASSERT(stat == 0)
            IF (ANY(permutation(:, axis) < 1) .OR. ANY(permutation(:, axis) > nat)) THEN
               CPABORT('Invalid seam permutation')
            END IF
            DO idx = 1, nat
               IF (COUNT(permutation(:, axis) == idx) /= 1) THEN
                  CPABORT('Non-bijective seam map')
               END IF
            END DO
         END DO
         DO mu = 1, dim
            DO nu = mu + 1, dim
               IF (ANY(permutation(permutation(:, mu), nu) /= permutation(permutation(:, nu), mu))) THEN
                  CPABORT('The declared torus atom permutations do not commute')
               END IF
            END DO
         END DO
         DO idx = 1, nframe
            READ (input, *, IOSTAT=stat) files(idx), points(idx), offsets(idx)
            CPASSERT(stat == 0)
            IF (points(idx) < 1 .OR. .NOT. ieee_is_finite(offsets(idx))) THEN
               CPABORT('Invalid frame record')
            END IF
         END DO
         report%nat = nat
         report%origin(:) = origin
         report%directions(:, 1:dim) = directions
         report%first_file = files(1)
         report%point = points(1)
      END IF
      DO vertex = 1, nvertex
         CALL unflat(vertex, shape, index)
         IF (physical) THEN
            na = flat(index, extended)
            CALL frame(na, a)
            ka(:) = origin + MATMUL(directions, REAL(index, dp)/REAL(shape, dp))
         END IF
         DO axis = 1, dim
            IF (physical) THEN
               next(:) = index
               next(axis) = next(axis) + 1
               nb = flat(next, extended)
               CALL frame(nb, b)
               kb(:) = origin + MATMUL(directions, REAL(next, dp)/REAL(shape, dp))
               CALL snapshot_overlap(a, b, ka, kb, raw)
               IF (next(axis) == SHAPE(axis)) THEN
                  zero(:) = next
                  zero(axis) = 0
                  nbase = flat(zero, extended)
                  CALL frame(nbase, base)
                  CALL check_snapshot_seam(base, b, permutation(:, axis), tol_metric)
                  CALL snapshot_overlap(b, base, kb, kb, sewing)
                  ! Both endpoint frames must span the same physical subspace, not merely have full rank.
                  link(:, :) = MATMUL(CONJG(TRANSPOSE(sewing)), sewing)
                  DO idx = 1, rank
                     link(idx, idx) = link(idx, idx) - 1.0_dp
                  END DO
                  err = MAXVAL(ABS(link))
                  seam_error = MAX(seam_error, err)
                  IF (err > tol_metric) THEN
                     CPABORT('Endpoint selected subspaces do not close')
                  END IF
                  CALL polar_link(sewing, link, local_min, stat, tol_sv)
                  IF (stat /= 0) THEN
                     CPABORT('Invalid endpoint sewing link')
                  END IF
                  raw(:, :) = MATMUL(raw, link)
               END IF
            ELSE
               DO j = 1, rank
                  DO i = 1, rank
                     READ (input, *, IOSTAT=stat) ka(1:2)
                     CPASSERT(stat == 0)
                     raw(i, j) = CMPLX(ka(1), ka(2), dp)
                  END DO
               END DO
            END IF
            CALL polar_link(raw, link, local_min, stat, tol_sv)
            IF (stat /= 0) THEN
               CPABORT('Singular/nonfinite mesh link; refine or check the subspace')
            END IF
            minimum = MIN(minimum, local_min)
            WRITE (store, IOSTAT=stat) link
            CPASSERT(stat == 0)
         END DO
      END DO
      DO
         READ (input, '(A)', IOSTAT=stat) header
         IF (stat < 0) EXIT
         IF (stat /= 0 .OR. LEN_TRIM(header) /= 0) THEN
            CPABORT('Unexpected trailing mesh data')
         END IF
      END DO
      CALL close_file(input)
      value = 0.0_dp
      compensation = 0.0_dp
      max_phase = 0.0_dp
      DO vertex = 1, nvertex
         CALL unflat(vertex, shape, index)
         p = 0
         DO mu = 1, dim
            DO nu = mu + 1, dim
               p = p + 1
               CALL read_link(vertex, mu, u(:, :, 1))
               next(:) = index
               next(mu) = MOD(next(mu) + 1, SHAPE(mu))
               CALL read_link(flat(next, shape), nu, u(:, :, 2))
               next(:) = index
               next(nu) = MOD(next(nu) + 1, SHAPE(nu))
               CALL read_link(flat(next, shape), mu, u(:, :, 3))
               CALL read_link(vertex, nu, u(:, :, 4))
               CALL link_plaquette(u(:, :, 1), u(:, :, 2), u(:, :, 3), u(:, :, 4), pqs(:, :, p))
            END DO
         END DO
         CALL curvature_density(pqs, dim, phase_limit, contribution, local_phase, stat)
         IF (stat /= 0) THEN
            CPABORT('Unresolved/nonunitary plaquette; refine the mesh')
         END IF
         max_phase = MAX(max_phase, local_phase)
         y = contribution - compensation
         t = value + y
         compensation = (t - value) - y
         value = t
      END DO
      CALL close_file(store)
      report%value = value
      WRITE (io_unit, '(A,I0)') 'PARAMETER_DIMENSION ', dim
      WRITE (io_unit, '(A,*(I0,1X))') 'MESH ', shape
      WRITE (io_unit, '(A,ES26.17)') 'CHERN_RAW ', value
      WRITE (io_unit, '(A,ES26.17)') 'MINIMUM_LINK_SINGULAR_VALUE ', minimum
      WRITE (io_unit, '(A,ES26.17)') 'MAXIMUM_PLAQUETTE_PHASE ', max_phase
      IF (physical) THEN
         WRITE (io_unit, '(A,ES26.17)') 'MAXIMUM_METRIC_ERROR ', max_metric
         WRITE (io_unit, '(A,ES26.17)') 'MAXIMUM_SEAM_ERROR ', seam_error
         WRITE (io_unit, '(A,ES26.17)') 'MINIMUM_DIRECT_GAP_HA ', min_gap
         IF (common_reference .AND. prefix) THEN
            WRITE (io_unit, '(A,ES26.17)') 'SAMPLED_INDIRECT_GAP_HA ', conduction_min - valence_max
         ELSE
            WRITE (io_unit, '(A)') 'INDIRECT_GAP_NOT_EVALUATED: no asserted common energy reference or not a lowest-band prefix'
         END IF
      END IF
      WRITE (io_unit, '(A)') 'MESH_CONVERGENCE_NOT_ESTABLISHED: repeat on finer meshes; retain unrounded C2 values.'
   CONTAINS
! **************************************************************************************************
!> \brief Read one frame and update validation diagnostics on its first use.
!> \param number One-based frame index in the extended mesh
!> \param state Physical frame with selected coefficients
! **************************************************************************************************
      SUBROUTINE frame(number, state)
      INTEGER, INTENT(IN)                                :: number
      TYPE(snapshot_type), INTENT(OUT)                   :: state

      INTEGER                                            :: iband

         CALL read_snapshot(TRIM(files(number)), points(number), state)
         IF (state%rank /= rank .OR. SIZE(state%atoms) /= nat) THEN
            CPABORT('Frame rank/atom count differs from manifest')
         END IF
         IF (.NOT. checked(number)) THEN
            CALL check_snapshot(state, tol_metric, tol_gap, err, gap)
            max_metric = MAX(max_metric, err)
            min_gap = MIN(min_gap, gap)
            checked(number) = .TRUE.
            DO iband = 1, rank
               IF (state%bands(iband) /= iband) prefix = .FALSE.
            END DO
            IF (rank >= SIZE(state%energies)) prefix = .FALSE.
            IF (prefix) THEN
               energy_shift = offsets(number)
               valence_max = MAX(valence_max, state%energies(rank) + energy_shift)
               conduction_min = MIN(conduction_min, state%energies(rank + 1) + energy_shift)
            END IF
         END IF
      END SUBROUTINE frame
! **************************************************************************************************
!> \brief Read a previously validated link from the scratch stream.
!> \param vertex_index One-based mesh vertex
!> \param direction One-based parameter axis
!> \param matrix Unitary selected-state link
! **************************************************************************************************
      SUBROUTINE read_link(vertex_index, direction, matrix)
      INTEGER, INTENT(IN)                                :: vertex_index, direction
      COMPLEX(KIND=dp), INTENT(OUT)                      :: matrix(:, :)

      INTEGER(KIND=int64)                                :: position

         position = (INT(vertex_index - 1, int64)*INT(dim, int64) + INT(direction - 1, int64))*INT(record_size, int64) + 1_int64
         READ (store, POS=position, IOSTAT=stat) matrix
         CPASSERT(stat == 0)
      END SUBROUTINE read_link
   END SUBROUTINE analyze_mesh

! **************************************************************************************************
!> \brief Compare two jointly refined meshes of the same declared family; never round the estimate.
!> \param coarse Coarser mesh report
!> \param fine Finer mesh report
!> \param tolerance Maximum inter-mesh difference and distance to an integer
!> \param io_unit Unit receiving convergence diagnostics
! **************************************************************************************************
   SUBROUTINE compare_refinements(coarse, fine, tolerance, io_unit)
      TYPE(phason_result_type), INTENT(IN) :: coarse, fine
      REAL(KIND=dp), INTENT(IN) :: tolerance
      INTEGER, INTENT(IN) :: io_unit
      TYPE(snapshot_type), ALLOCATABLE :: a, b
      COMPLEX(KIND=dp), ALLOCATABLE :: m(:, :), metric(:, :)
      INTEGER :: dim, i
      REAL(KIND=dp) :: difference, residual
      IF (.NOT. ieee_is_finite(tolerance)) THEN
         CPABORT('Nonfinite convergence tolerance')
      END IF
      IF (tolerance <= 0.0_dp .OR. tolerance >= 0.5_dp) THEN
         CPABORT('Invalid convergence tolerance')
      END IF
      dim = coarse%dimension
      IF (fine%dimension /= dim .OR. fine%rank /= coarse%rank .OR. fine%nat /= coarse%nat) THEN
         CPABORT('Incompatible refinement reports')
      END IF
      DO i = 1, dim
         IF (fine%shape(i) < 2*coarse%shape(i) .OR. MOD(fine%shape(i), coarse%shape(i)) /= 0) THEN
            CPABORT('Refine every mesh axis by an integer factor of at least two')
         END IF
      END DO
      IF (MAXVAL(ABS(coarse%origin - fine%origin)) > 1.e-10_dp .OR. &
          MAXVAL(ABS(coarse%directions - fine%directions)) > 1.e-10_dp) THEN
         CPABORT('Refinement changes the parameter origin or reciprocal periods')
      END IF
      IF (coarse%nat > 0) THEN
         ALLOCATE (a, b)
         CALL read_snapshot(TRIM(coarse%first_file), coarse%point, a)
         CALL read_snapshot(TRIM(fine%first_file), fine%point, b)
         ALLOCATE (m(a%rank, a%rank), metric(a%rank, a%rank))
         CALL snapshot_overlap(a, b, a%k, a%k, m)
         metric(:, :) = MATMUL(CONJG(TRANSPOSE(m)), m)
         DO i = 1, a%rank
            metric(i, i) = metric(i, i) - 1.0_dp
         END DO
         IF (MAXVAL(ABS(metric)) > 1.e-7_dp) THEN
            CPABORT('Refinement changes the initial physical subspace')
         END IF
      END IF
      difference = ABS(fine%value - coarse%value)
      residual = ABS(fine%value - ANINT(fine%value))
      WRITE (io_unit, '(A,ES26.17)') 'REFINEMENT_CHANGE ', difference
      WRITE (io_unit, '(A,ES26.17)') 'INTEGER_RESIDUAL ', residual
      IF (difference > tolerance .OR. residual > tolerance) THEN
         CPABORT('Chern mesh refinement has not converged to the requested tolerance')
      END IF
      WRITE (io_unit, '(A)') 'REFINEMENT_CHECK_PASSED: numerical check, not a certificate of a bulk gap.'
   END SUBROUTINE compare_refinements
END MODULE topology_phason_analysis

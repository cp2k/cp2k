!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Reader and physical moving-basis links for version-1 topology snapshots.
! **************************************************************************************************
MODULE topology_snapshot
   USE ai_moments,                      ONLY: cossin
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_snapshot'
   REAL(KIND=dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp, cell_tol = 1.e-10_dp
   TYPE :: snapshot_shell_type
      INTEGER :: first = 0, count = 0, lmin = 0, lmax = 0
      REAL(KIND=dp) :: center(3) = 0.0_dp, radius = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE :: exponent(:), radii(:), contraction(:, :)
   END TYPE snapshot_shell_type
   TYPE :: snapshot_atom_type
      INTEGER :: kind = 0, count = 0
      REAL(KIND=dp) :: center(3) = 0.0_dp
      TYPE(snapshot_shell_type), ALLOCATABLE :: shells(:)
   END TYPE snapshot_atom_type
   TYPE :: snapshot_type
      INTEGER :: nao = 0, rank = 0, nspin = 0, channel = 0
      INTEGER :: periodic(3) = 0
      REAL(KIND=dp) :: cell(3, 3) = 0.0_dp, inverse(3, 3) = 0.0_dp, k(3) = 0.0_dp
      TYPE(snapshot_atom_type), ALLOCATABLE :: atoms(:)
      INTEGER, ALLOCATABLE :: bands(:)
      REAL(KIND=dp), ALLOCATABLE :: energies(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: coefficients(:, :)
   END TYPE snapshot_type
   PUBLIC :: snapshot_type, read_snapshot, snapshot_overlap, check_snapshot, check_snapshot_seam
CONTAINS

! **************************************************************************************************
!> \brief Read a single point, retaining only its states; reject malformed/truncated exports.
!> \param filename Version-1 STATE_EXPORT file
!> \param point One-based point index
!> \param state Selected frame
! **************************************************************************************************
   SUBROUTINE read_snapshot(filename, point, state)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(IN)                                :: point
      TYPE(snapshot_type), INTENT(OUT)                   :: state

      CHARACTER(LEN=256)                                 :: header
      INTEGER                                            :: first, i, ia, ic, ios, ip, iset, ix, j, &
                                                            nat, nc, nk, np, ns, nt, offset, p, &
                                                            powers(3), refpowers(3), unit
      LOGICAL, ALLOCATABLE                               :: covered(:)
      REAL(KIND=dp)                                      :: det, k(3), pair(2)
      REAL(KIND=dp), ALLOCATABLE                         :: energies(:), row(:)

      CALL open_file(filename, unit_number=unit, file_status='OLD', file_action='READ')
      READ (unit, '(A)', IOSTAT=ios) header
      CPASSERT(ios == 0)
      IF (TRIM(header) /= 'CP2K_TOPOLOGY_STATE 1') THEN
         CPABORT('Unsupported topology snapshot version')
      END IF
      READ (unit, *, IOSTAT=ios) nat, state%nao, state%rank, nk, state%nspin, nt, state%channel
      CPASSERT(ios == 0)
      IF (MIN(nat, state%nao, state%rank, nk, nt, state%channel) < 1) THEN
         CPABORT('Invalid snapshot dimensions')
      END IF
      IF (state%rank > nt .OR. point < 1 .OR. point > nk) THEN
         CPABORT('Invalid snapshot point or rank')
      END IF
      IF (state%nspin /= 1 .AND. state%nspin /= 2) THEN
         CPABORT('Invalid snapshot spin components')
      END IF
      ALLOCATE (state%bands(state%rank), state%energies(nt), energies(nt), &
                state%coefficients(state%nao*state%nspin, state%rank), state%atoms(nat))
      READ (unit, *, IOSTAT=ios) state%bands
      CPASSERT(ios == 0)
      IF (ANY(state%bands < 1) .OR. ANY(state%bands > nt)) THEN
         CPABORT('Invalid snapshot band indices')
      END IF
      DO i = 2, state%rank
         IF (state%bands(i) <= state%bands(i - 1)) THEN
            CPABORT('Snapshot bands must be strictly ordered')
         END IF
      END DO
      DO j = 1, 3
         READ (unit, *, IOSTAT=ios) state%cell(:, j)
         CPASSERT(ios == 0)
      END DO
      IF (.NOT. ALL(ieee_is_finite(state%cell))) THEN
         CPABORT('Nonfinite snapshot cell')
      END IF
      ASSOCIATE (a => state%cell, b => state%inverse)
         b(1, :) = [a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2), a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3), a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)]
         b(2, :) = [a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3), a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1), a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)]
         b(3, :) = [a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1), a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2), a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)]
         det = DOT_PRODUCT(a(:, 1), b(1, :))
         IF (ABS(det) < cell_tol) THEN
            CPABORT('Singular snapshot cell')
         END IF
         b(:, :) = b/det
      END ASSOCIATE
      READ (unit, *, IOSTAT=ios) state%periodic
      CPASSERT(ios == 0)
      IF (ANY(state%periodic < 0) .OR. ANY(state%periodic > 1)) THEN
         CPABORT('Invalid snapshot periodicity')
      END IF
      offset = 0
      DO ia = 1, nat
         ASSOCIATE (atom => state%atoms(ia))
            READ (unit, *, IOSTAT=ios) ix, atom%kind, ns, atom%count
            CPASSERT(ios == 0)
            IF (ix /= ia .OR. MIN(atom%kind, ns, atom%count) < 1) THEN
               CPABORT('Invalid snapshot atom')
            END IF
            READ (unit, *, IOSTAT=ios) atom%center
            CPASSERT(ios == 0)
            IF (.NOT. ALL(ieee_is_finite(atom%center))) THEN
               CPABORT('Nonfinite atom center')
            END IF
            ALLOCATE (atom%shells(ns), covered(atom%count))
            covered(:) = .FALSE.
            DO iset = 1, ns
               ASSOCIATE (s => atom%shells(iset))
                  READ (unit, *, IOSTAT=ios) first, s%count, np, nc, s%lmin, s%radius
                  CPASSERT(ios == 0)
                  IF (MIN(first, s%count, np, nc) < 1 .OR. s%lmin < 0) THEN
                     CPABORT('Invalid Gaussian shell')
                  END IF
                  IF (.NOT. ieee_is_finite(s%radius) .OR. s%radius <= 0) THEN
                     CPABORT('Invalid shell screening radius')
                  END IF
                  IF (first + s%count - 1 > atom%count) THEN
                     CPABORT('Shell AO range exceeds atom')
                  END IF
                  IF (ANY(covered(first:first + s%count - 1))) THEN
                     CPABORT('Overlapping shell AO ranges')
                  END IF
                  covered(first:first + s%count - 1) = .TRUE.
                  s%lmax = 0
                  DO WHILE ((s%lmax + 1)*(s%lmax + 2)*(s%lmax + 3)/6 < nc)
                     s%lmax = s%lmax + 1
                  END DO
                  IF ((s%lmax + 1)*(s%lmax + 2)*(s%lmax + 3)/6 /= nc .OR. s%lmin > s%lmax) THEN
                     CPABORT('Invalid Cartesian Gaussian count')
                  END IF
                  CALL init_orbital_pointers(s%lmax + 1)
                  s%first = offset + first
                  s%center(:) = atom%center
                  ALLOCATE (s%exponent(np), s%radii(np), s%contraction(nc*np, s%count), row(s%count))
                  DO ip = 1, np
                     READ (unit, *, IOSTAT=ios) s%exponent(ip), s%radii(ip)
                     CPASSERT(ios == 0)
                     IF (.NOT. ieee_is_finite(s%exponent(ip)) .OR. .NOT. ieee_is_finite(s%radii(ip))) THEN
                        CPABORT('Nonfinite primitive Gaussian')
                     END IF
                     IF (s%exponent(ip) <= 0 .OR. s%radii(ip) <= 0) THEN
                        CPABORT('Invalid primitive Gaussian')
                     END IF
                     DO ic = 1, nc
                        READ (unit, *, IOSTAT=ios) powers, row
                        CPASSERT(ios == 0)
                        refpowers(:) = indco(:, ic)
                        IF (ANY(powers /= refpowers)) THEN
                           CPABORT('Unexpected Cartesian Gaussian ordering')
                        END IF
                        IF (.NOT. ALL(ieee_is_finite(row))) THEN
                           CPABORT('Nonfinite contraction')
                        END IF
                        s%contraction((ip - 1)*nc + ic, :) = row
                     END DO
                  END DO
                  DEALLOCATE (row)
               END ASSOCIATE
            END DO
            IF (.NOT. ALL(covered)) THEN
               CPABORT('Incomplete atom AO coverage')
            END IF
            DEALLOCATE (covered)
            offset = offset + atom%count
         END ASSOCIATE
      END DO
      IF (offset /= state%nao) THEN
         CPABORT('Snapshot AO count mismatch')
      END IF
      DO p = 1, nk
         READ (unit, *, IOSTAT=ios) ix, k
         CPASSERT(ios == 0)
         IF (ix /= p .OR. .NOT. ALL(ieee_is_finite(k))) THEN
            CPABORT('Invalid snapshot k-point')
         END IF
         READ (unit, *, IOSTAT=ios) energies
         CPASSERT(ios == 0)
         IF (.NOT. ALL(ieee_is_finite(energies))) THEN
            CPABORT('Nonfinite snapshot spectrum')
         END IF
         IF (ANY(energies(2:) < energies(:nt - 1))) THEN
            CPABORT('Unordered snapshot spectrum')
         END IF
         IF (p == point) THEN
            state%k(:) = k
            state%energies(:) = energies
         END IF
         DO j = 1, state%rank
            DO i = 1, state%nao*state%nspin
               READ (unit, *, IOSTAT=ios) pair
               CPASSERT(ios == 0)
               IF (.NOT. ALL(ieee_is_finite(pair))) THEN
                  CPABORT('Nonfinite state coefficient')
               END IF
               IF (p == point) state%coefficients(i, j) = CMPLX(pair(1), pair(2), dp)
            END DO
         END DO
      END DO
      DO
         READ (unit, '(A)', IOSTAT=ios) header
         IF (ios < 0) EXIT
         IF (ios /= 0 .OR. LEN_TRIM(header) /= 0) THEN
            CPABORT('Unexpected trailing snapshot content')
         END IF
      END DO
      CALL close_file(unit)
   END SUBROUTINE read_snapshot

! **************************************************************************************************
!> \brief Apply screened Gaussian cross-geometry operator by atom blocks, without a dense AO matrix.
!> \param a Left physical state
!> \param b Right physical state
!> \param ka Unwrapped left fractional k-point
!> \param kb Unwrapped right fractional k-point
!> \param overlap Selected-state overlap including all spinor components
! **************************************************************************************************
   SUBROUTINE snapshot_overlap(a, b, ka, kb, overlap)
      TYPE(snapshot_type), INTENT(IN)                    :: a, b
      REAL(KIND=dp), INTENT(IN)                          :: ka(3), kb(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: overlap(:, :)

      COMPLEX(KIND=dp) :: part(SIZE(overlap, 1), SIZE(overlap, 2))
      INTEGER                                            :: ia

      IF (ANY(SHAPE(overlap) /= [a%rank, b%rank])) THEN
         CPABORT('Incorrect overlap output dimensions')
      END IF
      IF (a%nspin /= b%nspin .OR. a%channel /= b%channel .OR. a%rank /= b%rank) THEN
         CPABORT('Changed spin or rank')
      END IF
      IF (ANY(a%bands /= b%bands)) THEN
         CPABORT('Changed selected band indices')
      END IF
      IF (SIZE(a%atoms) /= SIZE(b%atoms) .OR. a%nao /= b%nao) THEN
         CPABORT('Changed atom or AO count')
      END IF
      IF (MAXVAL(ABS(a%cell - b%cell)) > cell_tol .OR. ANY(a%periodic /= b%periodic)) THEN
         CPABORT('Changed cell')
      END IF
      IF (.NOT. ALL(ieee_is_finite(ka)) .OR. .NOT. ALL(ieee_is_finite(kb))) THEN
         CPABORT('Nonfinite k-point')
      END IF
      IF (MAXVAL(ABS(ka - a%k - ANINT(ka - a%k))) > cell_tol .OR. &
          MAXVAL(ABS(kb - b%k - ANINT(kb - b%k))) > cell_tol) THEN
         CPABORT('Frame k-point does not match request')
      END IF
      overlap(:, :) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(A,b,ka,kb) PRIVATE(ia,part) REDUCTION(+:overlap) SCHEDULE(DYNAMIC)
      DO ia = 1, SIZE(a%atoms)
         CALL atom_overlap(a, b, ia, ka, kb, part)
         overlap(:, :) = overlap + part
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE snapshot_overlap

! **************************************************************************************************
!> \brief Screened atom-row contribution using CP2K's cossin primitive integral recurrence.
!> \param a Left physical state
!> \param b Right physical state
!> \param ia Left atom whose AO rows are contracted
!> \param ka Unwrapped left fractional k-point
!> \param kb Unwrapped right fractional k-point
!> \param atom_link This atom's contribution to the selected-state overlap
! **************************************************************************************************
   SUBROUTINE atom_overlap(a, b, ia, ka, kb, atom_link)
      TYPE(snapshot_type), INTENT(IN)                    :: a, b
      INTEGER, INTENT(IN)                                :: ia
      REAL(KIND=dp), INTENT(IN)                          :: ka(3), kb(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: atom_link(:, :)

      COMPLEX(KIND=dp)                                   :: phase
      COMPLEX(KIND=dp), ALLOCATABLE                      :: block(:, :), oc(:, :, :)
      INTEGER                                            :: af, bf, hi(3), i, ib, image(3), j, k, &
                                                            lo(3), nc_a, nc_b, s, sa, sb
      REAL(KIND=dp)                                      :: bounds(3), cutoff, disp(3), q(3), rb(3)
      REAL(KIND=dp), ALLOCATABLE                         :: cosine(:, :), sine(:, :)

      q(:) = 2.0_dp*pi*MATMUL(TRANSPOSE(a%inverse), kb - ka)
      atom_link(:, :) = 0.0_dp
      DO sa = 1, SIZE(a%atoms(ia)%shells)
         ASSOCIATE (left => a%atoms(ia)%shells(sa))
            ALLOCATE (oc(left%count, b%rank, a%nspin))
            oc(:, :, :) = 0.0_dp
            nc_a = SIZE(left%contraction, 1)
            DO ib = 1, SIZE(b%atoms)
               DO sb = 1, SIZE(b%atoms(ib)%shells)
                  ASSOCIATE (right => b%atoms(ib)%shells(sb))
                     cutoff = left%radius + right%radius
                     disp(:) = MATMUL(a%inverse, right%center - left%center)
                     DO i = 1, 3
                        bounds(i) = cutoff*NORM2(a%inverse(i, :))
                     END DO
                     lo(:) = CEILING(-disp - bounds)
                     hi(:) = FLOOR(-disp + bounds)
                     WHERE (a%periodic == 0)
                        lo = 0
                        hi = 0
                     END WHERE
                     nc_b = SIZE(right%contraction, 1)
                     ALLOCATE (cosine(nc_a, nc_b), sine(nc_a, nc_b), block(left%count, right%count))
                     block(:, :) = 0.0_dp
                     DO k = lo(3), hi(3)
                        DO j = lo(2), hi(2)
                           DO i = lo(1), hi(1)
                              image(:) = [i, j, k]
                              rb(:) = right%center + MATMUL(a%cell, REAL(image, dp))
                              IF (NORM2(rb - left%center) > cutoff) CYCLE
                              CALL cossin(left%lmax, SIZE(left%exponent), left%exponent, left%radii, left%lmin, &
                                          right%lmax, SIZE(right%exponent), right%exponent, right%radii, right%lmin, &
                                          left%center, rb, q, cosine, sine)
                              phase = EXP(CMPLX(0.0_dp, 2.0_dp*pi*DOT_PRODUCT(kb, REAL(image, dp)), dp))
                              block(:, :) = block + phase*MATMUL(TRANSPOSE(left%contraction), &
                                                                 MATMUL(CMPLX(cosine, -sine, dp), right%contraction))
                           END DO
                        END DO
                     END DO
                     DO s = 1, a%nspin
                        bf = right%first + (s - 1)*b%nao
                        oc(:, :, s) = oc(:, :, s) + MATMUL(block, b%coefficients(bf:bf + right%count - 1, :))
                     END DO
                     DEALLOCATE (cosine, sine, block)
                  END ASSOCIATE
               END DO
            END DO
            DO s = 1, a%nspin
               af = left%first + (s - 1)*a%nao
               atom_link(:, :) = atom_link + MATMUL(CONJG(TRANSPOSE(a%coefficients(af:af + left%count - 1, :))), oc(:, :, s))
            END DO
            DEALLOCATE (oc)
         END ASSOCIATE
      END DO
   END SUBROUTINE atom_overlap

! **************************************************************************************************
!> \brief Check AO-metric normalization and separation at every selected/excluded boundary.
!> \param a Physical frame to validate
!> \param metric_tol Maximum entrywise error in C^dagger S C minus identity
!> \param gap_tol Minimum sampled selected/excluded-band separation in hartree
!> \param metric_error Observed normalization error
!> \param gap Smallest sampled selected/excluded-band separation in hartree
! **************************************************************************************************
   SUBROUTINE check_snapshot(a, metric_tol, gap_tol, metric_error, gap)
      TYPE(snapshot_type), INTENT(IN)                    :: a
      REAL(KIND=dp), INTENT(IN)                          :: metric_tol, gap_tol
      REAL(KIND=dp), INTENT(OUT)                         :: metric_error, gap

      INTEGER                                            :: i
      LOGICAL                                            :: selected(SIZE(a%energies))
      COMPLEX(KIND=dp)                                   :: metric(a%rank, a%rank)

      CALL snapshot_overlap(a, a, a%k, a%k, metric)
      DO i = 1, a%rank
         metric(i, i) = metric(i, i) - 1.0_dp
      END DO
      metric_error = MAXVAL(ABS(metric))
      IF (.NOT. ieee_is_finite(metric_error) .OR. metric_error > metric_tol) THEN
         CPABORT('Frame is not AO-metric normalized')
      END IF
      selected(:) = .FALSE.
      selected(a%bands) = .TRUE.
      gap = HUGE(1.0_dp)
      DO i = 1, SIZE(a%energies) - 1
         IF (selected(i) .NEQV. selected(i + 1)) gap = MIN(gap, a%energies(i + 1) - a%energies(i))
      END DO
      IF (gap == HUGE(1.0_dp) .OR. gap <= gap_tol) THEN
         CPABORT('No verified sampled subspace isolation')
      END IF
   END SUBROUTINE check_snapshot

! **************************************************************************************************
!> \brief Check explicitly prescribed atom permutation, periodic translations and basis at a seam.
!> \param a Start frame
!> \param b Endpoint frame
!> \param permutation One-based endpoint-atom to start-atom map
!> \param tolerance Absolute tolerance for geometry and basis agreement
! **************************************************************************************************
   SUBROUTINE check_snapshot_seam(a, b, permutation, tolerance)
      TYPE(snapshot_type), INTENT(IN)                    :: a, b
      INTEGER, INTENT(IN)                                :: permutation(:)
      REAL(KIND=dp), INTENT(IN)                          :: tolerance

      INTEGER                                            :: i, j, n, s
      REAL(KIND=dp)                                      :: delta(3), shift(3)

      n = SIZE(a%atoms)
      IF (SIZE(b%atoms) /= n .OR. SIZE(permutation) /= n) THEN
         CPABORT('Seam atom count mismatch')
      END IF
      IF (ANY(permutation < 1) .OR. ANY(permutation > n)) THEN
         CPABORT('Invalid seam permutation')
      END IF
      DO i = 1, n
         IF (COUNT(permutation == i) /= 1) THEN
            CPABORT('Seam atom map is not a permutation')
         END IF
      END DO
      IF (MAXVAL(ABS(a%cell - b%cell)) > tolerance .OR. ANY(a%periodic /= b%periodic)) THEN
         CPABORT('Seam cell mismatch')
      END IF
      DO i = 1, n
         j = permutation(i)
         ASSOCIATE (left => a%atoms(j), right => b%atoms(i))
            IF (left%kind /= right%kind .OR. left%count /= right%count) THEN
               CPABORT('Seam changes atom kind or AO count')
            END IF
            IF (SIZE(left%shells) /= SIZE(right%shells)) THEN
               CPABORT('Seam changes basis')
            END IF
            delta(:) = right%center - left%center
            shift(:) = ANINT(MATMUL(a%inverse, delta))
            WHERE (a%periodic == 0) shift = 0.0_dp
            IF (MAXVAL(ABS(delta - MATMUL(a%cell, shift))) > tolerance) THEN
               CPABORT('Geometry does not close at seam')
            END IF
            DO s = 1, SIZE(left%shells)
               ASSOCIATE (x => left%shells(s), y => right%shells(s))
                  IF (x%count /= y%count .OR. x%lmin /= y%lmin .OR. x%lmax /= y%lmax) THEN
                     CPABORT('Seam changes shell')
                  END IF
                  IF (ANY(SHAPE(x%contraction) /= SHAPE(y%contraction))) THEN
                     CPABORT('Seam changes contractions')
                  END IF
                  IF (SIZE(x%exponent) /= SIZE(y%exponent)) THEN
                     CPABORT('Seam changes primitives')
                  END IF
                  IF (MAXVAL(ABS(x%radii - y%radii)) > tolerance .OR. ABS(x%radius - y%radius) > tolerance) THEN
                     CPABORT('Seam changes Gaussian screening radii')
                  END IF
                  IF (MAXVAL(ABS(x%exponent - y%exponent)) > tolerance .OR. &
                      MAXVAL(ABS(x%contraction - y%contraction)) > tolerance) THEN
                     CPABORT('Seam changes physical basis')
                  END IF
               END ASSOCIATE
            END DO
         END ASSOCIATE
      END DO
   END SUBROUTINE check_snapshot_seam
END MODULE topology_snapshot

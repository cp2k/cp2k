!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Bounded geometry-only cache for native atom-grid interpolation.
MODULE qs_native_grid_cache
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE pw_grid_types,                   ONLY: pw_grid_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: native_grid_interp_offset_min = -6, native_grid_interp_offset_max = 6
   INTEGER, PARAMETER :: native_grid_interp_npts = native_grid_interp_offset_max - native_grid_interp_offset_min + 1
   INTEGER(KIND=int_8), PARAMETER :: cache_budget = 128_int_8*1024*1024

   TYPE native_grid_interpolation_stencil_type
      INTEGER :: relative_index(native_grid_interp_npts, 3) = 0
      REAL(dp) :: weight(native_grid_interp_npts, 3) = 0.0_dp
      LOGICAL :: valid(native_grid_interp_npts, 3) = .FALSE.
      LOGICAL :: active = .FALSE.
   END TYPE native_grid_interpolation_stencil_type

   TYPE stencil_entry_type
      TYPE(native_grid_interpolation_stencil_type) :: stencil
      REAL(dp) :: point(3) = 0.0_dp
      LOGICAL :: ready = .FALSE.
   END TYPE stencil_entry_type

   TYPE native_grid_cache_type
      PRIVATE
      TYPE(stencil_entry_type), ALLOCATABLE :: entries(:)
      REAL(dp) :: dh_inv(3, 3) = 0.0_dp, hmat(3, 3) = 0.0_dp
      INTEGER :: npts(3) = 0, periodic(3) = 0, nrows = -1
      LOGICAL :: wrap = .FALSE., prepared = .FALSE.
   END TYPE native_grid_cache_type

   PUBLIC :: native_grid_cache_type, native_grid_interpolation_stencil_type
   PUBLIC :: native_grid_interp_offset_min, native_grid_interp_offset_max, native_grid_interp_npts
   PUBLIC :: native_grid_interpolation_weights
   PUBLIC :: prepare_native_grid_cache, release_native_grid_cache, fetch_native_grid_stencil, store_native_grid_stencil
CONTAINS
! **************************************************************************************************
!> \brief C2 interpolation weights from two adjacent degree-11 Lagrange polynomials.
!>        A quintic smoothstep joins their values, slopes and curvatures at half-grid boundaries.
!>        The 13-node union retains polynomial reproduction, exact grid-node values and reflection.
!> \param fraction Coordinate relative to the nearest grid node, in [-1/2,1/2].
!> \param weights Value weights on nodes -6 through 6.
!> \param derivative_weights Optional first derivatives with respect to fraction.
!> \param second_derivative_weights Optional second derivatives with respect to fraction.
! **************************************************************************************************
   SUBROUTINE native_grid_interpolation_weights(fraction, weights, derivative_weights, second_derivative_weights)
      REAL(dp), INTENT(IN)                               :: fraction
      REAL(dp), INTENT(OUT)                              :: weights(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: derivative_weights(:), &
                                                            second_derivative_weights(:)

      INTEGER, PARAMETER :: n = native_grid_interp_npts - 1

      LOGICAL                                            :: first, second
      REAL(dp)                                           :: blend, blend_d, blend_dd, left(n), &
                                                            left_d(n), left_dd(n), right(n), &
                                                            right_d(n), right_dd(n), t

      CPASSERT(SIZE(weights) == native_grid_interp_npts)
      IF (PRESENT(derivative_weights)) THEN
         CPASSERT(SIZE(derivative_weights) == native_grid_interp_npts)
      END IF
      IF (PRESENT(second_derivative_weights)) THEN
         CPASSERT(SIZE(second_derivative_weights) == native_grid_interp_npts)
      END IF
      second = PRESENT(second_derivative_weights)
      first = PRESENT(derivative_weights) .OR. second
      ! Nearest-node selection can round a coordinate by one ulp across a half-grid
      ! boundary. The C2 polynomial continuation also handles that limiting case.
      CALL lagrange_polynomial_weights(fraction, native_grid_interp_offset_min, left, left_d, left_dd, first, second)
      CALL lagrange_polynomial_weights(fraction, native_grid_interp_offset_min + 1, right, right_d, right_dd, first, second)
      t = fraction + 0.5_dp
      blend = t**3*(10.0_dp + t*(-15.0_dp + 6.0_dp*t))
      blend_d = 30.0_dp*t**2*(1.0_dp - t)**2
      blend_dd = 60.0_dp*t*(1.0_dp - t)*(1.0_dp - 2.0_dp*t)
      weights = 0.0_dp
      weights(:n) = (1.0_dp - blend)*left
      weights(2:) = weights(2:) + blend*right
      IF (PRESENT(derivative_weights)) THEN
         derivative_weights = 0.0_dp
         derivative_weights(:n) = (1.0_dp - blend)*left_d - blend_d*left
         derivative_weights(2:) = derivative_weights(2:) + blend*right_d + blend_d*right
      END IF
      IF (second) THEN
         second_derivative_weights = 0.0_dp
         second_derivative_weights(:n) = (1.0_dp - blend)*left_dd - 2.0_dp*blend_d*left_d - blend_dd*left
         second_derivative_weights(2:) = second_derivative_weights(2:) + &
                                         blend*right_dd + 2.0_dp*blend_d*right_d + blend_dd*right
      END IF
   END SUBROUTINE native_grid_interpolation_weights

! **************************************************************************************************
!> \brief Lagrange polynomial weights and derivatives without division by the evaluation distance.
!> \param fraction Evaluation coordinate.
!> \param first_node First of the contiguous integer nodes.
!> \param weights Value weights.
!> \param derivative_weights First derivative weights, zero when not requested.
!> \param second_derivative_weights Second derivative weights, zero when not requested.
!> \param first Calculate first derivatives.
!> \param second Calculate second derivatives.
! **************************************************************************************************
   SUBROUTINE lagrange_polynomial_weights(fraction, first_node, weights, derivative_weights, &
                                          second_derivative_weights, first, second)
      REAL(dp), INTENT(IN)                               :: fraction
      INTEGER, INTENT(IN)                                :: first_node
      REAL(dp), INTENT(OUT)                              :: weights(:), derivative_weights(:), &
                                                            second_derivative_weights(:)
      LOGICAL, INTENT(IN)                                :: first, second

      INTEGER                                            :: inode, jnode, knode, lnode, n
      REAL(dp)                                           :: denominator, numerator, &
                                                            nodes(SIZE(weights))

      n = SIZE(weights)
      derivative_weights = 0.0_dp
      second_derivative_weights = 0.0_dp
      DO inode = 1, n
         nodes(inode) = REAL(first_node + inode - 1, dp)
      END DO
      DO inode = 1, n
         denominator = 1.0_dp
         numerator = 1.0_dp
         DO jnode = 1, n
            IF (jnode == inode) CYCLE
            denominator = denominator*(nodes(inode) - nodes(jnode))
            numerator = numerator*(fraction - nodes(jnode))
         END DO
         weights(inode) = numerator/denominator
         IF (.NOT. first) CYCLE
         DO knode = 1, n
            IF (knode == inode) CYCLE
            numerator = 1.0_dp
            DO jnode = 1, n
               IF (jnode == inode .OR. jnode == knode) CYCLE
               numerator = numerator*(fraction - nodes(jnode))
            END DO
            derivative_weights(inode) = derivative_weights(inode) + numerator/denominator
            IF (.NOT. second) CYCLE
            DO lnode = knode + 1, n
               IF (lnode == inode) CYCLE
               numerator = 1.0_dp
               DO jnode = 1, n
                  IF (jnode == inode .OR. jnode == knode .OR. jnode == lnode) CYCLE
                  numerator = numerator*(fraction - nodes(jnode))
               END DO
               second_derivative_weights(inode) = second_derivative_weights(inode) + 2.0_dp*numerator/denominator
            END DO
         END DO
      END DO
   END SUBROUTINE lagrange_polynomial_weights

! **************************************************************************************************
!> \brief Prepare outside parallel regions. Cache a bounded prefix and compute the rest normally.
!> \param cache per-QS-environment cache
!> \param grid current auxiliary grid
!> \param cell current cell and periodicity
!> \param nrows local atom-grid row count
!> \param wrap wrap the auxiliary cell in nonperiodic directions
!> \param max_bytes optional memory budget for testing, in bytes
! **************************************************************************************************
   SUBROUTINE prepare_native_grid_cache(cache, grid, cell, nrows, wrap, max_bytes)
      TYPE(native_grid_cache_type), INTENT(INOUT)        :: cache
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      TYPE(cell_type), INTENT(IN)                        :: cell
      INTEGER, INTENT(IN)                                :: nrows
      LOGICAL, INTENT(IN)                                :: wrap
      INTEGER(KIND=int_8), INTENT(IN), OPTIONAL          :: max_bytes

      INTEGER                                            :: count, ierr
      INTEGER(KIND=int_8)                                :: budget, bytes_per_entry
      LOGICAL                                            :: unchanged
      TYPE(stencil_entry_type)                           :: entry

      budget = cache_budget
      IF (PRESENT(max_bytes)) budget = MAX(0_int_8, max_bytes)
      bytes_per_entry = INT((STORAGE_SIZE(ENTRY) + 7)/8, int_8)
      count = INT(MIN(INT(MAX(0, nrows), int_8), budget/bytes_per_entry))
      unchanged = .FALSE.
      IF (ALLOCATED(cache%entries)) THEN
         unchanged = SIZE(cache%entries) == count
         IF (.NOT. unchanged) CALL release_native_grid_cache(cache)
      END IF
      IF (.NOT. ALLOCATED(cache%entries)) THEN
         ALLOCATE (cache%entries(count), STAT=ierr)
         IF (ierr /= 0) RETURN
      END IF
      unchanged = unchanged .AND. cache%prepared .AND. cache%nrows == nrows
      unchanged = unchanged .AND. ALL(cache%dh_inv == grid%dh_inv) .AND. ALL(cache%npts == grid%npts)
      unchanged = unchanged .AND. ALL(cache%hmat == cell%hmat) .AND. ALL(cache%periodic == cell%perd)
      unchanged = unchanged .AND. (cache%wrap .EQV. wrap)
      IF (.NOT. unchanged) cache%entries%ready = .FALSE.
      cache%dh_inv = grid%dh_inv
      cache%npts = grid%npts
      cache%hmat = cell%hmat
      cache%periodic = cell%perd
      cache%nrows = nrows
      cache%wrap = wrap
      cache%prepared = .TRUE.
   END SUBROUTINE prepare_native_grid_cache

! **************************************************************************************************
!> \brief Release geometry storage, including when the QS environment is only partly released.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE release_native_grid_cache(cache)
      TYPE(native_grid_cache_type), INTENT(INOUT)        :: cache

      IF (ALLOCATED(cache%entries)) DEALLOCATE (cache%entries)
      cache%prepared = .FALSE.
   END SUBROUTINE release_native_grid_cache

! **************************************************************************************************
!> \brief Read a stencil only on an exact coordinate match. Safe for concurrent readers.
!> \param cache ...
!> \param row local grid row
!> \param point current Cartesian coordinates, including atom motion and quadrature changes
!> \param stencil cached interpolation stencil on a hit
!> \return whether the cached stencil matches
! **************************************************************************************************
   FUNCTION fetch_native_grid_stencil(cache, row, point, stencil) RESULT(hit)
      TYPE(native_grid_cache_type), INTENT(IN)           :: cache
      INTEGER, INTENT(IN)                                :: row
      REAL(dp), INTENT(IN)                               :: point(3)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(OUT)                                     :: stencil
      LOGICAL                                            :: hit

      hit = .FALSE.
      IF (.NOT. cache%prepared .OR. .NOT. ALLOCATED(cache%entries)) RETURN
      IF (row < 1 .OR. row > SIZE(cache%entries)) RETURN
      IF (.NOT. cache%entries(row)%ready) RETURN
      IF (ANY(cache%entries(row)%point /= point)) RETURN
      stencil = cache%entries(row)%stencil
      hit = .TRUE.
   END FUNCTION fetch_native_grid_stencil

! **************************************************************************************************
!> \brief Store from the unique forward owner of a row, never from concurrent adjoint readers.
!> \param cache ...
!> \param row local grid row
!> \param point current Cartesian coordinates
!> \param stencil complete value stencil, not an indices-only stencil
! **************************************************************************************************
   SUBROUTINE store_native_grid_stencil(cache, row, point, stencil)
      TYPE(native_grid_cache_type), INTENT(INOUT)        :: cache
      INTEGER, INTENT(IN)                                :: row
      REAL(dp), INTENT(IN)                               :: point(3)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil

      IF (.NOT. cache%prepared .OR. .NOT. ALLOCATED(cache%entries)) RETURN
      IF (row < 1 .OR. row > SIZE(cache%entries)) RETURN
      cache%entries(row)%stencil = stencil
      cache%entries(row)%point = point
      cache%entries(row)%ready = .TRUE.
   END SUBROUTINE store_native_grid_stencil
END MODULE qs_native_grid_cache

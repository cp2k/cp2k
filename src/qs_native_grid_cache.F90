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

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: native_grid_interp_offset_min = -5, native_grid_interp_offset_max = 6
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
   PUBLIC :: prepare_native_grid_cache, release_native_grid_cache, fetch_native_grid_stencil, store_native_grid_stencil
CONTAINS
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

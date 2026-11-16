!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common selection and union operations for contiguous molecular-orbital windows.
! **************************************************************************************************
MODULE mo_window
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mo_window'

   TYPE mo_window_type
      INTEGER :: first_mo = 1
      INTEGER :: last_mo = 0
   END TYPE mo_window_type

   PUBLIC :: combine_mo_windows, determine_mo_window, mo_window_type

CONTAINS

! **************************************************************************************************
!> \brief Selects one contiguous MO window from an ordered reference spectrum.
!>
!> The input spectrum is supplied by the caller because different workflows may store more
!> entries than the number of physical MOs to which a cutoff applies. The returned bounds are
!> absolute, one-based indices in the caller's MO numbering.
!> \param eigenvalues Reference eigenvalues for one spin channel.
!> \param n_mo Number of physical MOs covered by the reference spectrum.
!> \param n_occ Number of occupied MOs in that channel.
!> \param cutoff_occ Occupied-state energy window; non-positive means no occupied cutoff.
!> \param cutoff_empty Empty-state energy window; non-positive means no empty-state cutoff.
!> \param window Selected absolute MO window (OUT).
! **************************************************************************************************
   SUBROUTINE determine_mo_window(eigenvalues, n_mo, n_occ, cutoff_occ, cutoff_empty, window)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      INTEGER, INTENT(IN)                                :: n_mo, n_occ
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_occ, cutoff_empty
      TYPE(mo_window_type), INTENT(OUT)                  :: window

      INTEGER                                            :: i

      IF (n_mo < 1 .OR. n_mo > SIZE(eigenvalues)) THEN
         CALL cp_abort(__LOCATION__, "determine_mo_window: invalid number of physical MOs")
      END IF
      IF (n_occ < 1 .OR. n_occ >= n_mo) THEN
         CALL cp_abort(__LOCATION__, "determine_mo_window: invalid occupied-state boundary")
      END IF

      window%first_mo = 1
      window%last_mo = n_mo

      IF (cutoff_occ > 0.0_dp .OR. cutoff_empty > 0.0_dp) THEN
         DO i = 2, n_mo
            IF (eigenvalues(i) < eigenvalues(i - 1)) THEN
               CALL cp_abort(__LOCATION__, &
                             "determine_mo_window: reference eigenvalues are not ascending. "// &
                             "Use the DFT/SCF energy axis; the G0W0 axis is not ordered.")
            END IF
         END DO

         IF (cutoff_occ > 0.0_dp) THEN
            DO i = 1, n_occ
               IF (eigenvalues(n_occ) - eigenvalues(i) <= cutoff_occ) THEN
                  window%first_mo = i
                  EXIT
               END IF
            END DO
         END IF

         IF (cutoff_empty > 0.0_dp) THEN
            DO i = n_occ + 1, n_mo
               IF (eigenvalues(i) - eigenvalues(n_occ + 1) > cutoff_empty) THEN
                  window%last_mo = i - 1
                  EXIT
               END IF
            END DO
         END IF
      END IF

   END SUBROUTINE determine_mo_window

! **************************************************************************************************
!> \brief Forms the smallest contiguous window covering all supplied spin windows.
!> \param spin_windows Per-spin absolute MO windows.
!> \param combined_window Union of the supplied windows (OUT).
!> \param windows_differ True if at least two input windows have different bounds (OUT).
! **************************************************************************************************
   SUBROUTINE combine_mo_windows(spin_windows, combined_window, windows_differ)
      TYPE(mo_window_type), DIMENSION(:), INTENT(IN)     :: spin_windows
      TYPE(mo_window_type), INTENT(OUT)                  :: combined_window
      LOGICAL, INTENT(OUT)                               :: windows_differ

      INTEGER                                            :: i

      IF (SIZE(spin_windows) < 1) THEN
         CALL cp_abort(__LOCATION__, "combine_mo_windows: no spin windows supplied")
      END IF

      combined_window%first_mo = spin_windows(1)%first_mo
      combined_window%last_mo = spin_windows(1)%last_mo
      windows_differ = .FALSE.

      DO i = 2, SIZE(spin_windows)
         combined_window%first_mo = MIN(combined_window%first_mo, spin_windows(i)%first_mo)
         combined_window%last_mo = MAX(combined_window%last_mo, spin_windows(i)%last_mo)
         windows_differ = windows_differ .OR. &
                          spin_windows(i)%first_mo /= spin_windows(1)%first_mo .OR. &
                          spin_windows(i)%last_mo /= spin_windows(1)%last_mo
      END DO

   END SUBROUTINE combine_mo_windows

END MODULE mo_window

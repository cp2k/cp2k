!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_frame_controls
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: frame_selected

CONTAINS

   PURE LOGICAL FUNCTION frame_selected(frame, first, last, stride)
      INTEGER, INTENT(IN)                                :: frame, first, last, stride

      frame_selected = frame >= first .AND. frame <= last .AND. MOD(frame - first, stride) == 0
   END FUNCTION frame_selected

END MODULE trajana_frame_controls

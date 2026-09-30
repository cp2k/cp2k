!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM main
   IMPLICIT NONE
   INTEGER :: values(4) = [2, 3, 2, 4]

   ! FINDLOC with DIM is used to convert active-space orbital indices.
   IF (FINDLOC(values, 2, dim=1) /= 1) ERROR STOP 1
   IF (FINDLOC(values, 5, dim=1) /= 0) ERROR STOP 2
END PROGRAM

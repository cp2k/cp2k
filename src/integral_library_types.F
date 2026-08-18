!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Library choices for electronic integral APIs.
! **************************************************************************************************
MODULE integral_library_types

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: library_native = 0, &
                                 library_libint = 1

   TYPE, PUBLIC :: integral_library_type
      INTEGER :: coulomb2_library = library_native
      INTEGER :: coulomb3_library = library_native
   END TYPE integral_library_type

   TYPE(integral_library_type), SAVE, PUBLIC :: active_integral_library = integral_library_type()

END MODULE integral_library_types

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Backend choices for electronic integral APIs.
! **************************************************************************************************
MODULE electronic_integral_backend_types

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: backend_native = 0, &
                                 backend_libint = 1

   TYPE, PUBLIC :: electronic_integral_backend_type
      INTEGER :: coulomb2_backend = backend_native
      INTEGER :: coulomb3_backend = backend_native
   END TYPE electronic_integral_backend_type

END MODULE electronic_integral_backend_types

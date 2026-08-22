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
   USE input_constants,                 ONLY: do_potential_coulomb
   USE kinds,                           ONLY: default_path_length,&
                                              dp

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: library_native = 0, &
                                 library_libint = 1

   TYPE, PUBLIC :: integral_library_type
      INTEGER :: coulomb2_library = library_native
      INTEGER :: coulomb3_library = library_native
   END TYPE integral_library_type

   TYPE, PUBLIC :: coulomb_operator_type
      INTEGER                        :: potential_type = do_potential_coulomb
      REAL(dp)                       :: omega = 0.0_dp ! SR: erfc(omega*r)/r
      REAL(dp)                       :: cutoff_radius = 0.0_dp ! TC cutoff/effective SR range
      CHARACTER(default_path_length) :: filename = ""
      REAL(dp)                       :: scale_coulomb = 0.0_dp ! Only for WFC methods
      REAL(dp)                       :: scale_longrange = 0.0_dp ! Only for WFC methods
   END TYPE coulomb_operator_type

   TYPE(integral_library_type), SAVE, PUBLIC :: active_integral_library = integral_library_type()

END MODULE integral_library_types

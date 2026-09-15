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

#if defined(__LIBINT)
   INTEGER, PARAMETER, PUBLIC :: default_coulomb_library = library_libint
   LOGICAL, PARAMETER         :: libint_available = .TRUE.
#else
   INTEGER, PARAMETER, PUBLIC :: default_coulomb_library = library_native
   LOGICAL, PARAMETER         :: libint_available = .FALSE.
#endif

   TYPE, PUBLIC :: integral_library_type
      INTEGER :: coulomb2_library = default_coulomb_library
      INTEGER :: coulomb3_library = default_coulomb_library
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

   PUBLIC :: integral_library_init

CONTAINS

! **************************************************************************************************
!> \brief Resolves requested Coulomb integral libraries against the compiled-in backends.
!> \param coulomb2_library requested two-center Coulomb integral library
!> \param coulomb3_library requested three-center Coulomb integral library
!> \param fallback_applied whether an unavailable Libint request was replaced by the native backend
! **************************************************************************************************
   SUBROUTINE integral_library_init(coulomb2_library, coulomb3_library, fallback_applied)
      INTEGER, INTENT(IN)                                :: coulomb2_library, coulomb3_library
      LOGICAL, INTENT(OUT)                               :: fallback_applied

      active_integral_library%coulomb2_library = coulomb2_library
      active_integral_library%coulomb3_library = coulomb3_library
      fallback_applied = .FALSE.

      IF (.NOT. libint_available) THEN
         IF (active_integral_library%coulomb2_library == library_libint) THEN
            active_integral_library%coulomb2_library = library_native
            fallback_applied = .TRUE.
         END IF
         IF (active_integral_library%coulomb3_library == library_libint) THEN
            active_integral_library%coulomb3_library = library_native
            fallback_applied = .TRUE.
         END IF
      END IF

   END SUBROUTINE integral_library_init

END MODULE integral_library_types

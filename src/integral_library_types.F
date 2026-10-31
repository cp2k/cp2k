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
   USE input_constants,                 ONLY: &
        do_potential_TShPSC, do_potential_coulomb, do_potential_gaussian, do_potential_long, &
        do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, do_potential_short, &
        do_potential_truncated
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
      REAL(dp)                       :: scale_gaussian = 0.0_dp ! Gaussian coefficient in mixed operators
   END TYPE coulomb_operator_type

   ! Enlarge the support slightly when screening the discontinuous truncated Coulomb operator.
   REAL(dp), PARAMETER, PUBLIC :: cutoff_screen_factor = 1.0001_dp

   TYPE(integral_library_type), SAVE, PUBLIC :: active_integral_library = integral_library_type()

   PUBLIC :: compare_potential_types, integral_library_init

CONTAINS

! **************************************************************************************************
!> \brief Compare the parameters defining two radial operators.
!> \param potential1 first operator
!> \param potential2 second operator
!> \return whether the operator selectors and their relevant parameters agree
!> \note Table filenames and effective short-range screening cutoffs are excluded.
! **************************************************************************************************
   PURE FUNCTION compare_potential_types(potential1, potential2) RESULT(equals)
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential1, potential2
      LOGICAL                                            :: equals

      equals = .FALSE.
      IF (potential1%potential_type /= potential2%potential_type) RETURN

      equals = .TRUE.
      SELECT CASE (potential1%potential_type)
      CASE (do_potential_short, do_potential_long, do_potential_gaussian)
         IF (potential1%omega /= potential2%omega) equals = .FALSE.
      CASE (do_potential_truncated, do_potential_TShPSC)
         IF (potential1%cutoff_radius /= potential2%cutoff_radius) equals = .FALSE.
      CASE (do_potential_mix_cl, do_potential_mix_cl_trunc)
         IF (potential1%omega /= potential2%omega) equals = .FALSE.
         IF (potential1%scale_coulomb /= potential2%scale_coulomb) equals = .FALSE.
         IF (potential1%scale_longrange /= potential2%scale_longrange) equals = .FALSE.
         IF (potential1%potential_type == do_potential_mix_cl_trunc) THEN
            IF (potential1%cutoff_radius /= potential2%cutoff_radius) equals = .FALSE.
         END IF
      CASE (do_potential_mix_lg)
         IF (potential1%omega /= potential2%omega) equals = .FALSE.
         IF (potential1%scale_longrange /= potential2%scale_longrange) equals = .FALSE.
         IF (potential1%scale_gaussian /= potential2%scale_gaussian) equals = .FALSE.
      END SELECT
   END FUNCTION compare_potential_types

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

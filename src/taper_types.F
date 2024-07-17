!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of the semi empirical parameter types.
!> \author Teodoro Laino [tlaino] - 10.2008 University of Zurich
! **************************************************************************************************
MODULE taper_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'taper_types'

! **************************************************************************************************
!> \brief Taper type
! **************************************************************************************************
   TYPE taper_type
      LOGICAL                               :: apply_taper = .FALSE.
      REAL(KIND=dp)                         :: r0 = -1.0_dp, rscale = -1.0_dp
   END TYPE taper_type

   PUBLIC :: taper_type, taper_create, taper_release, taper_eval, dtaper_eval

CONTAINS

! **************************************************************************************************
!> \brief Creates taper type
!> \param taper ...
!> \param rc ...
!> \param range ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE taper_create(taper, rc, range)
      TYPE(taper_type), INTENT(OUT)                      :: taper
      REAL(KIND=dp), INTENT(IN)                          :: rc, range

      IF (range > EPSILON(0.0_dp)) THEN
         taper%apply_taper = .TRUE.
         taper%r0 = 2.0_dp*rc - 20.0_dp*range
         taper%rscale = 1.0_dp/range
      ELSE
         taper%apply_taper = .FALSE.
      END IF

   END SUBROUTINE taper_create

! **************************************************************************************************
!> \brief Releases taper type
!> \param taper ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE taper_release(taper)
      TYPE(taper_type), INTENT(INOUT)                    :: taper

      MARK_USED(taper)

   END SUBROUTINE taper_release

! **************************************************************************************************
!> \brief Taper functions
!> \param taper ...
!> \param rij ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION taper_eval(taper, rij) RESULT(ft)
      TYPE(taper_type), INTENT(IN)                       :: taper
      REAL(KIND=dp), INTENT(IN)                          :: rij
      REAL(KIND=dp)                                      :: ft

      REAL(KIND=dp)                                      :: dr

      ft = 1._dp
      IF (taper%apply_taper) THEN
         dr = taper%rscale*(rij - taper%r0)
         ft = 0.5_dp*(1.0_dp - TANH(dr))
      END IF
   END FUNCTION taper_eval

! **************************************************************************************************
!> \brief Analytical derivatives for taper function
!> \param taper ...
!> \param rij ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION dtaper_eval(taper, rij) RESULT(dft)
      TYPE(taper_type), INTENT(IN)                       :: taper
      REAL(KIND=dp), INTENT(IN)                          :: rij
      REAL(KIND=dp)                                      :: dft

      REAL(KIND=dp)                                      :: dr

      dft = 0.0_dp
      IF (taper%apply_taper) THEN
         dr = taper%rscale*(rij - taper%r0)
         dft = -0.5_dp*(1.0_dp - TANH(dr)**2)*taper%rscale
      END IF
   END FUNCTION dtaper_eval

END MODULE taper_types

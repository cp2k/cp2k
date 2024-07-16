!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of the semi empirical multipole integral expansions types
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
MODULE semi_empirical_mpole_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_mpole_types'

! **************************************************************************************************
!> \brief Semi-empirical integral multipole expansion type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   TYPE semi_empirical_mpole_type
      LOGICAL, DIMENSION(3)                    :: task = .FALSE.
      INTEGER                                  :: indi = -1, indj = -1
      REAL(KIND=dp)                            :: c = -1.0_dp
      REAL(KIND=dp), DIMENSION(3)              :: d = -1.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)            :: qc = -1.0_dp ! quadrupole cartesian
      REAL(KIND=dp), DIMENSION(5)              :: qs = -1.0_dp ! quadrupole spherical
      ! alternative definition used in GKS integral routines
      REAL(KIND=dp)                            :: cs = -1.0_dp
      REAL(KIND=dp), DIMENSION(3)              :: ds = -1.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)            :: qq = -1.0_dp ! quadrupole cartesian
   END TYPE semi_empirical_mpole_type

! **************************************************************************************************
!> \brief Global Multipolar NDDO information type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   TYPE nddo_mpole_type
      REAL(KIND=dp), DIMENSION(:), POINTER     :: charge => NULL(), efield0 => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: dipole => NULL(), efield1 => NULL(), efield2 => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: quadrupole => NULL()
   END TYPE nddo_mpole_type

   PUBLIC :: semi_empirical_mpole_type, &
             semi_empirical_mpole_p_create, &
             semi_empirical_mpole_p_release, &
             nddo_mpole_type, &
             nddo_mpole_create, &
             nddo_mpole_release

CONTAINS

! **************************************************************************************************
!> \brief Allocate semi-empirical mpole type
!> \param mpole ...
!> \param ndim ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   SUBROUTINE semi_empirical_mpole_p_create(mpole, ndim)
      TYPE(semi_empirical_mpole_type), DIMENSION(:), &
         POINTER                                         :: mpole
      INTEGER, INTENT(IN)                                :: ndim

      INTEGER                                            :: i

      CPASSERT(.NOT. ASSOCIATED(mpole))
      ALLOCATE (mpole(ndim))
      DO i = 1, ndim
         CALL semi_empirical_mpole_create(mpole(i))
      END DO

   END SUBROUTINE semi_empirical_mpole_p_create

! **************************************************************************************************
!> \brief Deallocate the semi-empirical mpole type
!> \param mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   SUBROUTINE semi_empirical_mpole_p_release(mpole)
      TYPE(semi_empirical_mpole_type), DIMENSION(:), &
         POINTER                                         :: mpole

      INTEGER                                            :: i

      IF (ASSOCIATED(mpole)) THEN
         DO i = 1, SIZE(mpole)
            CALL semi_empirical_mpole_release(mpole(i))
         END DO
         DEALLOCATE (mpole)
      END IF

   END SUBROUTINE semi_empirical_mpole_p_release

! **************************************************************************************************
!> \brief Allocate semi-empirical mpole type
!> \param mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   SUBROUTINE semi_empirical_mpole_create(mpole)
      TYPE(semi_empirical_mpole_type), INTENT(INOUT)     :: mpole

      mpole%task = .FALSE.
      mpole%indi = 0
      mpole%indj = 0
      mpole%c = HUGE(0.0_dp)
      mpole%d = HUGE(0.0_dp)
      mpole%qc = HUGE(0.0_dp)
      mpole%qs = HUGE(0.0_dp)
      mpole%cs = HUGE(0.0_dp)
      mpole%ds = HUGE(0.0_dp)
      mpole%qq = HUGE(0.0_dp)
   END SUBROUTINE semi_empirical_mpole_create

! **************************************************************************************************
!> \brief Deallocate the semi-empirical mpole type
!> \param mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   SUBROUTINE semi_empirical_mpole_release(mpole)
      TYPE(semi_empirical_mpole_type), INTENT(INOUT)     :: mpole

      MARK_USED(mpole)

   END SUBROUTINE semi_empirical_mpole_release

! **************************************************************************************************
!> \brief Allocate NDDO multipole type
!> \param nddo_mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   SUBROUTINE nddo_mpole_create(nddo_mpole)
      TYPE(nddo_mpole_type), INTENT(OUT)                 :: nddo_mpole

      NULLIFY (nddo_mpole%charge)
      NULLIFY (nddo_mpole%dipole)
      NULLIFY (nddo_mpole%quadrupole)
      NULLIFY (nddo_mpole%efield0)
      NULLIFY (nddo_mpole%efield1)
      NULLIFY (nddo_mpole%efield2)
   END SUBROUTINE nddo_mpole_create

! **************************************************************************************************
!> \brief Deallocate NDDO multipole type
!> \param nddo_mpole ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! **************************************************************************************************
   SUBROUTINE nddo_mpole_release(nddo_mpole)
      TYPE(nddo_mpole_type), INTENT(INOUT)               :: nddo_mpole

      IF (ASSOCIATED(nddo_mpole%charge)) THEN
         DEALLOCATE (nddo_mpole%charge)
      END IF
      IF (ASSOCIATED(nddo_mpole%dipole)) THEN
         DEALLOCATE (nddo_mpole%dipole)
      END IF
      IF (ASSOCIATED(nddo_mpole%quadrupole)) THEN
         DEALLOCATE (nddo_mpole%quadrupole)
      END IF
      IF (ASSOCIATED(nddo_mpole%efield0)) THEN
         DEALLOCATE (nddo_mpole%efield0)
      END IF
      IF (ASSOCIATED(nddo_mpole%efield1)) THEN
         DEALLOCATE (nddo_mpole%efield1)
      END IF
      IF (ASSOCIATED(nddo_mpole%efield2)) THEN
         DEALLOCATE (nddo_mpole%efield2)
      END IF

   END SUBROUTINE nddo_mpole_release

END MODULE semi_empirical_mpole_types

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GRRM interface
!> \author JGH - 08.2019
! **************************************************************************************************
MODULE grrm_utils

   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grrm_utils'

   PUBLIC :: write_grrm

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Write GRRM interface file
!>
!> \param iounit ...
!> \param particles ...
!> \param energy ...
!> \param dipole ...
!> \param hessian ...
!> \param dipder ...
!> \param polar ...
! **************************************************************************************************
   SUBROUTINE write_grrm(iounit, particles, energy, dipole, hessian, dipder, polar)

      INTEGER, INTENT(IN)                                :: iounit
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particles
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: dipole
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: hessian, dipder
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: polar

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_grrm', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      INTEGER                                            :: i, j, natom, nc
      REAL(KIND=dp), DIMENSION(5)                        :: fz

      IF (iounit > 0) THEN
         natom = SIZE(particles)
         WRITE (iounit, "(A7)") "RESULTS"
         WRITE (iounit, "(A18)") "CURRENT COORDINATE"
         DO i = 1, natom
            WRITE (iounit, "(A,3F24.12)") TRIM(ADJUSTL(particles(i)%atomic_kind%element_symbol)), &
               particles(i)%r(1:3)*angstrom
         END DO
         WRITE (iounit, "(A8,3F18.12)") "ENERGY =", energy, zero, zero
         WRITE (iounit, "(A8,3F18.12)") "       =", zero, zero, zero
         WRITE (iounit, "(A8,F18.12)") "S**2   =", zero
         WRITE (iounit, "(A8)") "GRADIENT"
         DO i = 1, natom
            WRITE (iounit, "(F17.12)") particles(i)%f(1:3)
         END DO
         IF (PRESENT(dipole)) THEN
            WRITE (iounit, "(A8,3F18.12)") "DIPOLE =", dipole(1:3)
         ELSE
            WRITE (iounit, "(A8,3F18.12)") "DIPOLE =", zero, zero, zero
         END IF
         fz = zero
         WRITE (iounit, "(A7)") "HESSIAN"
         IF (PRESENT(hessian)) THEN
            nc = SIZE(hessian, 1)
            CPASSERT(nc == 3*natom)
            DO i = 1, nc, 5
               DO j = i, nc
                  WRITE (iounit, "(5(F13.9,1X))") hessian(j, i:MIN(j, i + 4))
               END DO
            END DO
         ELSE
            nc = 3*natom
            DO i = 1, nc, 5
               DO j = i, nc
                  WRITE (iounit, "(5(F13.9,1X))") fz(1:MIN(j - i + 1, 5))
               END DO
            END DO
         END IF
         WRITE (iounit, "(A18)") "DIPOLE DERIVATIVES"
         IF (PRESENT(dipder)) THEN
            DO i = 1, 3*natom
               WRITE (iounit, "(3(F17.12,7X))") dipder(1:3, i)
            END DO
         ELSE
            DO i = 1, 3*natom
               WRITE (iounit, "(3(F17.12,7X))") zero, zero, zero
            END DO
         END IF
         WRITE (iounit, "(A14)") "POLARIZABILITY"
         IF (PRESENT(polar)) THEN
            WRITE (iounit, "(1(F17.12,7X))") polar(1, 1)
            WRITE (iounit, "(2(F17.12,7X))") polar(2, 1), polar(2, 2)
            WRITE (iounit, "(3(F17.12,7X))") polar(3, 1), polar(3, 2), polar(3, 3)
         ELSE
            WRITE (iounit, "(1(F17.12,7X))") zero
            WRITE (iounit, "(2(F17.12,7X))") zero, zero
            WRITE (iounit, "(3(F17.12,7X))") zero, zero, zero
         END IF
      END IF

   END SUBROUTINE write_grrm

END MODULE grrm_utils

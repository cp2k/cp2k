!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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

      INTEGER                                            :: i, j, natom, nc
      REAL(KIND=dp)                                      :: zero = 0.0_dp
      REAL(KIND=dp), DIMENSION(5)                        :: fz

      IF (iounit > 0) THEN
         natom = SIZE(particles)
         WRITE (iounit, "(A7)") "RESULTS"
         WRITE (iounit, "(A18)") "CURRENT COORDINATE"
         DO i = 1, natom
            WRITE (iounit, "(A2,3F24.12)") particles(i)%atomic_kind%element_symbol, &
               particles(i)%r(1:3)*angstrom
         END DO
         WRITE (iounit, "(A7,3F17.12)") "ENERGY=", energy, zero, zero
         WRITE (iounit, "(A7,3F17.12)") "      =", zero, zero, zero
         WRITE (iounit, "(A7,F17.12)") "S**2  =", zero
         WRITE (iounit, "(A8)") "GRADIENT"
         DO i = 1, natom
            WRITE (iounit, "(F17.12)")-particles(i)%f(1:3)
         END DO
         IF (PRESENT(dipole)) THEN
            WRITE (iounit, "(A7,3F17.12)") "DIPOLE=", dipole(1:3)
         ELSE
            WRITE (iounit, "(A7,3F17.12)") "DIPOLE=", zero, zero, zero
         END IF
         fz = zero
         WRITE (iounit, "(A7)") "HESSIAN"
         IF (PRESENT(hessian)) THEN
            nc = SIZE(hessian, 1)
            CPASSERT(nc == 3*natom)
            DO i = 1, nc, 5
               DO j = i, nc
                  WRITE (iounit, "(5F17.12)") hessian(j, i:MIN(j, i+4))
               END DO
            END DO
         ELSE
            nc = 3*natom
            DO i = 1, nc, 5
               DO j = i, nc
                  WRITE (iounit, "(5F17.12)") fz(1:MIN(j-i+1, 5))
               END DO
            END DO
         END IF
         WRITE (iounit, "(A17)") "DIPOLEDERIVATIVES"
         IF (PRESENT(dipder)) THEN
            DO i = 1, 3*natom
               WRITE (iounit, "(3F17.12)") dipder(1:3, i)
            END DO
         ELSE
            DO i = 1, 3*natom
               WRITE (iounit, "(3F17.12)") zero, zero, zero
            END DO
         END IF
         WRITE (iounit, "(A14)") "POLARIZABILITY"
         IF (PRESENT(polar)) THEN
            WRITE (iounit, "(1F17.12)") polar(1, 1)
            WRITE (iounit, "(2F17.12)") polar(2, 1), polar(2, 2)
            WRITE (iounit, "(3F17.12)") polar(3, 1), polar(3, 2), polar(3, 3)
         ELSE
            WRITE (iounit, "(1F17.12)") zero
            WRITE (iounit, "(2F17.12)") zero, zero
            WRITE (iounit, "(3F17.12)") zero, zero, zero
         END IF
      END IF

   END SUBROUTINE write_grrm

END MODULE grrm_utils

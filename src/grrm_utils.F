!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GRRM interface
!> \author JGH - 08.2019
! **************************************************************************************************
MODULE grrm_utils

   USE cp_control_types,                ONLY: dft_control_type
   USE force_env_types,                 ONLY: force_env_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env
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
!> \param force_env ...
!> \param particles ...
!> \param energy ...
!> \param dipole ...
!> \param hessian ...
!> \param dipder ...
!> \param polar ...
!> \param fixed_atoms ...
! **************************************************************************************************
   SUBROUTINE write_grrm(iounit, force_env, particles, energy, dipole, hessian, dipder, polar, &
                         fixed_atoms)

      INTEGER, INTENT(IN)                                :: iounit
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particles
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: dipole
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: hessian, dipder
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: polar
      INTEGER, INTENT(IN), OPTIONAL                      :: fixed_atoms

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_grrm', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      INTEGER                                            :: i, j, natom, nc
      LOGICAL                                            :: nddo
      REAL(KIND=dp)                                      :: eout
      REAL(KIND=dp), DIMENSION(5)                        :: fz
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: qs_energy

      IF (iounit > 0) THEN
         ! the units depend on the qs method!
         CPASSERT(ASSOCIATED(force_env%qs_env))
         CALL get_qs_env(force_env%qs_env, dft_control=dft_control)
         nddo = dft_control%qs_control%semi_empirical
         IF (nddo) THEN
            CALL get_qs_env(force_env%qs_env, energy=qs_energy)
            eout = energy + qs_energy%core_self
         ELSE
            eout = energy
         END IF
         !
         natom = SIZE(particles)
         IF (PRESENT(fixed_atoms)) natom = natom - fixed_atoms
         WRITE (iounit, "(A7)") "RESULTS"
         WRITE (iounit, "(A18)") "CURRENT COORDINATE"
         DO i = 1, natom
            WRITE (iounit, "(A,3F24.12)") TRIM(ADJUSTL(particles(i)%atomic_kind%element_symbol)), &
               particles(i)%r(1:3)*angstrom
         END DO
         WRITE (iounit, "(A8,3F18.12)") "ENERGY =", eout, zero, zero
         WRITE (iounit, "(A8,3F18.12)") "       =", zero, zero, zero
         WRITE (iounit, "(A8,F18.12)") "S**2   =", zero
         WRITE (iounit, "(A8)") "GRADIENT"
         DO i = 1, natom
            WRITE (iounit, "(F17.12)") - particles(i)%f(1:3)
         END DO
         IF (PRESENT(dipole)) THEN
            WRITE (iounit, "(A8,3F18.12)") "DIPOLE =", dipole(1:3)
         ELSE
            WRITE (iounit, "(A8,3F18.12)") "DIPOLE =", zero, zero, zero
         END IF
         fz = zero
         WRITE (iounit, "(A7)") "HESSIAN"
         IF (PRESENT(hessian)) THEN
            nc = 3*natom
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

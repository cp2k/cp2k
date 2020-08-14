!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief SCINE interface
!> \author JGH - 01.2020
! **************************************************************************************************
MODULE scine_utils

   USE cell_types,                      ONLY: cell_type
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scine_utils'

   PUBLIC :: write_scine

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Write SCINE interface file
!>
!> \param iounit ...
!> \param force_env ...
!> \param particles ...
!> \param energy ...
!> \param hessian ...
! **************************************************************************************************
   SUBROUTINE write_scine(iounit, force_env, particles, energy, hessian)

      INTEGER, INTENT(IN)                                :: iounit
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particles
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: hessian

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_scine'
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      CHARACTER(LEN=20)                                  :: sunit
      INTEGER                                            :: i, j, natom, nc
      LOGICAL                                            :: nddo
      REAL(KIND=dp)                                      :: eout
      REAL(KIND=dp), DIMENSION(5)                        :: fz
      TYPE(cell_type), POINTER                           :: cell
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
         CALL get_qs_env(force_env%qs_env, cell=cell)
         natom = SIZE(particles)
         WRITE (iounit, "(A,A)") "title: ", "'SCINE Interface File'"
         sunit = "'hartree'"
         WRITE (iounit, "(A,F18.12,A,A)") "energy: [", eout, ", "//TRIM(sunit), " ]"
         sunit = "'angstrom'"
         WRITE (iounit, "(A)") "system:"
         WRITE (iounit, "(T2,A,A)") "unit: ", TRIM(sunit)
         WRITE (iounit, "(T2,A)") "cell: "
         WRITE (iounit, "(T4,A,F24.12,',',F24.12,',',F24.12,A)") &
            "- A:  [", cell%hmat(1:3, 1)*angstrom, " ]"
         WRITE (iounit, "(T4,A,F24.12,',',F24.12,',',F24.12,A)") &
            "- B:  [", cell%hmat(1:3, 2)*angstrom, " ]"
         WRITE (iounit, "(T4,A,F24.12,',',F24.12,',',F24.12,A)") &
            "- C:  [", cell%hmat(1:3, 3)*angstrom, " ]"
         WRITE (iounit, "(T2,A,L1,', ',L1,', ',L1,' ]')") "periodicity:  [ ", (cell%perd == 1)
         WRITE (iounit, "(T2,A)") "coordinates: "
         DO i = 1, natom
            WRITE (iounit, "(T4,A,A2,A,F24.12,',',F24.12,',',F24.12,A)") &
               "- ", TRIM(ADJUSTL(particles(i)%atomic_kind%element_symbol))//":", &
               "  [", particles(i)%r(1:3)*angstrom, " ]"
         END DO
         WRITE (iounit, "(A)") "gradient:"
         sunit = "'hartree/bohr'"
         WRITE (iounit, "(T2,A,A)") "unit: ", TRIM(sunit)
         WRITE (iounit, "(T2,A)") "values: "
         DO i = 1, natom
            WRITE (iounit, "(T4,A,A2,A,F24.12,',',F24.12,',',F24.12,A)") &
               "- ", TRIM(ADJUSTL(particles(i)%atomic_kind%element_symbol))//":", &
               "  [", particles(i)%f(1:3), " ]"
         END DO
         fz = zero
         IF (PRESENT(hessian)) THEN
            WRITE (iounit, "(A)") "hessian:"
            sunit = "'hartree/bohr^2'"
            WRITE (iounit, "(T2,A,A)") "unit: ", TRIM(sunit)
            DO i = 1, natom
               WRITE (iounit, "(T4,A)") TRIM(ADJUSTL(particles(i)%atomic_kind%element_symbol))//":"
               nc = 3*(i - 1) + 1
               WRITE (iounit, "(T6,'X:  [')")
               DO j = 1, nc, 5
                  IF (nc > j + 4) THEN
                     WRITE (iounit, "(T12,5(F20.12,','))") hessian(nc, j:j + 4)
                  ELSEIF (nc == j + 4) THEN
                     WRITE (iounit, "(T12,4(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 3) THEN
                     WRITE (iounit, "(T12,3(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 2) THEN
                     WRITE (iounit, "(T12,2(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 1) THEN
                     WRITE (iounit, "(T12,1(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j) THEN
                     WRITE (iounit, "(T12,F20.12,' ]')") hessian(nc, j:nc)
                  END IF
               END DO
               nc = 3*(i - 1) + 2
               WRITE (iounit, "(T6,'Y:  [')")
               DO j = 1, nc, 5
                  IF (nc > j + 4) THEN
                     WRITE (iounit, "(T12,5(F20.12,','))") hessian(nc, j:j + 4)
                  ELSEIF (nc == j + 4) THEN
                     WRITE (iounit, "(T12,4(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 3) THEN
                     WRITE (iounit, "(T12,3(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 2) THEN
                     WRITE (iounit, "(T12,2(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 1) THEN
                     WRITE (iounit, "(T12,1(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j) THEN
                     WRITE (iounit, "(T12,F20.12,' ]')") hessian(nc, j:nc)
                  END IF
               END DO
               nc = 3*(i - 1) + 3
               WRITE (iounit, "(T6,'Z:  [')")
               DO j = 1, nc, 5
                  IF (nc > j + 4) THEN
                     WRITE (iounit, "(T12,5(F20.12,','))") hessian(nc, j:j + 4)
                  ELSEIF (nc == j + 4) THEN
                     WRITE (iounit, "(T12,4(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 3) THEN
                     WRITE (iounit, "(T12,3(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 2) THEN
                     WRITE (iounit, "(T12,2(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j + 1) THEN
                     WRITE (iounit, "(T12,1(F20.12,','),F20.12,' ]')") hessian(nc, j:nc)
                  ELSEIF (nc == j) THEN
                     WRITE (iounit, "(T12,F20.12,' ]')") hessian(nc, j:nc)
                  END IF
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE write_scine

END MODULE scine_utils

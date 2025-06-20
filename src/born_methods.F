MODULE born_methods
   USE cp_log_handling, ONLY: cp_logger_type, cp_get_default_logger, &
                        cp_logger_get_default_io_unit
   USE kinds, ONLY: default_string_length, dp

   USE born_types, ONLY: born_samples_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'born_methods'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   PUBLIC :: do_born

CONTAINS
   SUBROUTINE do_born(born_samples, born_tensor)
      TYPE(born_samples_type), POINTER                      :: born_samples
      REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE         :: born_tensor
      !REAL(kind=dp), DIMENSION(3,3), ALLOCATABLE         :: polarizabilities !TODO:
      !REAL(kind=dp), DIMENSION(3,3), ALLOCATABLE         :: dielectric


      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'do_born'

      INTEGER                                            :: handle, i, j, log_unit, n
      TYPE(cp_logger_type), POINTER                      :: logger
      !TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      !TYPE(qs_environment_type), POINTER                 :: qs_env



      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      log_unit = cp_logger_get_default_io_unit()


      ALLOCATE(born_tensor(SIZE(born_samples%sample_no_field%forces, 1), 3, 3))

      DO j= 1, 3 ! axis force
         DO i = 1, 3 ! axis field
            DO n = 1, SIZE(born_samples%sample_no_field%forces, 1) ! atom
               born_tensor(n, i, j) = ( &
                  born_samples%samples_field(i, 1)%forces(n, j) - &
                  born_samples%samples_field(i, 2)%forces(n, j) &
               ) / 2 * born_samples%field_strength
            END DO
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE do_born
END MODULE

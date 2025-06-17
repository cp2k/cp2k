MODULE born_methods
   USE force_env_types, ONLY: force_env_get, force_env_type
   USE cp_log_handling, ONLY: cp_logger_type, cp_get_default_logger, &
                        cp_logger_get_default_io_unit
   USE cp_result_methods, ONLY: get_results
   USE cp_result_types, ONLY: cp_result_type
   USE cp_subsys_types, ONLY: cp_subsys_get, cp_subsys_type
   USE kinds, ONLY: default_string_length, dp
   USE particle_list_types, ONLY: particle_list_type
   USE physcon, ONLY: debye
   !USE qs_environment_types, ONLY: qs_environment_type


#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'born_methods'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   PUBLIC :: get_forces, do_born

CONTAINS
   SUBROUTINE get_forces(force_env, dipole, forces)
      TYPE(force_env_type), POINTER                      :: force_env
      REAL(kind=dp), DIMENSION(1:3)                      :: dipole
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE           :: forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_forces'

      CHARACTER(LEN=default_string_length) :: description
      INTEGER                                            :: handle, i, log_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_result_type), POINTER                  :: results
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(particle_list_type), POINTER                  :: particles
      !TYPE(qs_environment_type), POINTER                 :: qs_env

      ! TODO: do I need these?
      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      log_unit = cp_logger_get_default_io_unit()

      NULLIFY (cp_subsys)
      IF (ALLOCATED(forces)) THEN
         DEALLOCATE(forces)
      END IF

      CALL force_env_get(force_env, subsys=cp_subsys)
      !CALL force_env_get(force_env, qs_env=qs_env, subsys=cp_subsys)
      !CPASSERT(ASSOCIATED(qs_env))


      !CALL get_qs_env(qs_env, results=results)

      CALL cp_subsys_get(cp_subsys, particles=particles, results=results)
      CPASSERT(ASSOCIATED(particles))

      description = "[DIPOLE]"
      CALL get_results(results=results, description=description, values=dipole(1:3))
      dipole(1:3) = dipole(1:3) * debye

      !WRITE *, * 'get_forces: particles%n_els = ', particles%n_els
      !WRITE *, * 'get_forces: particles(1, 1:3) = ', particles%els(1)%f(1:3)

      ALLOCATE(forces(particles%n_els, 1:3))
      !forces(:, 1:3) = particles%els(:)%f(1:3)
      DO i = 1, particles%n_els
         forces(i, 1:3) = particles%els(i)%f(1:3)
      END DO


      CALL timestop(handle)

   END SUBROUTINE get_forces


   SUBROUTINE do_born(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'do_negf'

      INTEGER                                            :: handle, log_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      !TYPE(qs_environment_type), POINTER                 :: qs_env



      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      log_unit = cp_logger_get_default_io_unit()

      NULLIFY (cp_subsys)

      CALL force_env_get(force_env, subsys=cp_subsys)
      !CALL force_env_get(force_env, qs_env=qs_env, subsys=cp_subsys)


      CALL timestop(handle)
   END SUBROUTINE do_born
END MODULE

MODULE born_types
   USE cp_result_methods, ONLY: get_results
   USE cp_result_types, ONLY: cp_result_type
   USE cp_subsys_types, ONLY: cp_subsys_get, cp_subsys_type
   USE force_env_types, ONLY: force_env_get, force_env_type
   USE kinds, ONLY: dp, default_string_length
   USE particle_list_types, ONLY: particle_list_type
   USE physcon, ONLY: debye



#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'born_types'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   PUBLIC :: born_sample_type, born_samples_type, from_force_env

   ! TODO: different field strenghts for different axes?
   TYPE born_sample_type
      REAL(kind=dp) :: e_pot
      REAL(kind=dp), DIMENSION(1:3) :: dipole
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: forces
   END TYPE born_sample_type

   TYPE born_samples_type
      TYPE(born_sample_type) :: sample_no_field
      TYPE(born_sample_type), DIMENSION(1:3, 1:2) :: samples_field
      REAL(kind=dp) :: field_strength
   END TYPE born_samples_type


CONTAINS
   ! TODO: don't extract forces for no field
   SUBROUTINE from_force_env(force_env, born_sample)
      TYPE(force_env_type), POINTER :: force_env
      TYPE(born_sample_type) :: born_sample

      INTEGER :: i
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      CHARACTER(LEN=default_string_length) :: description
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(cp_result_type), POINTER                  :: results

      !ASSERT(ASSOCIATED(born_sample))

      !ALLOCATE(born_sample)

      !NULLIFY(cp_subsys)
      CALL force_env_get(force_env, subsys=cp_subsys, potential_energy=born_sample%e_pot)
      !CALL force_env_get(force_env, qs_env=qs_env, subsys=cp_subsys)
      !CPASSERT(ASSOCIATED(qs_env))


      !CALL get_qs_env(qs_env, results=results)

      CALL cp_subsys_get(cp_subsys, particles=particles, results=results)
      CPASSERT(ASSOCIATED(particles))


      ALLOCATE(born_sample%forces(particles%n_els, 1:3))
      !forces(:, 1:3) = particles%els(:)%f(1:3)
      DO i = 1, particles%n_els
         born_sample%forces(i, 1:3) = particles%els(i)%f(1:3)
      END DO

      description = "[DIPOLE]"
      CALL get_results(results=results, description=description, values=born_sample%dipole(1:3))
      !CALL get_results(results=results, description="[DIPOLE]", values=born_sample%dipole(1:3))
      born_sample%dipole(1:3) = born_sample%dipole(1:3) * debye



   END SUBROUTINE from_force_env



END MODULE born_types

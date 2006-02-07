!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/wiener_process [1.0] *
!!
!!   NAME
!!     wiener_process
!!
!!   FUNCTION
!!     Handling of the Wiener process currently employed in turn of the
!!     Langevin dynamics.
!!
!!   AUTHOR
!!     Matthias Krack (05.07.2005)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE wiener_process

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type,&
                                             init_local_particle_set
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type
  USE input_constants,                 ONLY: langevin_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             next_rng_seed
  USE particle_list_types,             ONLY: particle_list_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters in this module

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wiener_process'

  ! Public subroutines

  PUBLIC :: create_wiener_process

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE create_wiener_process(md_env,subsys_section,error)

    ! Purpose: Create a Wiener process for Langevin dynamics and initialize an
    !          independent random number generator for each atom in all force
    !          environment and all the subsystems/fragments therein.

    ! History: - Creation (06.07.2005,MK)

    ! *************************************************************************

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_wiener_process', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: name
    INTEGER :: i, iparticle, iparticle_kind, iparticle_local, istat, isubsys, &
      nparticle, nparticle_kind, nparticle_local
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: seed
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: work_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

! -------------------------------------------------------------------------

    failure = .FALSE.
    NULLIFY(work_section)
    CPPrecondition (ASSOCIATED(md_env),cp_failure_level,routineP,error,failure)

    CALL get_md_env(md_env=md_env,&
                    force_env_p=force_env_p,&
                    para_env=para_env,&
                    simpar=simpar)

    IF (simpar%ensemble == langevin_ensemble) THEN

      ! Load initial seed (not needed for a restart)

      initial_seed = next_rng_seed()

      DO i=1,SIZE(force_env_p)

        CALL force_env_get(force_env_p(i)%force_env,subsys=subsys)

        DO isubsys=1,SIZE(subsys)

          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             local_particles=local_particles,&
                             particles=particles,&
                             error=error)

          nparticle_kind = atomic_kinds%n_els
          nparticle = particles%n_els

          ! Allocate the (local) data structures for the Wiener process

          ALLOCATE(local_particles%local_particle_set(nparticle_kind),&
                   STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          DO iparticle_kind=1,nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            ALLOCATE (local_particles%local_particle_set(iparticle_kind)%&
                      &rng(nparticle_local),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            DO iparticle_local=1,nparticle_local
              NULLIFY (local_particles%local_particle_set(iparticle_kind)%&
                       &rng(iparticle_local)%stream)
            END DO
          END DO

          ! Each process generates all seeds. The seed generation should be
          ! quite fast and in this way a broadcast is avoided.

          ALLOCATE (seed(3,2,nparticle),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          seed(:,:,1) = initial_seed
          DO iparticle=2,nparticle
            seed(:,:,iparticle) = next_rng_seed(seed(:,:,iparticle-1))
          END DO

          ! Update initial seed

          initial_seed = next_rng_seed(seed(:,:,nparticle))

          ! Create a random number stream (Wiener process) for each particle

          DO iparticle_kind=1,nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%&
                          &array(iparticle_local)
              WRITE (UNIT=name,FMT="(A,I6)") "RNG for particle",iparticle
              CALL create_rng_stream(local_particles%&
                                     &local_particle_set(iparticle_kind)%&
                                     &rng(iparticle_local)%stream,&
                                     name=name,&
                                     extended_precision=.TRUE.,&
                                     seed=seed(:,:,iparticle))
            END DO
          END DO

          DEALLOCATE (seed,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          !
          ! Possibly restart rng
          !
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(subsys_section,"RNG_INIT",error=error)
          CALL init_local_particle_set(local_particles, atomic_kinds%n_els, particles%n_els,&
               work_section=work_section, error=error)

        END DO

      END DO

    END IF

  END SUBROUTINE create_wiener_process

END MODULE wiener_process

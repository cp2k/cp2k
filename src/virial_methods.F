!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH [04042007] code refactoring
! *****************************************************************************
MODULE virial_methods

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: virial_evaluate, virial_pair_force, virial_update

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_methods'

CONTAINS
! *****************************************************************************
!> \brief Updates the virial given the virial and subsys
!> \param virial ...
!> \param subsys ...
!> \param para_env ...
!> \param error ...
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE virial_update(virial, subsys, para_env, error)
    TYPE(virial_type), INTENT(INOUT)         :: virial
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_update', &
      routineP = moduleN//':'//routineN

    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles

    CALL cp_subsys_get(subsys, local_particles=local_particles, atomic_kinds=atomic_kinds,&
         particles=particles, error=error)

    CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles,&
         virial, para_env%group, error)

  END SUBROUTINE virial_update

! *****************************************************************************
!> \brief Computes the kinetic part of the pressure tensor and updates
!>      the full VIRIAL (PV)
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param virial ...
!> \param igroup ...
!> \param error ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE virial_evaluate ( atomic_kind_set, particle_set, local_particles,&
       virial, igroup, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: igroup
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      nnodes, nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: mass, mfl
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    IF ( virial%pv_availability ) THEN
       CALL timeset(routineN,handle)
       NULLIFY(atomic_kind)
       mfl = 0.0_dp
       nparticle_kind = SIZE ( atomic_kind_set )
       virial%pv_kinetic = 0.0_dp
       DO i = 1, 3
          DO j = 1, i
             nnodes = 0
             DO iparticle_kind=1,nparticle_kind
                atomic_kind => atomic_kind_set(iparticle_kind)
                CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass )
                nparticle_local = local_particles%n_el(iparticle_kind)
                DO iparticle_local=1,nparticle_local
                   nnodes = nnodes + 1
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   virial%pv_kinetic(i,j) = virial%pv_kinetic(i,j) + &
                        mass * particle_set(iparticle)%v(i)*particle_set(iparticle)%v(j)
                END DO
             END DO
             virial%pv_kinetic(j,i) = virial%pv_kinetic(i,j)
          END DO
       END DO
       mfl = REAL( 9 * nnodes, KIND=dp) * 2.0_dp * 1.e-6_dp

       CALL mp_sum(virial%pv_kinetic,igroup)

       ! total virial
       virial%pv_total = virial%pv_virial + virial%pv_kinetic + virial%pv_constraint

       CALL timestop(handle)
    ENDIF

  END SUBROUTINE virial_evaluate

! *****************************************************************************
!> \brief Computes the contribution to the stress tensor from two-body
!>      pair-wise forces
!> \param pv_virial ...
!> \param f0 ...
!> \param force ...
!> \param rab ...
!> \param error ...
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE virial_pair_force ( pv_virial, f0, force, rab, error)

    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_virial
    REAL(KIND=dp)                            :: f0
    REAL(KIND=dp), DIMENSION(3)              :: force, rab
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_pair_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

    DO i=1,3
       DO j=1,3
          pv_virial(i,j) = pv_virial(i,j) + f0 * force(i) * rab(j)
       END DO
    END DO

  END SUBROUTINE virial_pair_force

END MODULE virial_methods


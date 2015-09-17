!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Basic container type for QM/MM.
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmm_types
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_env_release,&
                                             fist_environment_type
  USE kinds,                           ONLY: dp
  USE qmmm_types_low,                  ONLY: qmmm_env_qm_release,&
                                             qmmm_env_qm_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_types'

  PUBLIC :: qmmm_env_type, qmmm_env_retain, qmmm_env_release, qmmm_env_get

  TYPE qmmm_env_type
     INTEGER                                                 :: ref_count = 1
     TYPE(qs_environment_type), POINTER                      :: qs_env => Null()
     TYPE(fist_environment_type), POINTER                    :: fist_env => Null()
     TYPE(qmmm_env_qm_type), POINTER                         :: qm => Null()
  END TYPE qmmm_env_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qmmm_env_get(qmmm_env,subsys,potential_energy,kinetic_energy,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(qs_energy_type), POINTER            :: qs_energy

    failure=.FALSE.
    NULLIFY(qs_energy, thermo)

    CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(kinetic_energy)) THEN
        CALL fist_env_get(qmmm_env%fist_env,thermo=thermo,error=error)
        kinetic_energy = thermo%kin
    END IF
    IF (PRESENT(subsys)) THEN
        CALL fist_env_get(qmmm_env%fist_env,subsys=subsys,error=error)
    ENDIF
    IF (PRESENT(potential_energy)) THEN
         ! get the underlying energies from primary subsys.  This is the only subsys
         ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
         CALL fist_env_get(qmmm_env%fist_env, thermo=thermo, error=error)
         CALL get_qs_env(qmmm_env%qs_env,energy=qs_energy,error=error)
         potential_energy = thermo%pot + qs_energy%total
    ENDIF
  END SUBROUTINE qmmm_env_get


! *****************************************************************************
!> \brief releases the given qmmm_env (see doc/ReferenceCounting.html)
!> \param qmmm_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_release(qmmm_env,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(qmmm_env)) THEN
       CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
       qmmm_env%ref_count=qmmm_env%ref_count-1
       IF (qmmm_env%ref_count==0) THEN
          CALL qs_env_release(qmmm_env%qs_env, error)
          CALL fist_env_release(qmmm_env%fist_env, error)
          CALL qmmm_env_qm_release(qmmm_env%qm, error)
          DEALLOCATE(qmmm_env)
       END IF
    END IF
    NULLIFY(qmmm_env)
  END SUBROUTINE qmmm_env_release

! *****************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_retain(qmmm_env,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
    qmmm_env%ref_count=qmmm_env%ref_count+1
  END SUBROUTINE qmmm_env_retain

END MODULE qmmm_types

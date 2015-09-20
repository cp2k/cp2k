!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Basic container type for QM/MM with force mixing.
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_types
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE kinds,                           ONLY: dp
  USE qmmm_types,                      ONLY: qmmm_env_get,&
                                             qmmm_env_release,&
                                             qmmm_env_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_types'

  PUBLIC :: qmmmx_env_type, qmmmx_env_get, qmmmx_env_retain, qmmmx_env_release

  TYPE qmmmx_env_type
     INTEGER                                                 :: ref_count = 1
     TYPE(qmmm_env_type), POINTER                            :: core => Null()
     TYPE(qmmm_env_type), POINTER                            :: ext => Null()
  END TYPE qmmmx_env_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
! *****************************************************************************
  SUBROUTINE qmmmx_env_get(qmmmx_env,subsys,potential_energy,kinetic_energy)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(qmmmx_env))
    CPASSERT(qmmmx_env%ref_count>0)

    ! get the underlying energies from primary subsys.  This is the only subsys
    ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
    CALL qmmm_env_get(qmmmx_env%ext,&
                      kinetic_energy=kinetic_energy,&
                      potential_energy=potential_energy,&
                      subsys=subsys)

  END SUBROUTINE qmmmx_env_get


! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmmx_env_retain(qmmmx_env)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(qmmmx_env))
    CPASSERT(qmmmx_env%ref_count>0)
    qmmmx_env%ref_count = qmmmx_env%ref_count+1
  END SUBROUTINE qmmmx_env_retain


! *****************************************************************************
!> \brief releases the given qmmmx_env (see doc/ReferenceCounting.html)
!> \param qmmmx_env the object to release
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmmx_env_release(qmmmx_env)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(qmmmx_env)) THEN
       CPASSERT(qmmmx_env%ref_count>0)
       qmmmx_env%ref_count = qmmmx_env%ref_count-1
       IF (qmmmx_env%ref_count==0) THEN
          CALL qmmm_env_release(qmmmx_env%core)
          CALL qmmm_env_release(qmmmx_env%ext)
          DEALLOCATE(qmmmx_env)
       END IF
    END IF
    NULLIFY(qmmmx_env)
  END SUBROUTINE qmmmx_env_release

END MODULE qmmmx_types

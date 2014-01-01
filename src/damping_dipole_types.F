!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Rodolphe Vuilleumier (29.12.2009)
! *****************************************************************************
MODULE damping_dipole_types

  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'damping_dipole_types'

! *** Global public parameters

  INTEGER, PUBLIC, PARAMETER :: no_damping=-1,&
                                tang_toennies=1

! *** Define the damping types ***
! *****************************************************************************
  TYPE damping_info_type
    CHARACTER (LEN=default_string_length)   :: atm_name1,atm_name2
    CHARACTER (LEN=default_string_length)   :: dtype
    INTEGER                                 :: order
    REAL(KIND=dp)                           :: bij,cij
  END TYPE damping_info_type
! *****************************************************************************
  TYPE damping_type
    INTEGER :: itype
    INTEGER :: order
    REAL(KIND=dp) :: bij,cij
  END TYPE damping_type

  TYPE damping_p_type
    TYPE(damping_type), DIMENSION(:), POINTER :: damp
  END TYPE

! *****************************************************************************

! *** Public data types ***

  PUBLIC :: damping_info_type, damping_type

! *** Public subroutines ***

  PUBLIC :: damping_p_type, damping_p_create, damping_p_release

CONTAINS

! *****************************************************************************
!> \brief Creates Data-structure that contains damping information
!> \author Rodolphe Vuilleumier
! *****************************************************************************
  SUBROUTINE damping_p_create(damping,nkinds,error)
    TYPE(damping_p_type), POINTER            :: damping
    INTEGER, INTENT(IN)                      :: nkinds
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'damping_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

     failure=.FALSE.
     IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(damping),cp_failure_level,routineP,error,failure)
       ALLOCATE ( damping, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( damping%damp( nkinds ), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, nkinds
          CALL init_damping ( damping%damp(i) )
       END DO

     END IF

  END SUBROUTINE damping_p_create
! *****************************************************************************
!> \brief Release Data-structure that contains damping information
!> \author Rodolphe Vuilleumier [RV]
! *****************************************************************************
  SUBROUTINE damping_p_release(damping,error)
    TYPE(damping_p_type), POINTER            :: damping
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'damping_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(ASSOCIATED(damping)) THEN
      IF (ASSOCIATED(damping%damp)) THEN
         DEALLOCATE(damping%damp,stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
      DEALLOCATE(damping,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    NULLIFY(damping)

  END SUBROUTINE damping_p_release

! *****************************************************************************
  SUBROUTINE init_damping(damping)
    TYPE(damping_type)                       :: damping

    damping%itype=no_damping
    damping%order=1
    damping%bij=HUGE(0.0_dp)
    damping%cij=0.0_dp

  END SUBROUTINE init_damping

! *****************************************************************************
END MODULE damping_dipole_types

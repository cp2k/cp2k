!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Rodolphe Vuilleumier (29.12.2009)
! *****************************************************************************
MODULE damping_dipole_types

  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "../common/cp_common_uses.f90"

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
!> \param damping ...
!> \param nkinds ...
!> \author Rodolphe Vuilleumier
! *****************************************************************************
  SUBROUTINE damping_p_create(damping,nkinds)
    TYPE(damping_p_type), POINTER            :: damping
    INTEGER, INTENT(IN)                      :: nkinds

    CHARACTER(LEN=*), PARAMETER :: routineN = 'damping_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

     failure=.FALSE.
     CPASSERT(.NOT.ASSOCIATED(damping))
     ALLOCATE ( damping)
     ALLOCATE ( damping%damp( nkinds ))
     DO i = 1, nkinds
        CALL init_damping ( damping%damp(i) )
     END DO


  END SUBROUTINE damping_p_create
! *****************************************************************************
!> \brief Release Data-structure that contains damping information
!> \param damping ...
!> \author Rodolphe Vuilleumier [RV]
! *****************************************************************************
  SUBROUTINE damping_p_release(damping)
    TYPE(damping_p_type), POINTER            :: damping

    CHARACTER(len=*), PARAMETER :: routineN = 'damping_p_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(ASSOCIATED(damping)) THEN
      IF (ASSOCIATED(damping%damp)) THEN
         DEALLOCATE(damping%damp)
      END IF
      DEALLOCATE(damping)
    END IF
    NULLIFY(damping)

  END SUBROUTINE damping_p_release

! *****************************************************************************
!> \brief ...
!> \param damping ...
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

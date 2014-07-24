!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief an exclusion type
! *****************************************************************************
MODULE exclusion_types

#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \brief A type used to store lists of exclusions and onfos
!> \par History
!>      12.2010 created [Joost VandeVondele]
! *****************************************************************************
  TYPE exclusion_type
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_vdw
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_ei
     INTEGER, POINTER, DIMENSION(:)        :: list_onfo
  END TYPE

  PUBLIC :: exclusion_type,&
            exclusion_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'exclusion_types'

CONTAINS

! *****************************************************************************
!> \brief Release exclusion type
!> \param exclusions ...
!> \param error ...
!> \par History
!>      12.2010 created [Teodoro Laino] - teodoro.laino@gmail.com
!> \author teo
! *****************************************************************************
 SUBROUTINE exclusion_release(exclusions, error)
    TYPE(exclusion_type), DIMENSION(:), &
      POINTER                                :: exclusions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exclusion_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, stat
    LOGICAL                                  :: failure

   failure=.FALSE.
   IF (ASSOCIATED(exclusions)) THEN
      DO iatom=1,SIZE(exclusions)
         IF (ASSOCIATED(exclusions(iatom)%list_exclude_vdw,&
              exclusions(iatom)%list_exclude_ei)) THEN
            DEALLOCATE(exclusions(iatom)%list_exclude_vdw,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ELSE
            IF (ASSOCIATED(exclusions(iatom)%list_exclude_vdw)) THEN
               DEALLOCATE(exclusions(iatom)%list_exclude_vdw,STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END IF
            IF (ASSOCIATED(exclusions(iatom)%list_exclude_ei)) THEN
               DEALLOCATE(exclusions(iatom)%list_exclude_ei,STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END IF
         ENDIF
         IF (ASSOCIATED(exclusions(iatom)%list_onfo)) THEN
           DEALLOCATE(exclusions(iatom)%list_onfo,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ENDIF
      ENDDO
      DEALLOCATE(exclusions,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   ENDIF
 END SUBROUTINE exclusion_release

END MODULE exclusion_types

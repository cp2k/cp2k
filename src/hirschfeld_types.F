!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief The types needed for the calculation of Hirschfeld charges and
!>        related functions
!> \par History
!>      11.2014 created [JGH]
!> \author JGH
! *****************************************************************************
MODULE hirschfeld_types
  USE f77_blas
  USE input_constants,                 ONLY: shape_function_gaussian
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hirschfeld_types'

  PUBLIC :: hirschfeld_type
  PUBLIC :: create_hirschfeld_type, release_hirschfeld_type
  PUBLIC :: get_hirschfeld_info, set_hirschfeld_info

! *****************************************************************************
!> \brief quantities needed for a Hischfeld based partitioning of real space
!> \author JGH
! *****************************************************************************
  TYPE hirschfeld_type
     LOGICAL                       :: iterative
     INTEGER                       :: shape_function_type
     INTEGER                       :: ref_charge
     TYPE(shape_fn),DIMENSION(:),&
        POINTER                    :: kind_shape_fn
     REAL(KIND=dp),DIMENSION(:),&
        POINTER                    :: charges
     TYPE(pw_p_type), POINTER      :: fnorm
  END TYPE hirschfeld_type

  TYPE shape_fn
     INTEGER                       :: numexp
     REAL(KIND=dp),DIMENSION(:),&
        POINTER                    :: zet
     REAL(KIND=dp),DIMENSION(:),&
        POINTER                    :: coef
  END TYPE shape_fn

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param hirschfeld_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_hirschfeld_type(hirschfeld_env,error)
    TYPE(hirschfeld_type), POINTER           :: hirschfeld_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hirschfeld_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(hirschfeld_env)) THEN
       CALL release_hirschfeld_type(hirschfeld_env,error)
    END IF

    ALLOCATE(hirschfeld_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    hirschfeld_env%iterative = .FALSE.
    hirschfeld_env%shape_function_type = shape_function_gaussian
    NULLIFY(hirschfeld_env%kind_shape_fn)
    NULLIFY(hirschfeld_env%charges)
    NULLIFY(hirschfeld_env%fnorm)

  END SUBROUTINE create_hirschfeld_type

! *****************************************************************************
!> \brief ...
!> \param hirschfeld_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_hirschfeld_type(hirschfeld_env,error)
    TYPE(hirschfeld_type), POINTER           :: hirschfeld_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_hirschfeld_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure
    TYPE(shape_fn), DIMENSION(:), POINTER    :: kind_shape

    failure = .FALSE.

    IF(ASSOCIATED(hirschfeld_env)) THEN

       IF(ASSOCIATED(hirschfeld_env%kind_shape_fn)) THEN
          kind_shape => hirschfeld_env%kind_shape_fn
          DO ikind=1,SIZE(kind_shape)
             IF(ASSOCIATED(hirschfeld_env%kind_shape_fn(ikind)%zet)) THEN
                DEALLOCATE(kind_shape(ikind)%zet, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF(ASSOCIATED(hirschfeld_env%kind_shape_fn(ikind)%coef)) THEN
                DEALLOCATE(kind_shape(ikind)%coef, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
          DEALLOCATE(kind_shape, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF(ASSOCIATED(hirschfeld_env%charges)) THEN
          DEALLOCATE(hirschfeld_env%charges, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(hirschfeld_env%fnorm)) THEN
          CALL pw_release(hirschfeld_env%fnorm%pw,error=error)
          DEALLOCATE(hirschfeld_env%fnorm,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF

       DEALLOCATE(hirschfeld_env, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

  END SUBROUTINE release_hirschfeld_type

! *****************************************************************************
!> \brief ...
!> \param hirschfeld_env ...
!> \param shape_function_type ...
!> \param iterative ...
!> \param ref_charge ...
!> \param fnorm ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_hirschfeld_info(hirschfeld_env,shape_function_type,iterative,&
             ref_charge,fnorm,error)
    TYPE(hirschfeld_type), POINTER           :: hirschfeld_env
    INTEGER, INTENT(OUT), OPTIONAL           :: shape_function_type
    LOGICAL, INTENT(OUT), OPTIONAL           :: iterative
    INTEGER, INTENT(OUT), OPTIONAL           :: ref_charge
    TYPE(pw_p_type), OPTIONAL, POINTER       :: fnorm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_hirschfeld_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(hirschfeld_env),cp_failure_level,routineP,error,failure)

    IF(PRESENT(shape_function_type)) THEN
       shape_function_type = hirschfeld_env%shape_function_type
    END IF
    IF(PRESENT(iterative)) THEN
       iterative = hirschfeld_env%iterative
    END IF
    IF(PRESENT(ref_charge)) THEN
       ref_charge = hirschfeld_env%ref_charge
    END IF
    IF(PRESENT(fnorm)) THEN
       fnorm => hirschfeld_env%fnorm
    END IF

  END SUBROUTINE get_hirschfeld_info

! *****************************************************************************
!> \brief ...
!> \param hirschfeld_env ...
!> \param shape_function_type ...
!> \param iterative ...
!> \param ref_charge ...
!> \param fnorm ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_hirschfeld_info(hirschfeld_env,shape_function_type,iterative,&
             ref_charge,fnorm,error)
    TYPE(hirschfeld_type), POINTER           :: hirschfeld_env
    INTEGER, INTENT(IN), OPTIONAL            :: shape_function_type
    LOGICAL, INTENT(IN), OPTIONAL            :: iterative
    INTEGER, INTENT(IN), OPTIONAL            :: ref_charge
    TYPE(pw_p_type), OPTIONAL, POINTER       :: fnorm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_hirschfeld_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(hirschfeld_env),cp_failure_level,routineP,error,failure)

    IF(PRESENT(shape_function_type)) THEN
       hirschfeld_env%shape_function_type = shape_function_type
    END IF
    IF(PRESENT(iterative)) THEN
       hirschfeld_env%iterative = iterative
    END IF
    IF(PRESENT(ref_charge)) THEN
       hirschfeld_env%ref_charge = ref_charge
    END IF
    IF(PRESENT(fnorm)) THEN
       hirschfeld_env%fnorm => fnorm
    END IF

  END SUBROUTINE set_hirschfeld_info
! *****************************************************************************

END MODULE hirschfeld_types

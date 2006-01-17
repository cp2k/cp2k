!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_environment_types [1.0] *
!!
!!   NAME
!!     tbmd_environment_types
!!
!!   FUNCTION
!!     The tight-binding environment type, which contains everything needed 
!!     to actually do tight-binding.
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_environment_types

#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='tbmd_environment_types'
  INTEGER, PRIVATE, save :: last_tbmd_id=0

  ! *** Public data types ***
  PUBLIC :: tbmd_environment_type, &
            tbmd_environment_p_type

  ! *** Public subroutines ***
  PUBLIC :: tbmd_env_retain, &
            tbmd_env_release, &
            tbmd_env_set,&
            tbmd_env_get, &
            tbmd_env_clear, &
            tbmd_env_create
!!***
!****************************************************************************

!!****s* tbmd_environment_types/tbmd_environment_type [1.0] *
!!
!!   NAME
!!     tbmd_environment_type
!!
!!   FUNCTION
!!     The TBMD Environment
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE tbmd_environment_type
    INTEGER                                          :: id_nr, &
                                                        ref_count

  END TYPE tbmd_environment_type
!!***
!****************************************************************************

!!****s* tbmd_environment_types/tbmd_environment_p_type [1.0] *
!!
!!   NAME
!!     tbmd_environment_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - tbmd_env_ptr: Pointer to the actual tbmd environment
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE tbmd_environment_p_type
    TYPE(tbmd_environment_type), POINTER :: tbmd_env_ptr
  END TYPE tbmd_environment_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* tbmd_environment_types/tbmd_env_retain [1.0] *
!!
!!   NAME
!!     tbmd_env_retain
!!
!!   FUNCTION
!!     Retains a tbmd environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_retain(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------
   
    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_env%ref_count = tbmd_env%ref_count+1
    END IF
  END SUBROUTINE tbmd_env_retain
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_release [1.0] *
!!
!!   NAME
!!     tbmd_env_release
!!
!!   FUNCTION
!!     Releases the given tbmd environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_release(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(tbmd_env)) THEN
      CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_env%ref_count = tbmd_env%ref_count-1
      IF (tbmd_env%ref_count<1) THEN
        !IF (ASSOCIATED(tbmd_env%tbmd_force)) THEN
        !  CALL tbmd_force_release(tbmd_env%tbmd_force, error=error)
        !END IF

        DEALLOCATE(tbmd_env, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(tbmd_env)
  END SUBROUTINE tbmd_env_release
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_get [1.0] *
!!
!!   NAME
!!     tbmd_env_get
!!
!!   FUNCTION
!!     Returns various attributes of the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The enquired tbmd environment
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_get(tbmd_env, id_nr, error)

    TYPE(tbmd_environment_type), POINTER         :: tbmd_env
    INTEGER, INTENT(OUT), OPTIONAL               :: id_nr
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                      :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = tbmd_env%id_nr
      !IF (PRESENT(tbmd_force)) tbmd_force => tbmd_env%tbmd_force
    END IF

  END SUBROUTINE tbmd_env_get
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_set [1.0] *
!!
!!   NAME
!!     tbmd_env_set
!!
!!   FUNCTION
!!     Sets various attributes of the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The enquired tbmd environment
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_set(tbmd_env, id_nr, error)

    TYPE(tbmd_environment_type), POINTER        :: tbmd_env
    INTEGER, INTENT(IN), OPTIONAL               :: id_nr
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                     :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) tbmd_env%id_nr = id_nr
      !IF (PRESENT(tbmd_force)) tbmd_env%tbmd_force => tbmd_force
    END IF
  END SUBROUTINE tbmd_env_set
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_clear [1.0] *
!!
!!   NAME
!!     tbmd_env_clear
!!
!!   FUNCTION
!!     Reinitializes the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The tbmd environment to be reinitialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_clear(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      !CALL tbmd_force_clear(tbmd_env%tbmd_force, error=error)
     ELSE
      CALL tbmd_env_create(tbmd_env, error=error)
    END IF
  END SUBROUTINE tbmd_env_clear
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_create [1.0] *
!!
!!   NAME
!!     tbmd_env_create
!!
!!   FUNCTION
!!     Creates the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The tbmd environment to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_create(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(tbmd_env, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
   
    IF (.NOT. failure) THEN
      tbmd_env%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_env%id_nr = last_tbmd_id
    END IF
  END SUBROUTINE tbmd_env_create
!***************************************************************************

END MODULE tbmd_environment_types

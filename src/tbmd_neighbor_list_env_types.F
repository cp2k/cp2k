!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_neighbor_list_env_types [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_types
!!
!!   FUNCTION
!!     The TBMD neighbor-list environment type
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
MODULE tbmd_neighbor_list_env_types
  USE fist_neighbor_list_types,           ONLY: neighbor_list_set_p_type, &
                                                deallocate_neighbor_list_set, &
                                                neighbor_list_set_type
  USE kinds,                              ONLY: dp, &
                                                default_string_length

#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "tbmd_neighbor_list_env_types"
  INTEGER, PRIVATE, save :: last_tbmd_id=0

  ! *** Public data types ***
  PUBLIC :: tbmd_neighbor_list_env_type, &
            tbmd_neighbor_list_env_p_type

! *** Public subroutines ***
  PUBLIC :: tbmd_neighbor_list_env_retain, &
            tbmd_neighbor_list_env_release, &
            tbmd_neighbor_list_env_get, &
            tbmd_neighbor_list_env_set, &
            tbmd_neighbor_list_env_clear, &
            tbmd_neighbor_list_env_create

!!***
! *****************************************************************************

!!****s* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_type [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_type
!!
!!   FUNCTION
!!     The TBMD neighbor-list environment
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - unit_type:
!!     - counter:
!!     - last_update:
!!     - num_update:
!!     - verlet_skin:
!!     - lup:
!!     - aup:
!!     - r_last_update:
!!     - rlist_cutsq:
!!     - tbmd_neighbor_list:
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE tbmd_neighbor_list_env_type
    INTEGER                                          :: id_nr, &
                                                        ref_count
    CHARACTER(len = default_string_length)           :: unit_type
    INTEGER                                          :: counter, &
                                                        last_update, &
                                                        num_update
    REAL(KIND=dp)                                    :: verlet_skin, &
                                                        lup, &
                                                        aup
    REAL(KIND=dp), DIMENSION(:,:), POINTER           :: r_last_update, &
                                                        rlist_cutsq
    TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                     :: tbmd_neighbor_list
  END TYPE tbmd_neighbor_list_env_type
!!***
!****************************************************************************

!!****s* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_p_type [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - tbmd_neighbor_list_env_ptr: Pointer to the actual tbmd neighbor-list 
!!                                   environment
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE tbmd_neighbor_list_env_p_type
    TYPE(tbmd_neighbor_list_env_type), POINTER :: tbmd_neighbor_list_env_ptr
  END TYPE tbmd_neighbor_list_env_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_retain [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_retain
!!
!!   FUNCTION
!!     Retains a tbmd neighbor-list environment 
!!        (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_neighbor_list_env: The tbmd neighbor list environment to retain
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
  SUBROUTINE tbmd_neighbor_list_env_retain(tbmd_neighbor_list_env, error)

    TYPE(tbmd_neighbor_list_env_type), POINTER :: tbmd_neighbor_list_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_neighbor_list_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                    :: failure

!   ------------------------------------------------------------------------
   
    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_neighbor_list_env), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_neighbor_list_env%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_neighbor_list_env%ref_count = tbmd_neighbor_list_env%ref_count+1
    END IF
  END SUBROUTINE tbmd_neighbor_list_env_retain
!***************************************************************************

!!****f* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_release [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_release
!!
!!   FUNCTION
!!     Releases the given tbmd neighbor-list environment 
!!        (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_neighbor_list_env: The tbmd neighbor list environment to release
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
  SUBROUTINE tbmd_neighbor_list_env_release(tbmd_neighbor_list_env, error)

    TYPE(tbmd_neighbor_list_env_type), POINTER :: tbmd_neighbor_list_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_neighbor_list_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                    :: stat, &
                                                  i
    LOGICAL                                    :: failure
    TYPE(neighbor_list_set_type), POINTER      :: neighbor_list_set_ptr

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(tbmd_neighbor_list_env)) THEN
      CPPrecondition(tbmd_neighbor_list_env%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_neighbor_list_env%ref_count = tbmd_neighbor_list_env%ref_count-1
      IF (tbmd_neighbor_list_env%ref_count<1) THEN
        IF (ASSOCIATED(tbmd_neighbor_list_env%tbmd_neighbor_list)) THEN
          DO i=1,SIZE(tbmd_neighbor_list_env%tbmd_neighbor_list)
            neighbor_list_set_ptr => tbmd_neighbor_list_env%tbmd_neighbor_list(i)%neighbor_list_set
            CALL deallocate_neighbor_list_set(neighbor_list_set_ptr)
          END DO
        END IF
        DEALLOCATE(tbmd_neighbor_list_env, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(tbmd_neighbor_list_env)
  END SUBROUTINE tbmd_neighbor_list_env_release
!***************************************************************************

!!****f* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_get [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_get
!!
!!   FUNCTION
!!     Returns various attributes of the tbmd neighbor-list environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_neighbor_list_env: The enquired tbmd neighbor-list environment
!!     - id_nr: unique identifier
!!     - unit_type:
!!     - counter:
!!     - last_update:
!!     - num_update:
!!     - verlet_skin:
!!     - lup:
!!     - aup:
!!     - r_last_update:
!!     - rlist_cutsq:
!!     - tbmd_neighbor_list:
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_neighbor_list_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_neighbor_list_env_get(tbmd_neighbor_list_env, id_nr, &
                                        unit_type, counter, last_update, &
                                        num_update, verlet_skin, lup, aup, &
                                        r_last_update, rlist_cutsq, &
                                        tbmd_neighbor_list, error)

    TYPE(tbmd_neighbor_list_env_type), POINTER   :: tbmd_neighbor_list_env
    INTEGER, INTENT(OUT), OPTIONAL               :: id_nr
    CHARACTER(len = default_string_length), &
         INTENT(OUT), OPTIONAL                   :: unit_type
    INTEGER, INTENT(OUT), OPTIONAL               :: counter, &
                                                    last_update, &
                                                    num_update
    REAL(KIND=dp), INTENT(OUT), OPTIONAL         :: verlet_skin, &
                                                    lup, &
                                                    aup
    REAL(KIND=dp), DIMENSION(:,:), POINTER, &
         OPTIONAL                                :: r_last_update, &
                                                    rlist_cutsq
    TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER, OPTIONAL                       :: tbmd_neighbor_list
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_neighbor_list_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                      :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(tbmd_neighbor_list_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(tbmd_neighbor_list_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = tbmd_neighbor_list_env%id_nr
      IF (PRESENT(unit_type)) unit_type = tbmd_neighbor_list_env%unit_type
      IF (PRESENT(counter)) counter = tbmd_neighbor_list_env%counter
      IF (PRESENT(last_update)) THEN
        last_update = tbmd_neighbor_list_env%last_update
      END IF
      IF (PRESENT(num_update)) THEN
        num_update = tbmd_neighbor_list_env%num_update
      END IF
      IF (PRESENT(verlet_skin)) THEN
        verlet_skin = tbmd_neighbor_list_env%verlet_skin
      END IF
      IF (PRESENT(lup)) lup = tbmd_neighbor_list_env%lup
      IF (PRESENT(aup)) aup = tbmd_neighbor_list_env%aup
      IF (PRESENT(r_last_update)) THEN
        r_last_update => tbmd_neighbor_list_env%r_last_update
      END IF
      IF (PRESENT(rlist_cutsq)) THEN
        rlist_cutsq => tbmd_neighbor_list_env%rlist_cutsq
      END IF
      IF (PRESENT(tbmd_neighbor_list)) THEN
        tbmd_neighbor_list => tbmd_neighbor_list_env%tbmd_neighbor_list
      END IF
    END IF

  END SUBROUTINE tbmd_neighbor_list_env_get
!***************************************************************************

!!****f* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_set [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_set
!!
!!   FUNCTION
!!     Sets various attributes of the tbmd neighbor-list environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_neighbor_list_env: The enquired tbmd neighbor-list environment
!!     - id_nr: unique identifier
!!     - unit_type:
!!     - counter:
!!     - last_update:
!!     - num_update:
!!     - verlet_skin:
!!     - lup:
!!     - aup:
!!     - r_last_update:
!!     - rlist_cutsq:
!!     - tbmd_neighbor_list:
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_neighbor_list_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_neighbor_list_env_set(tbmd_neighbor_list_env, id_nr, &
                                        unit_type, counter, last_update, &
                                        num_update, verlet_skin, lup, aup, &
                                        r_last_update, rlist_cutsq, &
                                        tbmd_neighbor_list, error)

    TYPE(tbmd_neighbor_list_env_type), POINTER   :: tbmd_neighbor_list_env
    INTEGER, INTENT(OUT), OPTIONAL               :: id_nr
    CHARACTER(len = default_string_length), &
         INTENT(OUT), OPTIONAL                   :: unit_type
    INTEGER, INTENT(OUT), OPTIONAL               :: counter, &
                                                    last_update, &
                                                    num_update
    REAL(KIND=dp), INTENT(OUT), OPTIONAL         :: verlet_skin, &
                                                    lup, &
                                                    aup
    REAL(KIND=dp), DIMENSION(:,:), POINTER, &
         OPTIONAL                                :: r_last_update, &
                                                    rlist_cutsq
    TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER, OPTIONAL                       :: tbmd_neighbor_list
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_neighbor_list_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                      :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_neighbor_list_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(tbmd_neighbor_list_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) tbmd_neighbor_list_env%id_nr = id_nr
      IF (PRESENT(unit_type)) tbmd_neighbor_list_env%unit_type = unit_type
      IF (PRESENT(counter)) tbmd_neighbor_list_env%counter = counter
      IF (PRESENT(last_update)) THEN
        tbmd_neighbor_list_env%last_update = last_update
      END IF
      IF (PRESENT(num_update)) THEN
        tbmd_neighbor_list_env%num_update = num_update
      END IF
      IF (PRESENT(verlet_skin)) THEN
        tbmd_neighbor_list_env%verlet_skin = verlet_skin
      END IF
      IF (PRESENT(lup)) tbmd_neighbor_list_env%lup = lup
      IF (PRESENT(aup)) tbmd_neighbor_list_env%aup = aup
      IF (PRESENT(r_last_update)) THEN
        tbmd_neighbor_list_env%r_last_update => r_last_update
      END IF
      IF (PRESENT(rlist_cutsq)) THEN
        tbmd_neighbor_list_env%rlist_cutsq => rlist_cutsq
      END IF
      IF (PRESENT(tbmd_neighbor_list)) THEN
        tbmd_neighbor_list_env%tbmd_neighbor_list => tbmd_neighbor_list
      END IF
    END IF
  END SUBROUTINE tbmd_neighbor_list_env_set
!***************************************************************************

!!****f* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_clear [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_clear
!!
!!   FUNCTION
!!     Reinitializes the tbmd neighbor-list environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_neighbor_list_env: The tbmd neighbor- list environment to 
!!                               be reinitialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_neighbor_list_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_neighbor_list_env_clear(tbmd_neighbor_list_env, error)

    TYPE(tbmd_neighbor_list_env_type), POINTER :: tbmd_neighbor_list_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_neighbor_list_env_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                    :: failure
    INTEGER                                    :: i
    TYPE(neighbor_list_set_type), POINTER      :: neighbor_list_set_ptr    

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_neighbor_list_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      !id_nr = 0
      tbmd_neighbor_list_env%unit_type          = ""
      tbmd_neighbor_list_env%counter            = 0
      tbmd_neighbor_list_env%last_update        = 0
      tbmd_neighbor_list_env%num_update         = 0
      tbmd_neighbor_list_env%verlet_skin        = 0.0_dp
      tbmd_neighbor_list_env%lup                = 0.0_dp
      tbmd_neighbor_list_env%aup                = 0.0_dp
      tbmd_neighbor_list_env%r_last_update(:,:) = 0.0_dp
      tbmd_neighbor_list_env%rlist_cutsq(:,:)   = 0.0_dp
      IF (ASSOCIATED(tbmd_neighbor_list_env%tbmd_neighbor_list)) THEN
        DO i=1,SIZE(tbmd_neighbor_list_env%tbmd_neighbor_list)
            neighbor_list_set_ptr => tbmd_neighbor_list_env%tbmd_neighbor_list(i)%neighbor_list_set
          CALL deallocate_neighbor_list_set(neighbor_list_set_ptr)
        END DO
      END IF
     ELSE
      CALL tbmd_neighbor_list_env_create(tbmd_neighbor_list_env, error=error)
    END IF
  END SUBROUTINE tbmd_neighbor_list_env_clear
!***************************************************************************

!!****f* tbmd_neighbor_list_env_types/tbmd_neighbor_list_env_create [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_env_create
!!
!!   FUNCTION
!!     Creates the tbmd nighbor-list environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_neighbor_list_env: The tbmd neighbor-list environment to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_neighbor_list_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_neighbor_list_env_create(tbmd_neighbor_list_env, error)

    TYPE(tbmd_neighbor_list_env_type), POINTER :: tbmd_neighbor_list_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_neighbor_list_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                    :: stat
    LOGICAL                                    :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(tbmd_neighbor_list_env, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    NULLIFY(tbmd_neighbor_list_env%tbmd_neighbor_list)

    IF (.NOT. failure) THEN
      tbmd_neighbor_list_env%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_neighbor_list_env%id_nr = last_tbmd_id
      CALL tbmd_neighbor_list_env_clear(tbmd_neighbor_list_env, error=error)
    END IF
  END SUBROUTINE tbmd_neighbor_list_env_create
!***************************************************************************

END MODULE tbmd_neighbor_list_env_types

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_bas_projection_types [1.0] *
!!
!!   NAME
!!     pao_bas_projection_types
!!
!!   FUNCTION
!!     represent a simple projection (simply keeps some indexes)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 3.2002 created [fawzi]
!!     - 4.2002 added create, dealloc, write
!!
!!   SOURCE
!****************************************************************************
MODULE pao_bas_projection_types
  USE cp_array_i_utils,                ONLY: cp_1d_i_write
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_type
  USE cp_linked_list_int,              ONLY: cp_create,&
                                             cp_dealloc,&
                                             cp_insert,&
                                             cp_next,&
                                             cp_sll_int_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_bas_projection'
  
  PUBLIC :: pao_bas_projection_type
  PUBLIC :: cp_create, cp_dealloc, cp_init, cp_dealloc_ref,&
       cp_valid, cp_get, cp_get_proj_indexes, cp_get_non_proj_indexes,&
       cp_write
  PUBLIC :: cp_proj_create, cp_proj_dealloc, cp_proj_init,&
       cp_proj_dealloc_ref, cp_proj_valid, cp_proj_get,&
       cp_proj_get_proj_indexes, cp_proj_get_non_proj_indexes,&
       cp_proj_write , cp_proj_init2, cp_proj_create2

  INTERFACE cp_create
     MODULE PROCEDURE cp_proj_create , cp_proj_create2
  END INTERFACE
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_proj_dealloc
  END INTERFACE
  INTERFACE  cp_init
     MODULE PROCEDURE cp_proj_init , cp_proj_init2
  END INTERFACE
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_proj_dealloc_ref
  END INTERFACE
  INTERFACE cp_valid
     MODULE PROCEDURE cp_proj_valid
  END INTERFACE
  INTERFACE cp_get
     MODULE PROCEDURE cp_proj_get
  END INTERFACE
  INTERFACE cp_get_proj_indexes
     MODULE PROCEDURE cp_proj_get_proj_indexes
  END INTERFACE
  INTERFACE cp_get_non_proj_indexes
     MODULE PROCEDURE cp_proj_get_non_proj_indexes
  END INTERFACE
  INTERFACE cp_write
     MODULE PROCEDURE cp_proj_write
  END INTERFACE

!!***
!****************************************************************************

!!****s* pao_types/pao_bas_projection_type [1.0] *
!!
!!   NAME
!!     pao_bas_projection_type
!!
!!   FUNCTION
!!     represent a simple projection (ie a subset of indexes)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - initialized: true if initialized (bug catcher)
!!     - proj_indexes: the local indexes of the full basis that are kept
!!       by the projection
!!     - non_proj_indexes: the local indexes of the full basis that are
!!       not in the min_bas (i.e. the indexes of the excluded basis)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  TYPE pao_bas_projection_type
     PRIVATE
     LOGICAL :: initialized
     INTEGER, DIMENSION(:), POINTER :: proj_indexes
     INTEGER, DIMENSION(:), POINTER :: non_proj_indexes
  END TYPE pao_bas_projection_type
!!***
  !****************************************************************************

CONTAINS

!!****f* pao_bas_projection/cp_proj_create [1.0] *
!!
!!   NAME
!!     cp_proj_create
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_create(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Pointer:: proj
!!       Integer, Dimension(:), Intent (IN):: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_create
!!
!!   FUNCTION
!!     allocates and initianizes a pao_bas_projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     see cp_proj_init
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_proj_create(proj, proj_indexes, size_full_bas, &
     error)
  TYPE(pao_bas_projection_type), POINTER :: proj
  INTEGER, DIMENSION(:), INTENT(in) ::proj_indexes
  INTEGER, INTENT(in) :: size_full_bas
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: stat
  CHARACTER(len=*), PARAMETER :: routineN='cp_proj_create',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  ALLOCATE(proj,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.not.failure) THEN
     CALL cp_proj_init(proj, proj_indexes, size_full_bas, error)
  END IF
END SUBROUTINE cp_proj_create
!***************************************************************************

!!****f* pao_bas_projection/cp_proj_create2 [1.0] *
!!
!!   NAME
!!     cp_proj_create2
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_create2(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Pointer:: proj
!!       Type(cp_sll_int_type), Pointer:: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_create2
!!
!!   FUNCTION
!!     allocates and initianizes a pao_bas_projection (from a linked list)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     see cp_proj_init2
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_proj_create2(proj, proj_indexes, size_full_bas, &
     error)
  TYPE(pao_bas_projection_type), POINTER :: proj
  TYPE(cp_sll_int_type), POINTER ::proj_indexes
  INTEGER, INTENT(in) :: size_full_bas
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: stat
  CHARACTER(len=*), PARAMETER :: routineN='cp_proj_create2',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  ALLOCATE(proj,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.not.failure) THEN
     CALL cp_proj_init2(proj, proj_indexes, size_full_bas, error)
  END IF
END SUBROUTINE cp_proj_create2
!***************************************************************************

!!****f* pao_bas_projection/cp_proj_dealloc [1.0] *
!!
!!   NAME
!!     cp_proj_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_dealloc(proj, error)
!!       Type(pao_bas_projection_type), Pointer:: proj
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_dealloc
!!
!!   FUNCTION
!!     deallocates a projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to dealloc
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_proj_dealloc(proj,error)
  TYPE(pao_bas_projection_type), POINTER :: proj
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: stat
  CHARACTER(len=*), PARAMETER :: routineN='cp_proj_dealloc',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  IF (ASSOCIATED(proj)) THEN
     CALL cp_proj_dealloc_ref(proj,error=error)
     DEALLOCATE(proj,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  END IF
END SUBROUTINE cp_proj_dealloc
!***************************************************************************

!!****f* pao_bas_projection/cp_proj_init [1.0] *
!!
!!   NAME
!!     cp_proj_init
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_init(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Intent (OUT):: proj
!!       Integer, Dimension(:), Intent (IN):: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_init
!!
!!   FUNCTION
!!     initializes the given projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to be initialized
!!     - proj_indexes: the indexes of the full basis that are kept
!!     - size_full_bas: the size of the full basis
!!     - non_proj_indexes: the index of the excluded basis
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_proj_init(proj, proj_indexes, size_full_bas, &
       error)
    TYPE(pao_bas_projection_type), INTENT(out) :: proj
    INTEGER, DIMENSION(:), INTENT(in) ::proj_indexes
    INTEGER, INTENT(in) :: size_full_bas
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure,found
    INTEGER :: stat, i, j, exc_pos
    CHARACTER(len=*), PARAMETER :: routineN='cp_proj_init',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ALL(proj_indexes>0),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(proj_indexes)<=size_full_bas,cp_failure_level,routineP,error,failure)
    ALLOCATE(proj%proj_indexes(SIZE(proj_indexes)),&
         proj%non_proj_indexes(size_full_bas-SIZE(proj_indexes)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !do i=2,size(proj_indexes)
    !   CPPrecondition(proj_indexes(i-1)<proj_indexes(i),cp_failure_level,routineP,error,failure)
    !end do
    IF (.not.failure) THEN
       proj%proj_indexes=proj_indexes
       exc_pos=0
       DO i=1,size_full_bas
          found=.FALSE.
          DO j=1,SIZE(proj%proj_indexes)
             IF (i==proj%proj_indexes(j)) THEN
                found=.TRUE.
                EXIT
             END IF
          END DO
          IF (.not.found) THEN
             CPInvariant(exc_pos<SIZE(proj%non_proj_indexes),cp_failure_level,routineP,error,failure)
             IF (.not.failure) THEN
                exc_pos=exc_pos+1
                proj%non_proj_indexes(exc_pos)=i
             END IF
          END IF
       END DO
       proj%initialized=.TRUE.
    END IF
  END SUBROUTINE cp_proj_init
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_init2 [1.0] *
!!
!!   NAME
!!     cp_proj_init2
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_init2(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Intent (OUT):: proj
!!       Type(cp_sll_int_type), Pointer:: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_init2
!!
!!   FUNCTION
!!     initializes the given projection using a linked list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to be initialized
!!     - proj_indexes: the indexes of the full basis that are kept
!!     - size_full_bas: the size of the full basis
!!     - non_proj_indexes: the index of the excluded basis
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_proj_init2(proj, proj_indexes, size_full_bas, &
       error)
    TYPE(pao_bas_projection_type), INTENT(out) :: proj
    TYPE(cp_sll_int_type), POINTER ::proj_indexes
    INTEGER, INTENT(in) :: size_full_bas
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure,found
    INTEGER :: stat, size_min_bas, i, j, exc_pos, el_att
    CHARACTER(len=*), PARAMETER :: routineN='cp_proj_init2',&
         routineP=moduleN//':'//routineN
    TYPE(cp_sll_int_type), POINTER :: iterator
    failure=.FALSE.

    iterator => proj_indexes
    size_min_bas=0
    DO
       IF (.not.cp_next(iterator, el_att=el_att, error=error)) EXIT
       size_min_bas=size_min_bas+1
       CPPrecondition(el_att>0,cp_failure_level,routineP,error,failure)
    END DO
    CPPrecondition(size_min_bas<=size_full_bas,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(size_min_bas>0,cp_warning_level,routineP,error)
    ALLOCATE(proj%proj_indexes(size_min_bas),&
         proj%non_proj_indexes(size_full_bas-size_min_bas), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       iterator => proj_indexes
       DO i=1,size_min_bas
          IF (.not.cp_next(iterator, el_att=el_att, error=error)) EXIT
          proj%proj_indexes(i)=el_att
       END DO
       exc_pos=0
       DO i=1,size_full_bas
          found=.FALSE.
          DO j=1,SIZE(proj%proj_indexes)
             IF (i==proj%proj_indexes(j)) THEN
                found=.TRUE.
                EXIT
             END IF
          END DO
          IF (.not.found) THEN
             CPInvariant(exc_pos<SIZE(proj%non_proj_indexes),cp_failure_level,routineP,error,failure)
             IF (.not.failure) THEN
                exc_pos=exc_pos+1
                proj%non_proj_indexes(exc_pos)=i
             END IF
          END IF
       END DO
       proj%initialized=.TRUE.
    END IF
  END SUBROUTINE cp_proj_init2
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_proj_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_dealloc_ref(proj, error)
!!       Type(pao_bas_projection_type), Intent (INOUT):: proj
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory allocated by this object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to be deallocated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_proj_dealloc_ref(proj, error)
    TYPE(pao_bas_projection_type), INTENT(inout) :: proj
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: handle, stat
    CHARACTER(len=*), PARAMETER :: routineN='cp_proj_dealloc_ref',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    !call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DEALLOCATE(proj%proj_indexes, proj%non_proj_indexes, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       proj%initialized=.FALSE.
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_proj_dealloc_ref
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_valid [1.0] *
!!
!!   NAME
!!     cp_proj_valid
!!
!!   SYNOPSIS
!!     Function cp_proj_valid(proj, error) Result(res)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Logical:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_proj_valid
!!
!!   FUNCTION
!!     returns true if the projection is valid (minimal tests)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to test
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_proj_valid(proj,error) RESULT(res)
    TYPE(pao_bas_projection_type), INTENT(in) :: proj
    LOGICAL ::res
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_proj_valid',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(proj%proj_indexes),cp_warning_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(proj%non_proj_indexes),cp_warning_level,routineP,error,failure)
    res=.not.failure
  END FUNCTION cp_proj_valid
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_get [1.0] *
!!
!!   NAME
!!     cp_proj_get
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_get(proj, size_min_bas, size_full_bas,&
!!         size_excl_bas, proj_indexes, non_proj_indexes, error)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Integer, Optional, Intent (OUT):: size_excl_bas, size_full_bas,&
!!         size_min_bas
!!       Integer, Optional, Dimension(:), Pointer:: proj_indexes,&
!!         non_proj_indexes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_get
!!
!!   FUNCTION
!!     returns various attributes about the given projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection object
!!     - size_min_bas: the size of the minimal basis
!!     - size_full_bas: the size of the full basis
!!     - size_excl_bas: the size of the excluded basis
!!     - proj_indexes: the local indexes of the full basis that are kept
!!       by the projection
!!     - non_proj_indexes: the local indexes of the full basis that are
!!       not in the min_bas (i.e. the indexes of the excluded basis)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_proj_get(proj, size_min_bas, size_full_bas, &
       size_excl_bas, proj_indexes, non_proj_indexes, error)
    TYPE(pao_bas_projection_type), INTENT(in) :: proj
    INTEGER, OPTIONAL, INTENT(out) :: size_excl_bas,size_full_bas,&
         size_min_bas
    INTEGER, OPTIONAL, DIMENSION(:), POINTER :: proj_indexes,&
         non_proj_indexes
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(size_min_bas)) size_min_bas=SIZE(proj%proj_indexes)
       IF (PRESENT(size_full_bas)) &
            size_full_bas=SIZE(proj%proj_indexes)+SIZE(proj%non_proj_indexes)
       IF (PRESENT(size_excl_bas)) size_excl_bas=SIZE(proj%non_proj_indexes)
       IF (PRESENT(proj_indexes)) proj_indexes => proj%proj_indexes
       IF (PRESENT(non_proj_indexes)) non_proj_indexes => proj%non_proj_indexes
    END IF
  END SUBROUTINE cp_proj_get
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_get_proj_indexes [1.0] *
!!
!!   NAME
!!     cp_proj_get_proj_indexes
!!
!!   SYNOPSIS
!!     Function cp_proj_get_proj_indexes(proj, error) Result(res)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Integer, Dimension(:), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_proj_get_proj_indexes
!!
!!   FUNCTION
!!     returns the projection indexes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to get the info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_proj_get_proj_indexes(proj, error) RESULT(res)
    TYPE(pao_bas_projection_type), INTENT(in) :: proj
    INTEGER, DIMENSION(:), POINTER :: res
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_proj_get_proj_indexes',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    IF (cp_debug.and.debug_this_module) THEN
       CPPrecondition(cp_proj_valid(proj,error=error),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       res => proj%proj_indexes
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION cp_proj_get_proj_indexes
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_get_non_proj_indexes [1.0] *
!!
!!   NAME
!!     cp_proj_get_non_proj_indexes
!!
!!   SYNOPSIS
!!     Function cp_proj_get_non_proj_indexes(proj, error) Result(res)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Integer, Dimension(:), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_proj_get_non_proj_indexes
!!
!!   FUNCTION
!!     returns the indexes excluded by th projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to get the info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_proj_get_non_proj_indexes(proj, error) RESULT(res)
    TYPE(pao_bas_projection_type), INTENT(in) :: proj
    INTEGER, DIMENSION(:), POINTER :: res
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_proj_get_non_proj_indexes',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    IF (cp_debug.and.debug_this_module) THEN
       CPPrecondition(cp_proj_valid(proj,error=error),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       res => proj%non_proj_indexes
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION cp_proj_get_non_proj_indexes
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_write [1.0] *
!!
!!   NAME
!!     cp_proj_write
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_write(proj, unit_nr, error)
!!       Type(pao_bas_projection_type), Intent (IN), Target:: proj
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_write
!!
!!   FUNCTION
!!     writes a description of the projection to the given unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to output
!!     - unit_nr: the unit number to write to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_proj_write(proj,unit_nr,error)
  TYPE(pao_bas_projection_type), INTENT(in), TARGET :: proj
  INTEGER, INTENT(in) :: unit_nr
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_proj_write',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     WRITE(unit=unit_nr,fmt="(t2,a)")'# projection indexes'
     CALL cp_1d_i_write(proj%proj_indexes,unit_nr=unit_nr,error=error)
     WRITE(unit=unit_nr,fmt="(t2,a)")'# non projection indexes'
     CALL cp_1d_i_write(proj%non_proj_indexes,unit_nr=unit_nr,error=error)
  END IF
END SUBROUTINE cp_proj_write
!***************************************************************************

END MODULE pao_bas_projection_types

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_para_env [1.0] *
!!
!!   NAME
!!     cp_para_env
!!
!!   FUNCTION
!!     type to store parallelization informations (at the moment assumes 1d
!!     position and uses mpi)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_para_env
  USE cp_para_types,                   ONLY: cp_para_env_type, cp_para_cart_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_environ
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_para_env'

  PUBLIC :: cp_para_env_retain, cp_para_env_release, cp_para_env_create, &
      cp_para_env_duplicate,&
      cp_para_env_write
!***
!****************************************************************************
CONTAINS

!!****f* cp_para_env/cp_para_env_create [1.0] *
!!
!!   NAME
!!     cp_para_env_create
!!
!!   FUNCTION
!!     creates a new para environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: the new parallel environment
!!     - group: the id of the actual mpi_group
!!     - source: the id of the special (master) processor (defaults to 0)
!!     - mepos: the id of the actual processor
!!     - num_pe: the number of processors in the group
!!     - owns_group: if the group is owned by this object (defaults to true)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_para_env_create(para_env, group, source,mepos, num_pe,&
      owns_group,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: group
    INTEGER, INTENT(in), OPTIONAL            :: source, mepos, num_pe
    LOGICAL, INTENT(in), OPTIONAL            :: owns_group
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(para_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
      para_env%group=group
      para_env%source=0
      para_env%ref_count=1
      para_env%owns_group=.TRUE.
      IF (PRESENT(source)) para_env%source=source
      IF (PRESENT(owns_group)) para_env%owns_group=owns_group
      IF (.NOT.(PRESENT(mepos).AND.PRESENT(num_pe))) THEN
         CALL cp_para_env_update(para_env,error=error)
      ELSE
         para_env%mepos=mepos
         para_env%num_pe=num_pe
      END IF
    END IF
  END SUBROUTINE cp_para_env_create
!***************************************************************************

!!****f* cp_para_env/cp_para_env_duplicate [1.0] *
!!
!!   NAME
!!     cp_para_env_duplicate
!!
!!   FUNCTION
!!     creates a new para environment by duplicating an environment
!!     (and also the undelying mpi group)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: the new group
!!     - template: the environment to copy
!!     - source: the source of the new para_env, defaults to template%source
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_para_env_duplicate(para_env, TEMPLATE, source,error)
    TYPE(cp_para_env_type), POINTER          :: para_env, TEMPLATE
    INTEGER, INTENT(in), OPTIONAL            :: source
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_duplicate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(TEMPLATE),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
      CPPrecondition(template%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
      ALLOCATE(para_env,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
      para_env%ref_count=1
      CALL mp_comm_dup(template%group,para_env%group)
      para_env%source=template%source
      para_env%owns_group=.TRUE.
      IF (PRESENT(source)) para_env%source=source
      CALL cp_para_env_update(para_env,error=error)
    END IF
  END SUBROUTINE cp_para_env_duplicate
!***************************************************************************

!!****f* cp_para_env/cp_para_env_retain [1.0] *
!!
!!   NAME
!!     cp_para_env_retain
!!
!!   FUNCTION
!!     retains the para object (to be called when you want to keep a
!!     shared copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: the new group
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_para_env_retain(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
      para_env%ref_count=para_env%ref_count+1
    END IF
  END SUBROUTINE cp_para_env_retain
!***************************************************************************

!!****f* cp_para_env/cp_para_env_release [1.0] *
!!
!!   NAME
!!     cp_para_env_release
!!
!!   FUNCTION
!!     releases the para object (to be called when you don't want anymore
!!     the shared copy of this object)
!!
!!   NOTES
!!     to avoid circular dependencies cp_log_handling has a private copy
!!     of this method (see cp_log_handling:my_cp_para_env_release)!
!!
!!   INPUTS
!!     - para_env: the new group
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_para_env_release(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(para_env)) THEN
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
      para_env%ref_count=para_env%ref_count-1
      IF (para_env%ref_count<1) THEN
         IF (para_env%owns_group) THEN
            CALL mp_comm_free(para_env%group)
         END IF
         DEALLOCATE(para_env,stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF
    NULLIFY(para_env)
  END SUBROUTINE cp_para_env_release
!***************************************************************************

!!****f* cp_para_env/cp_para_env_update [1.0] *
!!
!!   NAME
!!     cp_para_env_update
!!
!!   FUNCTION
!!     gets again the position and size of the group from the mpi_group
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: the new group
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_para_env_update(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_update', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
      CALL mp_environ(taskid=para_env%mepos,numtask=para_env%num_pe,&
           groupid=para_env%group)
    END IF
  END SUBROUTINE cp_para_env_update
!***************************************************************************
!!****f* cp_para_env/cp_para_env_write [1.0] *
!!
!!   NAME
!!     cp_para_env_write
!!
!!   FUNCTION
!!     writes a description of the parallel environment to the given unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: the parallel environment to output
!!     - unit_nr: the unit where to output
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_para_env_write(para_env,unit_nr,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iostat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(para_env)) THEN
     WRITE (unit=unit_nr,&
          fmt="(' <cp_para_env>:{ owns_group=',l1,',')",&
          iostat=iostat) para_env%owns_group
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) para_env%group, para_env%ref_count
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  mepos=',i8,',')",&
          iostat=iostat) para_env%mepos
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  source=',i8,',')",&
          iostat=iostat) para_env%source
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  num_pe=',i8,'}')",&
          iostat=iostat) para_env%num_pe
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  ELSE
     WRITE (unit=unit_nr,&
          fmt="(a)", iostat=iostat) ' <cp_para_env>:*null* '
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE cp_para_env_write
!***************************************************************************

!!****f* cp_para_env/cp_cart_create *
!!
!!   NAME
!!     cp_cart_create
!!
!!   FUNCTION
!!     creates a cart (multidimensional parallel environment)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cart: the cart environment to create
!!     - group: the mpi communicator
!!     - ndims: the number of dimensions of the cart
!!     - owns_group: if this object owns the underlying cart (and should
!!       free it)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_cart_create(cart, group, ndims,owns_group, error)
    TYPE(cp_para_cart_type), POINTER         :: cart
    INTEGER, INTENT(in)                      :: group
    INTEGER, INTENT(in) :: ndims
    LOGICAL, INTENT(in), OPTIONAL :: owns_group
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cart_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(cart),cp_failure_level,routineP,error,failure)
    ALLOCATE(cart,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       cart%owns_group=.true.
       IF (PRESENT(owns_group)) cart%owns_group=owns_group
       cart%ndims=ndims
       cart%group=group
       
       ALLOCATE(cart%source(ndims),cart%periodic(ndims),cart%mepos(ndims),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       
       cart%source=0
       cart%mepos=0
       cart%periodic=.false.
       cart%ref_count=1
       cart%ntask=1
       CALL cp_cart_update(cart,error=error)
    END IF
  END SUBROUTINE cp_cart_create
!***************************************************************************

!!****f* cp_para_env/cp_cart_update *
!!
!!   NAME
!!     cp_cart_update
!!
!!   FUNCTION
!!     updates the information about the given cart
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cart: the cart to update
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_cart_update(cart,error)
  TYPE(cp_para_cart_type), pointer :: cart
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_cart_update',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(ASSOCIATED(cart),cp_failure_level,routineP,error,failure)
  CPPrecondition(cart%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL mp_environ( cart%group, cart%ndims, cart%num_pe, task_coor=cart%mepos, &
          periods=cart%periodic)
     CALL mp_environ( numtask=cart%ntask, taskid=cart%rank, groupid=cart%group)
  END IF
END SUBROUTINE cp_cart_update
!***************************************************************************

!!****f* cp_para_env/cp_cart_release *
!!
!!   NAME
!!     cp_cart_release
!!
!!   FUNCTION
!!     releases the given cart
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cart: the cart to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_cart_release(cart,error)
  TYPE(cp_para_cart_type), pointer :: cart
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_cart_release',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  
  IF (ASSOCIATED(cart)) THEN
     CPPrecondition(cart%ref_count>0,cp_failure_level,routineP,error,failure)
     cart%ref_count=cart%ref_count-1
     IF (cart%ref_count==0) THEN
        if (cart%owns_group) then
            CALL mp_comm_free(cart%group)
         END IF
        deallocate(cart%source,cart%periodic,cart%mepos,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(cart,stat=stat)
     END IF
  END IF
  nullify(cart)
END SUBROUTINE cp_cart_release
!***************************************************************************

!!****f* cp_para_env/cp_cart_retain *
!!
!!   NAME
!!     cp_cart_retain
!!
!!   FUNCTION
!!     retains the given cart
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cart: the cart to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_cart_retain(cart,error)
  TYPE(cp_para_cart_type), pointer :: cart
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_cart_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(ASSOCIATED(cart),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(cart%ref_count>0,cp_failure_level,routineP,error,failure)
     cart%ref_count=cart%ref_count+1
  END IF
END SUBROUTINE cp_cart_retain
!***************************************************************************

END MODULE cp_para_env

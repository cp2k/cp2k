!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief type to store parallelization informations (at the moment assumes 1d
!>      position and uses mpi)
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_para_env
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE f77_blas
  USE message_passing,                 ONLY: mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_environ
#include "cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_para_env'

  PUBLIC :: cp_para_env_retain, cp_para_env_release, cp_para_env_create,&
            cp_para_env_write
  PUBLIC :: cp_cart_create, cp_cart_release,&
            cp_cart_write
!***
CONTAINS

! *****************************************************************************
!> \brief creates a new para environment
!> \param para_env the new parallel environment
!> \param group the id of the actual mpi_group
!> \param source the id of the special (master) processor (defaults to 0)
!> \param mepos the id of the actual processor
!> \param num_pe the number of processors in the group
!> \param owns_group if the group is owned by this object (defaults to true)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_para_env_create(para_env, group, source,mepos, num_pe,&
      owns_group,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: group
    INTEGER, INTENT(in), OPTIONAL            :: source, mepos, num_pe
    LOGICAL, INTENT(in), OPTIONAL            :: owns_group
    TYPE(cp_error_type), INTENT(inout)       :: error

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
      para_env%ionode=para_env%mepos==para_env%source
    END IF
  END SUBROUTINE cp_para_env_create

! *****************************************************************************
!> \brief creates a new para environment by duplicating an environment
!>      (and also the undelying mpi group)
!> \param para_env the new group
!> \param TEMPLATE ...
!> \param source the source of the new para_env, defaults to template%source
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param template the environment to copy UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_para_env_duplicate(para_env, TEMPLATE, source,error)
    TYPE(cp_para_env_type), POINTER          :: para_env, TEMPLATE
    INTEGER, INTENT(in), OPTIONAL            :: source
    TYPE(cp_error_type), INTENT(inout)       :: error

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
      para_env%ionode=para_env%mepos==para_env%source
    END IF
  END SUBROUTINE cp_para_env_duplicate

! *****************************************************************************
!> \brief retains the para object (to be called when you want to keep a
!>      shared copy of this object)
!> \param para_env the new group
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_para_env_retain(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief releases the para object (to be called when you don't want anymore
!>      the shared copy of this object)
!> \param para_env the new group
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      to avoid circular dependencies cp_log_handling has a private copy
!>      of this method (see cp_log_handling:my_cp_para_env_release)!
! *****************************************************************************
  SUBROUTINE cp_para_env_release(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief gets again the position and size of the group from the mpi_group
!> \param para_env the new group
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_para_env_update(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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
      para_env%ionode=para_env%mepos==para_env%source
    END IF
  END SUBROUTINE cp_para_env_update

! *****************************************************************************
!> \brief writes a description of the parallel environment to the given unit
!> \param para_env the parallel environment to output
!> \param unit_nr the unit where to output
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_para_env_write(para_env,unit_nr,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief creates a cart (multidimensional parallel environment)
!> \param cart the cart environment to create
!> \param group the mpi communicator
!> \param ndims the number of dimensions of the cart
!> \param owns_group if this object owns the underlying cart (and should
!>        free it)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_cart_create(cart, group, ndims,owns_group, error)
    TYPE(cp_para_cart_type), POINTER         :: cart
    INTEGER, INTENT(in)                      :: group, ndims
    LOGICAL, INTENT(in), OPTIONAL            :: owns_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cart_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(cart),cp_failure_level,routineP,error,failure)
    ALLOCATE(cart,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       cart%owns_group=.TRUE.
       IF (PRESENT(owns_group)) cart%owns_group=owns_group
       cart%ndims=ndims
       cart%group=group

       ALLOCATE(cart%source(ndims),cart%periodic(ndims),cart%mepos(ndims),&
            cart%num_pe(ndims),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

       cart%source=0
       cart%mepos=0
       cart%periodic=.FALSE.
       cart%ref_count=1
       cart%ntask=1
       CALL cp_cart_update(cart,error=error)
    END IF
  END SUBROUTINE cp_cart_create

! *****************************************************************************
!> \brief updates the information about the given cart
!> \param cart the cart to update
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_cart_update(cart,error)
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cart_update', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(cart),cp_failure_level,routineP,error,failure)
  CPPrecondition(cart%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL mp_environ( cart%group, cart%ndims, cart%num_pe, task_coor=cart%mepos, &
          periods=cart%periodic)
     CALL mp_environ( numtask=cart%ntask, taskid=cart%rank, groupid=cart%group)
  END IF
END SUBROUTINE cp_cart_update

! *****************************************************************************
!> \brief releases the given cart
!> \param cart the cart to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_cart_release(cart,error)
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cart_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(cart)) THEN
     CPPrecondition(cart%ref_count>0,cp_failure_level,routineP,error,failure)
     cart%ref_count=cart%ref_count-1
     IF (cart%ref_count==0) THEN
        IF (cart%owns_group) THEN
            CALL mp_comm_free(cart%group)
         END IF
        DEALLOCATE(cart%source,cart%periodic,cart%mepos,cart%num_pe,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(cart,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(cart)
END SUBROUTINE cp_cart_release

! *****************************************************************************
!> \brief retains the given cart
!> \param cart the cart to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_cart_retain(cart,error)
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cart_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(cart),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(cart%ref_count>0,cp_failure_level,routineP,error,failure)
     cart%ref_count=cart%ref_count+1
  END IF
END SUBROUTINE cp_cart_retain

! *****************************************************************************
!> \brief writes out the informations about the cart
!> \param cart the parallel environment to output
!> \param unit_nr the unit where to output
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_cart_write(cart,unit_nr,error)
    TYPE(cp_para_cart_type), POINTER         :: cart
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cart_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iostat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(cart)) THEN
     WRITE (unit=unit_nr,&
          fmt="(' <cp_cart>:{ owns_group=',l1,',')",&
          iostat=iostat) cart%owns_group
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) cart%group, cart%ref_count
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  ndims=',i8,',')",&
          iostat=iostat) cart%ndims
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  rank=',i8,', ntask=',i8,',')",&
          iostat=iostat) cart%rank, cart%ntask
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  mepos=',10i8)",&
          iostat=iostat) cart%mepos
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  source=',10i8,',')",&
          iostat=iostat) cart%source
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  num_pe=',10i8,'}')",&
          iostat=iostat) cart%num_pe
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  periodic=',10l2,'}')",&
          iostat=iostat) cart%periodic
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  ELSE
     WRITE (unit=unit_nr,&
          fmt="(a)", iostat=iostat) ' <cp_cart>:*null* '
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE cp_cart_write

END MODULE cp_para_env

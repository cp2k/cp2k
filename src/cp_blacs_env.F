!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_blacs_env [1.0] *
!!
!!   NAME
!!     cp_blacs_env
!!
!!   FUNCTION
!!     methods related to the blacs parallel environment
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
module cp_blacs_env
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use kinds, only: int_size,&
                   wp => dp,&
                   wp_size => dp_size
  use cp_para_types, only: cp_blacs_env_type
  USE message_passing,     ONLY: mp_bcast,&
                                 mp_max,&
                                 mp_sum,&
                                 mp_sync
  USE termination,         ONLY: stop_memory,&
                                 stop_program
  use timings, only: timeset, timestop
  USE global_types,        ONLY: global_environment_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_blacs_env'

  public :: cp_blacs_env_type ! make it accessible only through cp_para_types?
  public :: cp_blacs_env_create, cp_blacs_env_retain, cp_blacs_env_release
  public :: start_blacs, finish_blacs, cp_blacs_env_write,&
       cp_blacs_env_from_globenv
!***
!****************************************************************************
contains

!!****f* qs_blacs/cp_blacs_env_create [1.0] *
!!
!!   NAME
!!     cp_blacs_env_create
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_create(blacs_env, group, owns_group, source,&
!!         error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Integer, Intent (IN):: group
!!       Logical, Intent (IN), Optional:: owns_group
!!       Integer, Dimension(2), Intent (IN), Optional:: source
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_create
!!
!!   FUNCTION
!!     allocates and initializes a type that represent a blacs context
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the type to initialize
!!     - group: the id of the blacs context
!!     - owns_group: if this type owns the group and should release it
!!       (via blacs_gridexit) when deallocated
!!     - source: position of a special processor (defaults to (0,0))
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
subroutine cp_blacs_env_create(blacs_env,group,owns_group,source,error)
  type(cp_blacs_env_type), pointer :: blacs_env
  integer, intent(in) :: group
  logical, intent(in), optional :: owns_group
  integer, dimension(2), intent(in), optional :: source
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_blacs_env_create',&
        routineP=moduleN//':'//routineN
  integer :: stat
  failure=.false.

#ifdef __parallel
#ifndef __BLACS
  call cp_assert(.false.,cp_failure_level, cp_assertion_failed,&
       fromWhere=routineP,message="to use the blacs environment "//&
       "(used by full matrices) you need the blacs library "//&
       CPSourceFileRef,&
       error=error)
#endif
#endif

  allocate(blacs_env,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     blacs_env%group=group
     blacs_env%owns_group=.true.
     blacs_env%ref_count=1
     blacs_env%mepos=0
     blacs_env%num_pe=1
     blacs_env%source=0
     blacs_env%my_pid=0
     blacs_env%n_pid=1
     if (present(owns_group)) blacs_env%owns_group=owns_group
     if (present(source)) then
        blacs_env%source=source
     end if
#ifdef __BLACS
     call blacs_pinfo(blacs_env%my_pid,blacs_env%n_pid)
     call blacs_gridinfo(blacs_env%group,blacs_env%num_pe(1),&
          blacs_env%num_pe(2),blacs_env%mepos(1),blacs_env%mepos(2))
#endif
  end if
end subroutine cp_blacs_env_create
!***************************************************************************

!!****f* qs_blacs/cp_blacs_env_retain [1.0] *
!!
!!   NAME
!!     cp_blacs_env_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_retain(blacs_env, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_retain
!!
!!   FUNCTION
!!     retains the given blacs env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the blacs env to retain
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
subroutine cp_blacs_env_retain(blacs_env,error)
  type(cp_blacs_env_type), pointer :: blacs_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_blacs_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.false.
  CPPrecondition(associated(blacs_env),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPrecondition(blacs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     blacs_env%ref_count=blacs_env%ref_count+1
  end if
end subroutine cp_blacs_env_retain
!***************************************************************************

!!****f* qs_blacs/cp_blacs_env_release [1.0] *
!!
!!   NAME
!!     cp_blacs_env_release
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_release(blacs_env, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_release
!!
!!   FUNCTION
!!     releases the given blacs_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the blacs env to relase
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
subroutine cp_blacs_env_release(blacs_env,error)
  type(cp_blacs_env_type), pointer :: blacs_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_blacs_env_release',&
        routineP=moduleN//':'//routineN
  integer :: stat
  failure=.false.
  
  if (associated(blacs_env)) then
     CPPrecondition(blacs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     blacs_env%ref_count=blacs_env%ref_count-1
     if (blacs_env%ref_count<1) then
        if (blacs_env%owns_group) then
#ifdef __BLACS
           call blacs_gridexit(blacs_env%group)
#endif
        end if
        deallocate(blacs_env,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     end if
  end if
  nullify(blacs_env)
end subroutine cp_blacs_env_release
!***************************************************************************

!!****f* cp_blacs_env/cp_blacs_env_write [1.0] *
!!
!!   NAME
!!     cp_blacs_env_write
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_write(blacs_env, unit_nr, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_write
!!
!!   FUNCTION
!!     writes the description of the given blacs env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the blacs environment to write
!!     - unit_nr: the unit number where to write the description of the
!!       blacs environment
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
subroutine cp_blacs_env_write(blacs_env, unit_nr, error)
  type(cp_blacs_env_type), pointer :: blacs_env
  integer, intent(in) :: unit_nr
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_blacs_env_write',&
        routineP=moduleN//':'//routineN
  integer :: iostat
  failure=.false.
  
  if (associated(blacs_env)) then
     write (unit=unit_nr,&
          fmt="(' <cp_blacs_env>:{ owns_group=',l1,',')",&
          iostat=iostat) blacs_env%owns_group
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) blacs_env%group, blacs_env%ref_count
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  mepos=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%mepos(1), blacs_env%mepos(2)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  source=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%source(1), blacs_env%source(2)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  num_pe=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%num_pe(1), blacs_env%num_pe(2)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  my_pid=',i10,', n_pid=',i10,' }')",&
          iostat=iostat) blacs_env%my_pid, blacs_env%n_pid
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  else
     write (unit=unit_nr,&
          fmt="(a)", iostat=iostat) ' <cp_blacs_env>:*null* '
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  end if
end subroutine cp_blacs_env_write
!***************************************************************************

!!****f* cp_blacs_env/cp_blacs_env_from_globenv [1.0] *
!!
!!   NAME
!!     cp_blacs_env_from_globenv
!!
!!   SYNOPSIS
!!     Function cp_blacs_env_from_globenv(global_env, error) Result(res)
!!       Type(cp_blacs_env_type), Pointer:: res
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_blacs_env_from_globenv
!!
!!   FUNCTION
!!     returns the blacs environment stored in the given global_env.
!!     if null issues a warning and initializes it
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - global_env: the global environment
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
function cp_blacs_env_from_globenv(global_env,error) result(res)
  type(cp_blacs_env_type), pointer :: res
  type(global_environment_type), intent(in),target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_blacs_env_from_globenv',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(global_env%blacs_env),cp_failure_level,routineP,error,failure)
  !if (failure) then
  !   call start_blacs(global_env) ! needs intent inout
  !end if
  res => global_env%blacs_env
end function cp_blacs_env_from_globenv
!***************************************************************************

  SUBROUTINE start_blacs(globenv)

!   Purpose: Initialize a BLACS process grid. The BLACS context is returned.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE start_blacs (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: context,group,ipe,istat,mypcol,mype,myprow,npcol,npe,nprow,&
               output_unit
    LOGICAL :: ionode
    type(cp_error_type) :: error

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol,pcon,prow

!   ---------------------------------------------------------------------------

#ifdef __parallel
#ifndef __BLACS
  call cp_assert(.false.,cp_failure_level, cp_assertion_failed,&
       fromWhere=routine,message="to use the blacs environment "//&
       "(used by full matrices) you need the blacs library "//&
       CPSourceFileRef,&
       error=error)
#endif
#endif

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    context = 0
    npcol = 0
    nprow = 0
#if defined(__BLACS)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_get(-1,0,context)

    IF (nprow*npcol /= npe) THEN
      DO ipe=CEILING(SQRT(REAL(npe,wp))),npe
        IF (MODULO(npe,ipe) == 0) THEN
          nprow = ipe
          npcol = npe/nprow
          EXIT
        END IF
      END DO
    END IF

    CALL blacs_gridinit(context,"Row-major",nprow,npcol)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    IF (globenv%print%blacs_info) THEN

      ALLOCATE (prow(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"prow",npe*int_size)
      prow(:) = 0
      ALLOCATE (pcol(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcol",npe*int_size)
      pcol(:) = 0
      ALLOCATE (pcon(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon",npe*int_size)
      pcon(:) = 0

      prow(mype) = myprow
      pcol(mype) = mypcol
      pcon(mype) = context

      CALL mp_sum(prow,group)
      CALL mp_sum(pcol,group)
      CALL mp_sum(pcon,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,(T3,A,T32,I6))")&
          "BLACS INFORMATION (BLACS started)",&
          "Number of processes:         ",nprow*npcol,&
          "Number of process rows:      ",nprow,&
          "Number of process columns:   ",npcol
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
          " PE   process row   process column   BLACS context"
        WRITE (UNIT=output_unit,FMT="(I5,T14,I6,T31,I6,T41,I12)")&
          (ipe,prow(ipe),pcol(ipe),pcon(ipe),ipe=0,npe-1)
      END IF

      DEALLOCATE (prow,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"prow")
      DEALLOCATE (pcol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcol")
      DEALLOCATE (pcon,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon")

    END IF

#endif
    call cp_blacs_env_create(globenv%blacs_env,group = context)

  END SUBROUTINE start_blacs

! *****************************************************************************

  SUBROUTINE finish_blacs(globenv)

!   Purpose: Release the resources of a BLACS context.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE finish_blacs (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: group,ipe,istat,mype,npe,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcon
    type(cp_blacs_env_type), pointer :: context

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    context => globenv%blacs_env
#if defined(__BLACS)

    IF (globenv%print%blacs_info) THEN

      mype=context%n_pid
      npe=context%my_pid

      ALLOCATE (pcon(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon",npe*int_size)
      pcon(:) = 0
      
      pcon(mype) = context%group

      CALL mp_sum(pcon,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "BLACS INFORMATION (BLACS finished)"
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
          " PE   BLACS context"
        WRITE (UNIT=output_unit,FMT="(I5,T10,I12)")&
          (ipe,pcon(ipe),ipe=0,npe-1)
      END IF

      DEALLOCATE (pcon,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon")

    END IF

    CALL mp_sync(group)

    call cp_blacs_env_release(globenv%blacs_env)
    !    CALL blacs_gridexit(context)

    !globenv%context = 0
    !globenv%nprow = 0
    !globenv%npcol = 0

#endif
  END SUBROUTINE finish_blacs

! *****************************************************************************

end module cp_blacs_env

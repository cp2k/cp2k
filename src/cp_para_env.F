!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
module cp_para_env
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
      cp_note_level
  use cp_error_handling, only: cp_debug, cp_internal_error, cp_assertion_failed,&
      cp_assert, cp_unimplemented_error, cp_error_type
  use kinds, only: wp=>dp
  use timings, only: timeset, timestop
  use cp_para_types, only: cp_para_env_type
  use message_passing, only: mp_comm_dup, mp_comm_free, mp_environ
  use global_types, only: global_environment_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_para_env'

  public :: cp_para_env_type
  public :: cp_para_retain, cp_para_release, cp_para_create, &
      cp_para_duplicate, cp_para_env_from_globenv,&
      cp_para_env_write
!***
!****************************************************************************
contains

!!****f* cp_para_env/cp_para_create [1.0] *
!!
!!   NAME
!!     cp_para_create
!!
!!   SYNOPSIS
!!     Subroutine cp_para_create(para_env, group, source, mepos, num_pe,&
!!         owns_group, error)
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Intent (IN):: group
!!       Integer, Intent (IN), Optional:: source, mepos, num_pe
!!       Logical, Intent (IN), Optional:: owns_group
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_para_create
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
  subroutine cp_para_create(para_env, group, source,mepos, num_pe,&
      owns_group,error)
    type(cp_para_env_type), pointer :: para_env
    integer, intent(in) :: group
    integer, intent(in), optional :: source,mepos, num_pe
    logical, intent(in), optional :: owns_group
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_para_create',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    allocate(para_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
      para_env%group=group
      para_env%source=0
      para_env%ref_count=1
      para_env%owns_group=.true.
      if (present(source)) para_env%source=source
      if (present(owns_group)) para_env%owns_group=owns_group
      if (.not.(present(mepos).and.present(num_pe))) then
         call cp_para_update(para_env,error=error)
      else
         para_env%mepos=mepos
         para_env%num_pe=num_pe
      end if
    end if
  end subroutine cp_para_create
!***************************************************************************

!!****f* cp_para_env/cp_para_duplicate [1.0] *
!!
!!   NAME
!!     cp_para_duplicate
!!
!!   SYNOPSIS
!!     Subroutine cp_para_duplicate(para_env, template, source, error)
!!       Type(cp_para_env_type), Pointer:: para_env, template
!!       Integer, Intent (IN), Optional:: source
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_para_duplicate
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
  subroutine cp_para_duplicate(para_env, template, source,error)
    type(cp_para_env_type), pointer :: para_env, template
    integer, intent(in), optional :: source
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_para_duplicate',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    CPPrecondition(associated(template),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
      CPPrecondition(template%ref_count>0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
      allocate(para_env,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
      para_env%ref_count=1
      call mp_comm_dup(template%group,para_env%group)
      para_env%source=template%source
      para_env%owns_group=.true.
      if (present(source)) para_env%source=source
      call cp_para_update(para_env,error=error)
    end if
  end subroutine cp_para_duplicate
!***************************************************************************

!!****f* cp_para_env/cp_para_retain [1.0] *
!!
!!   NAME
!!     cp_para_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_para_retain(para_env, error)
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_para_retain
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
  subroutine cp_para_retain(para_env, error)
    type(cp_para_env_type), pointer :: para_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_para_retain',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    CPPrecondition(associated(para_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
      para_env%ref_count=para_env%ref_count+1
    end if
  end subroutine cp_para_retain
!***************************************************************************

!!****f* cp_para_env/cp_para_release [1.0] *
!!
!!   NAME
!!     cp_para_release
!!
!!   SYNOPSIS
!!     Subroutine cp_para_release(para_env, error)
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_para_release
!!
!!   FUNCTION
!!     releases the para object (to be called when you don't want anymore
!!     the shared copy of this object)
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
  subroutine cp_para_release(para_env, error)
    type(cp_para_env_type), pointer :: para_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_para_release',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    if (associated(para_env)) then
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
      para_env%ref_count=para_env%ref_count-1
      if (para_env%ref_count<1) then
         if (para_env%owns_group) then
            call mp_comm_free(para_env%group)
         end if
         deallocate(para_env,stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      end if
    end if
    nullify(para_env)
  end subroutine cp_para_release
!***************************************************************************

!!****f* cp_para_env/cp_para_update [1.0] *
!!
!!   NAME
!!     cp_para_update
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
  subroutine cp_para_update(para_env, error)
    type(cp_para_env_type), pointer :: para_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_para_update',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    CPPrecondition(associated(para_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
      call mp_environ(taskid=para_env%mepos,numtask=para_env%num_pe,&
           groupid=para_env%group)
    end if
  end subroutine cp_para_update
!***************************************************************************

!!****f* cp_para_env/cp_para_env_from_globenv [1.0] *
!!
!!   NAME
!!     cp_para_env_from_globenv
!!
!!   FUNCTION
!!     returns the parallel environment stored in the global environment
!!     gives a warning if not initialized
!!
!!   NOTES
!!     should try to generate a para env??
!!
!!   INPUTS
!!     - global_env: the global environment that contains the para_env
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
function cp_para_env_from_globenv(global_env, error) result(res)
  type(cp_para_env_type),pointer :: res
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_para_env_from_globenv',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(global_env%para_env),cp_failure_level,routineP,error,failure)
  res => global_env%para_env
end function cp_para_env_from_globenv
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
subroutine cp_para_env_write(para_env,unit_nr,error)
  type(cp_para_env_type), pointer :: para_env
  integer, intent(in) :: unit_nr
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_blacs_env_write',&
        routineP=moduleN//':'//routineN
  integer :: iostat
  failure=.false.
  
  if (associated(para_env)) then
     write (unit=unit_nr,&
          fmt="(' <cp_para_env>:{ owns_group=',l1,',')",&
          iostat=iostat) para_env%owns_group
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) para_env%group, para_env%ref_count
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  mepos=',i8,',')",&
          iostat=iostat) para_env%mepos
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  source=',i8,',')",&
          iostat=iostat) para_env%source
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     write (unit=unit_nr,fmt="('  num_pe=',i8,'}')",&
          iostat=iostat) para_env%num_pe
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  else
     write (unit=unit_nr,&
          fmt="(a)", iostat=iostat) ' <cp_para_env>:*null* '
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  end if
end subroutine cp_para_env_write
!***************************************************************************

end module cp_para_env

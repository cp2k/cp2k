!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_mpi_env [1.0] *
!!
!!   NAME
!!     cp_mpi_env
!!
!!   FUNCTION
!!     type to store mpi informations
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
module cp_mpi_env
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_mpi_env'

!***
!****************************************************************************

!!****s* cp_mpi_env/cp_mpi_env_type [1.0] *
!!
!!   NAME
!!     cp_mpi_env_type
!!
!!   FUNCTION
!!     stores all the informations relevanto to an mpi environement
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     mepos: id of the actual processor
!!     group: id of the actual mpi group
!!     num_pe: number of processors in the group
!!     source: id of a special processor (for example the one for i-o)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  type cp_mpi_env_type
     integer :: mepos, group, source, num_pe
  end type cp_mpi_env_type
!!***
  !****************************************************************************

  public :: cp_mpi_env_type
  public :: cp_mpi_get_from_glob
contains

!!****f* cp_mpi_env/cp_mpi_get_from_glob [1.0] *
!!
!!   NAME
!!     cp_mpi_get_from_glob
!!
!!   SYNOPSIS
!!     Function cp_mpi_get_from_glob(global_env, error) Result(res)
!!       Type(cp_mpi_env_type):: res
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_mpi_get_from_glob
!!
!!   FUNCTION
!!     returns the mpi environment of the global_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - global_env: the global_env that contains the default mpi info
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  function cp_mpi_get_from_glob(global_env,error) result(res)
    type(cp_mpi_env_type) ::res
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_mpi_get_from_glob',&
         routineP=moduleN//':'//routineN
    failure=.false.

    res%mepos=global_env%mepos
    res%group=global_env%group
    res%source=global_env%source
    res%num_pe=global_env%num_pe
  end function cp_mpi_get_from_glob
!***************************************************************************

end module cp_mpi_env

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_qs_env_methods [1.0] *
!!
!!   NAME
!!     pao_qs_env_methods
!!
!!   FUNCTION
!!     pao related methods that work on the qs_env methods
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
module pao_qs_env_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types
  use pao_bas_projection_types
  use qs_environment_types
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind,&
       get_atomic_kind_set
  use particle_types
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_qs_env_methods'

  public :: cp_get_b_projection
  ! underlying routines
  public :: pao_qs_env_get_b_projection, pao_qs_env_get
  interface cp_get_b_projection
     module procedure pao_qs_env_get_b_projection
  end interface
!***
!****************************************************************************
contains

!!****f* pao_qs_env_methods/pao_qs_env_get_b_projection [1.0] *
!!
!!   NAME
!!     pao_qs_env_get_b_projection
!!
!!   SYNOPSIS
!!     Function pao_qs_env_get_p_index(pao_qs_env, atom_nr, i, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_env_type), Intent (IN):: pao_env
!!       Integer, Intent (IN):: atom_nr, i
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_env_get_p_index
!!
!!   FUNCTION
!!     returns the base projection projection on atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     qs_env: the qs environement you get the data from
!!     atom_nr: the number of the atom
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 1.2002 created [fawzi]
!!     - 19.04.02 used new projection
!!     - 23.04.02 moved to qs
!!
!!*** **********************************************************************
  function pao_qs_env_get_b_projection(qs_env, atom_nr, error) result(res)
    type(pao_bas_projection_type), pointer ::res
    type(qs_environment_type), intent(in) :: qs_env
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_qs_env_get_b_projection',&
        routineP=moduleN//':'//routineN
    type(particle_type), dimension(:), pointer :: particle_set
    failure=.false.

    call get_qs_env(qs_env, particle_set=particle_set)
    call get_atomic_kind(particle_set(atom_nr)%atomic_kind,pao_bas_proj=res)
  end function pao_qs_env_get_b_projection
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_get [1.0] *
!!
!!   NAME
!!     pao_qs_env_get
!!
!!   FUNCTION
!!     returns pao related info from the qs_environment_type
!!
!!   NOTES
!!     use this function rather than get_atomic_kind_set: some values
!!     might be cached in pao_env in the future
!!
!!   INPUTS
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
subroutine pao_qs_env_get(qs_env, global_env, pao_max_min_bas,&
     pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
     pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas, error)
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  integer, intent(out), optional :: pao_max_min_bas,&
     pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
     pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_qs_env_get',&
        routineP=moduleN//':'//routineN
  type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
  failure=.false.
  
  call get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
  CPPrecondition(associated(atomic_kind_set),cp_failure_level,routineP,error,failure)
  call get_atomic_kind_set(atomic_kind_set,pao_max_min_bas=pao_max_min_bas,&
       pao_max_full_bas=pao_max_full_bas,&
       pao_max_excl_bas=pao_max_excl_bas, pao_tot_min_bas=pao_tot_min_bas,&
       pao_tot_full_bas=pao_tot_full_bas, pao_tot_excl_bas=pao_tot_excl_bas,&
       pao_max_l_angles=pao_max_l_angles)
end subroutine pao_qs_env_get
!***************************************************************************

end module pao_qs_env_methods

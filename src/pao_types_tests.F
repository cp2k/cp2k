!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_types_tests [1.0] *
!!
!!   NAME
!!     pao_types_tests
!!
!!   FUNCTION
!!     tests the pao types
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module pao_types_tests
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use pao_utils, only:  cp_get_random_integer
  use cp_matrix_utils, only: cp_sup_norm_matrix
  use pao_types
  use pao_proj_methods
  use pao_env_methods
  use pao_matrix_methods
  use pao_l_angles_methods
  use pao_glob_angles_methods
  implicit none

  logical, private, parameter :: debug_this_module=.true.

!!***
!****************************************************************************
contains

!!****f* pao_types_tests/cp_test_pao_types [1.0] *
!!
!!   NAME
!!     cp_test_pao_types
!!
!!   SYNOPSIS
!!     Subroutine cp_test_pao_types(g_env, failed, error)
!!       Type(global_environment_type), Intent (IN):: g_env
!!       Logical, Intent (INOUT):: failed
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_test_pao_types
!!
!!   FUNCTION
!!     does some tests on the pao_types module
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     g_env: global environement (for the future)
!!     failure: an in out value that become true if some tests fail
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_test_pao_types(g_env, failed,error)
    type(global_environment_type), intent(in) ::g_env
    logical, intent(inout) :: failed
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_types_tests:cp_test_pao_types'
    type(pao_env_type) :: pao_env
    type(pao_glob_angles_type) :: glob_angles
    type(pao_local_angles_type) :: local_angles
    type(pao_matrix_type) :: pao_t
    integer, parameter :: atomNr=3
    integer, dimension(atomNr) :: minB, fullB
    integer :: i

    failure=.false.
    if (cp_debug.and.debug_this_module) then
       call cp_log(cp_error_get_logger(error),cp_note_level,routineP, &
            'began pao_types tests')
    end if
    do i=1,atomNr
       fullB(i)=cp_get_random_integer(from=1,to=30,error=error)
    end do
    do i=1,atomNr
       minB(i)=cp_get_random_integer(from=1,to=fullB(i),error=error)
    end do
    call cp_init(pao_env, nr_min_bas=minB, nr_full_bas=fullB, &
         global_env=g_env,error=error)
    call cp_init(glob_angles, pao_env=pao_env)
    call cp_init(local_angles, glob_angles=glob_angles, &
         atom_nr=cp_get_random_integer(from=1,to=atomNr))
    call cp_init(pao_t,pao_env=pao_env, min_bas_cols=.false.,&
         min_bas_rows=.false.)
    call cp_set(pao_env, angles_att=glob_angles)

    CPAssert(cp_validate(pao_env,error=error),cp_warning_level,routineP,error,failure)
    CPAssert(cp_validate(local_angles,error=error),cp_warning_level,routineP,error,failure)
    call cp_dealloc_ref(pao_env,error=error)
    call cp_dealloc_ref(local_angles,error=error)

    call cp_dealloc_ref(glob_angles,error=error)
    call cp_dealloc_ref(pao_t, error=error)
    if (failure) then
       call cp_error_message(cp_warning_level,routineP,&
            "cp_test_pao_types FAILED!",error)
    else if (cp_debug .and. debug_this_module) then
       call cp_error_message(cp_warning_level,routineP,&
            "cp_test_pao_types SUCCESS!",error) ! bug in print utils
       !print *, "cp_test_pao_types SUCCESS!"
    end if
    failed=failure
  end subroutine cp_test_pao_types
  !***************************************************************************

!!****f* pao_types_tests/cp_init_test_pao_env [1.0] *
!!
!!   NAME
!!     cp_init_test_pao_env
!!
!!   SYNOPSIS
!!     Subroutine cp_init_test_pao_env(pao_env, glob_env, atom_nr, error)
!!       Type(pao_env_type), Intent (OUT):: pao_env
!!       Type(global_environment_type), Target, Intent (IN):: glob_env
!!       Integer, Optional:: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_init_test_pao_env
!!
!!   FUNCTION
!!     initializes a valid pao env for testing purposes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the pao environement that gets initialized
!!     glob_env:glob_env pointer for pao_env
!!     atom_nr: the number of atoms
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_init_test_pao_env(pao_env, glob_env, atom_nr, error)
    type(pao_env_type), intent(out) :: pao_env
    type(global_environment_type), target, intent(in) :: glob_env
    integer, optional :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_types_tests:cp_init_test_pao_env'
    integer, dimension(:), allocatable :: minB, fullB
    integer :: n_atom,i, stat
    type(pao_matrix_type), pointer :: u
    nullify(u)
    failure=.false.

    n_atom=3
    if (present(atom_nr)) n_atom=atom_nr
    CPPrecondition(n_atom>0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       allocate(minB(n_atom),fullB(n_atom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          do i=1,n_atom
             fullB(i)=cp_get_random_integer(from=1,to=30,error=error)
          end do
          do i=1,n_atom
             minB(i)=cp_get_random_integer(from=1,to=fullB(i),error=error)
          end do
          call cp_init(pao_env, nr_min_bas=minB, nr_full_bas=fullB, &
               global_env=glob_env,error=error)
          call cp_set(pao_env,error=error) ! just to be sure
          CPPostcondition(cp_valid(pao_env,error=error),cp_warning_level,routineP,error,failure)
       end if
       deallocate(minB,fullB,stat=stat) !just to be sure...
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
  end subroutine cp_init_test_pao_env
  !***************************************************************************

end module pao_types_tests


!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_pao_types_tests [1.0] *
!!
!!   NAME
!!     qs_pao_types_tests
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
module qs_pao_types_tests
  use qs_log_handling
  use qs_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use qs_pao_utils, only:  qs_get_random_integer, qs_sup_norm_matrix,&
       qs_set_to_identity
  use qs_pao_types
  implicit none

  logical, private, parameter :: debug_this_module=.true.

!!***
!****************************************************************************
contains

!!****f* qs_pao_types_tests/qs_test_pao_types [1.0] *
!!
!!   NAME
!!     qs_test_pao_types
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  subroutine qs_test_pao_types(g_env, failed,error)
    type(global_environment_type), intent(in) ::g_env
    logical, intent(inout) :: failed
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types_tests:qs_test_pao_types'
    type(qs_pao_env) :: pao_env
    type(qs_pao_glob_angles) :: glob_angles
    type(qs_pao_local_angles) :: local_angles
    type(qs_pao_matrix) :: pao_t
    integer, parameter :: atomNr=3
    integer, dimension(atomNr) :: minB, fullB
    integer :: i

    failure=.false.
    if (qs_debug.and.debug_this_module) then
       call qs_log(qs_get_error_logger(error),qs_note_level,routineN, &
            'began pao_types tests')
    end if
    forall (i=1:atomNr)
       fullB(i)=qs_get_random_integer(from=1,to=30,error=error)
    end forall
    forall (i=1:atomNr)
       minB(i)=qs_get_random_integer(from=1,to=fullB(i),error=error)
    end forall
    call qs_init(pao_env, nr_min_bas=minB, nr_full_bas=fullB, &
         global_env=g_env,error=error)
    call qs_init(glob_angles, pao_env=pao_env)
    call qs_init(local_angles, glob_angles=glob_angles, &
         atom_nr=qs_get_random_integer(from=1,to=atomNr))
    call qs_init(pao_t,pao_env=pao_env, min_bas_cols=.false.,&
         min_bas_rows=.false.)
    call qs_set(pao_env, subspace_att=pao_t, angles_att=glob_angles)

    QSAssert(qs_validate(pao_env,error=error),qs_warning_level,routineN,error,failure)
    QSAssert(qs_validate(local_angles,error=error),qs_warning_level,routineN,error,failure)
    call qs_dealloc_ref(pao_env,error=error)
    call qs_dealloc_ref(local_angles,error=error)

    call qs_dealloc_ref(glob_angles,error=error)
    call qs_dealloc_ref(pao_t, error=error)
    if (failure) then
       call qs_error_message(qs_warning_level,routineN,&
            "qs_test_pao_types FAILED!",error)
    else if (qs_debug .and. debug_this_module) then
       call qs_error_message(qs_warning_level,routineN,&
            "qs_test_pao_types SUCCESS!",error) ! bug in print utils
       !print *, "qs_test_pao_types SUCCESS!"
    end if
    failed=failure
  end subroutine qs_test_pao_types
!!***
  !***************************************************************************

!!****f* qs_pao_types_tests/qs_init_test_pao_env [1.0] *
!!
!!   NAME
!!     qs_init_test_pao_env
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  subroutine qs_init_test_pao_env(pao_env, glob_env, atom_nr, error)
    type(qs_pao_env), intent(out) :: pao_env
    type(global_environment_type), target, intent(in) :: glob_env
    integer, optional :: atom_nr
    type(qs_error), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_types_tests:qs_init_test_pao_env'
    integer, dimension(:), allocatable :: minB, fullB
    integer :: n_atom,i, stat
    type(qs_pao_matrix), pointer :: u
    nullify(u)
    failure=.false.

    n_atom=3
    if (present(atom_nr)) n_atom=atom_nr
    QSPrecondition(n_atom>0,qs_failure_level,routineN,error,failure)
    if (.not.failure) then
       allocate(minB(n_atom),fullB(n_atom),stat=stat)
       QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
       if (.not. failure) then
          forall (i=1:n_atom)
             fullB(i)=qs_get_random_integer(from=1,to=30,error=error)
          end forall
          forall (i=1:n_atom)
             minB(i)=qs_get_random_integer(from=1,to=fullB(i),error=error)
          end forall
          call qs_init(pao_env, nr_min_bas=minB, nr_full_bas=fullB, &
               global_env=glob_env,error=error)
          call qs_set(pao_env,error=error) ! just to be sure
          call qs_get(pao_env, subspace_att=u,error=error)
          call qs_set_to_identity(qs_get_matrix(u),error=error)
          QSPostcondition(qs_valid(pao_env,error=error),qs_warning_level,routineN,error,failure)
       end if
       deallocate(minB,fullB,stat=stat) !just to be sure...
       QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
    end if
  end subroutine qs_init_test_pao_env
!!***
  !***************************************************************************

end module qs_pao_types_tests


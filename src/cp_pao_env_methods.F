!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_env_methods [1.0] *
!!
!!   NAME
!!     cp_pao_env_methods
!!
!!   FUNCTION
!!     methods to access and modify the pao environement
!!
!!   NOTES
!!     the get methods are defined in cp_pao_types, because of the circular
!!     dependencies in the types
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_pao_env_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use atoms, only: atom_info, atom_info_type
  use atomic_kinds, only: kind_info,nkind
  use cp_pao_types
  use cp_pao_utils
  use cp_pao_matrix_methods
  use cp_pao_glob_angles_methods
  use matrix_types
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  logical, private, parameter :: standalone=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_env_methods'

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set, cp_did_change

! core functions
  public :: cp_validate

! underlying routines
  public :: cp_pao_env_init, cp_pao_env_dealloc_ref, cp_pao_env_validate,&
       cp_pao_env_set, cp_pao_env_did_change

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure cp_pao_env_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure cp_pao_env_dealloc_ref
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure cp_pao_env_validate
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure cp_pao_env_set
  end interface

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  interface cp_did_change
     module procedure cp_pao_env_did_change
  end interface

! +++ end of the core methods +++

!!***
!****************************************************************************

contains

! ====== core methods =====

!!****f* cp_pao_env_methods/cp_pao_env_init [1.0] *
!!
!!   NAME
!!     cp_pao_env_init
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_init(pao_env, nr_min_bas, nr_full_bas,&
!!         global_env, full_overlap_sm, full_hamiltonian_sm,&
!!         wanted_subspace_pos_precision, wanted_in_subspace_precision,&
!!         projection, error)
!!       Type(cp_pao_env_type), Intent (OUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Integer, Dimension(:), Intent (IN):: nr_min_bas, nr_full_bas
!!       Type(real_matrix_type), Target, Intent (IN), Optional::&
!!         full_hamiltonian_sm, full_overlap_sm
!!       Real(Kind=wp), Optional, Intent (IN)::&
!!         wanted_in_subspace_precision, wanted_subspace_pos_precision
!!       Type(cp_1d_i_array), Dimension(:), Optional, Intent (IN)::&
!!         projection
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_init
!!
!!   FUNCTION
!!     initialize the pao environement
!!
!!   NOTES
!!     every environement HAS to be initialized with this routine before 
!!     using.
!!     Unless said all arguments are COPIED if nedeed, not shared.
!!     global_env, is SHARED.
!!
!!   INPUTS
!!     pao_env: the environement that has to be initialized
!!     projection: the projection, if given, the projection subarrays 
!!          are shared
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     the other arguments are actually a shortcut not to call
!!     cp_pao_env_set, for their description see cp_pao_env_set.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_env_init(pao_env,nr_min_bas, nr_full_bas,global_env,&
       full_overlap_sm, full_hamiltonian_sm,&
       wanted_subspace_pos_precision,&
       wanted_in_subspace_precision,&
       projection, error)
    type(cp_pao_env_type), intent(out) :: pao_env
    type(global_environment_type), intent(in), target :: global_env
    integer, dimension(:), intent(in) :: nr_min_bas, nr_full_bas
    type(real_matrix_type), target, intent(in), optional :: &
         full_hamiltonian_sm, full_overlap_sm
    real(kind=wp), optional, intent(in) :: wanted_in_subspace_precision, &
         wanted_subspace_pos_precision
    type(cp_1d_i_array), dimension(:), optional, intent(in) :: projection
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_env_methods:cp_pao_env_init'
    type(cp_error_type) :: iError
    integer :: stat, i
    failure=.false.
    call cp_init(iError, template_error=error)

    nullify(pao_env%nr_min_bas, pao_env%global_env, pao_env%nr_full_bas)
    nullify(pao_env%projection, pao_env%angles_begin_at)
    pao_env%min_d_m_valid=.false.
    pao_env%n_atoms=-1
! global_env
    pao_env%global_env => global_env
    call cp_assert(associated(pao_env%global_env),cp_warning_level,&
         cp_assertion_failed,routineP,"glob env not associated in "//&
         CPSourceFileRef,error)
! min_bas
    call cp_guarantee_size(pao_env%nr_min_bas,size(nr_min_bas), &
         error=iError)
    CPPostcondition(cp_get_error_level(iError)<cp_failure_level,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       pao_env%nr_min_bas=nr_min_bas
    end if
! full_bas
    call cp_guarantee_size(pao_env%nr_full_bas,size(nr_full_bas), &
         error=iError)
    CPPostcondition(cp_get_error_level(iError)<cp_failure_level,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       pao_env%nr_full_bas=nr_full_bas
    end if
    call cp_propagate_error(iError,routineP,&
         CPSourceFileRef,&
         error,failure,cp_failure_level)
! check min_bas/full_bas
    CPPrecondition(associated(pao_env%nr_min_bas),cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(pao_env%nr_full_bas),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_assert(size(pao_env%nr_full_bas)==size(pao_env%nr_min_bas),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "min_bas and full bas must have the same size in "//&
            CPSourceFileRef,&
            error,failure)
    end if
! angles_begin_at
    if (.not.failure) then
       pao_env%n_atoms=size(pao_env%nr_full_bas)
       allocate(pao_env%angles_begin_at(pao_env%n_atoms+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       if (pao_env%n_atoms>0) pao_env%angles_begin_at(1)=0
       do i=1,pao_env%n_atoms
          pao_env%angles_begin_at(i+1)=pao_env%angles_begin_at(i)+&
               (pao_env%nr_full_bas(i)-pao_env%nr_min_bas(i))*&
               pao_env%nr_min_bas(i)
       end do
    end if
! others
    pao_env%wanted_subspace_pos_precision=1.0e-4_wp
    pao_env%wanted_in_subspace_precision=1.0e-4_wp
    pao_env%initialized=.true.
! projection
    if (present(projection)) then
       allocate(pao_env%projection(size(projection)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          pao_env%projection=projection
       end if
       pao_env%dealloc_projections=.false.
    else
       pao_env%dealloc_projections=.true.
       call cp_pao_env_init_projection(pao_env,error)
    end if
! n_projection
! call cp_pao_env_init_n_projection(pao_env,error)
    allocate(pao_env%min_density_m,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    allocate(pao_env%full_overlap_m,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    allocate(pao_env%full_hamiltonian_m,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    allocate(pao_env%angles_att,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
! full_h
    if (.not.failure) then
       call cp_init(pao_env%full_hamiltonian_m, matrix=full_hamiltonian_sm, &
            min_bas_cols=.false., min_bas_rows=.false., pao_env=pao_env,&
            should_dealloc_matrix=.not.present(full_hamiltonian_sm), error=error)
    end if
! full_s
    if (.not.failure) then
       call cp_init(pao_env%full_overlap_m, matrix=full_overlap_sm,&
            min_bas_cols=.false., min_bas_rows=.false., pao_env=pao_env,&
            should_dealloc_matrix=.not.present(full_overlap_sm), error=error)
    end if
! min_d
    if (.not.failure) then
       call cp_init(pao_env%min_density_m, min_bas_rows=.true.,&
            min_bas_cols=.true., pao_env=pao_env,&
            should_dealloc_matrix=.true., error=error)
    end if
! angles_att
    if (.not.failure) then
       call cp_init(pao_env%angles_att, pao_env=pao_env,&
            error=error)
    end if
    call cp_pao_env_set(pao_env,&
         wanted_subspace_pos_precision=wanted_subspace_pos_precision,&
         wanted_in_subspace_precision=wanted_in_subspace_precision,&
         error=error)
    if (cp_debug.and.debug_this_module) then
       call cp_assert(cp_validate(pao_env,test_pos_att=.false.,error=error),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "pao_env is not valid after init in "//&
            CPSourceFileRef,&
            error,failure)
       call cp_assert(cp_valid(pao_env,error=error),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "pao_env is not valid after init in "//&
            CPSourceFileRef,&
            error,failure)
    end if
    call cp_dealloc_ref(iError,error)
  end subroutine cp_pao_env_init
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_pao_env_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_dealloc_ref(pao_env, error)
!!       Type(cp_pao_env_type), Intent (INOUT):: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_dealloc_ref
!!
!!   FUNCTION
!!     deallocate the memory that the pao environement has allocated,
!!     but NOT the pao environement itself.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the environement to be deallocated
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
  subroutine cp_pao_env_dealloc_ref(pao_env, error)
    type(cp_pao_env_type), intent(inout) :: pao_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_env_methods:cp_pao_env_dealloc_ref'
    integer :: stat,i
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(pao_env%angles_att)) then
          call cp_dealloc(pao_env%angles_att, error=error)
       end if
       call cp_dealloc(pao_env%min_density_m,&
            error=error)
       call cp_dealloc(pao_env%full_hamiltonian_m,&
            error=error)
       call cp_dealloc(pao_env%full_overlap_m,&
            error=error)
       nullify(pao_env%global_env)
       if (associated(pao_env%projection)) then
          if (pao_env%dealloc_projections) then
             do i=1,size(pao_env%projection)
                if (associated(pao_env%projection(i)%array)) then
                   deallocate(pao_env%projection(i)%array, stat=stat)
                   CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
                end if
             end do
          end if
          deallocate(pao_env%projection,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       if (associated(pao_env%nr_min_bas)) then
          deallocate(pao_env%nr_min_bas, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure) 
       end if
       if (associated(pao_env%nr_full_bas)) then
          deallocate(pao_env%nr_full_bas, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure) 
       end if
       if (associated(pao_env%angles_begin_at)) then
          deallocate(pao_env%angles_begin_at,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       pao_env%initialized=.false.
    end if failureIf
  end subroutine cp_pao_env_dealloc_ref
!***************************************************************************


!!****f* cp_pao_env_methods/cp_pao_env_validate [1.0] *
!!
!!   NAME
!!     cp_pao_env_validate
!!
!!   SYNOPSIS
!!     Function cp_pao_env_validate(pao_env, test_pos_att, error)&
!!         Result(res)
!!       Logical:: res
!!       Type(cp_pao_env_type), Intent (IN):: pao_env
!!       Logical, Optional:: test_pos_att
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_env_validate
!!
!!   FUNCTION
!!     returns true if the pao environement is valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     pao_env: the environement to be tested
!!     test_pos_att: if it is false does not test the validity of the
!!          actual angles and unitary transformation. defaults to true
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_env_validate(pao_env, test_pos_att, error) result(res)
    logical ::res
    type(cp_pao_env_type), intent(in) :: pao_env
    logical ,optional ::test_pos_att
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, test_pos
    character(len=*), parameter :: routineP='cp_pao_env_methods:cp_pao_env_validate'
    type(cp_logger), pointer :: logger
    integer :: i,m_bas_att,f_bas_att,j
    failure=.false.
    logger => cp_get_error_logger(error)
    if (present(test_pos_att)) then
       test_pos=test_pos_att
    else
       test_pos=.true.
    end if
    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       CPAssert(associated(pao_env%global_env),cp_warning_level,routineP,error,failure)
       CPAssert(pao_env%wanted_subspace_pos_precision<1.0e-2_wp,cp_warning_level,routineP,error,failure)
       CPAssert(pao_env%wanted_in_subspace_precision<1.0e-2_wp,cp_warning_level,routineP,error,failure)
       CPAssert(associated(pao_env%nr_min_bas), cp_warning_level,routineP,error,failure)
       CPAssert(associated(pao_env%nr_full_bas), cp_warning_level,routineP,error,failure)
       CPAssert(size(pao_env%nr_min_bas)>0,cp_warning_level,routineP,error,failure)
       call cp_assert(size(pao_env%nr_min_bas)==size(pao_env%nr_full_bas),&
            cp_warning_level,cp_assertion_failed,routineP,&
            "min_bas and full_s must have the same length in"//&
            CPSourceFileRef,&
            error,failure)
       CPAssert(size(pao_env%nr_min_bas)==pao_env%n_atoms,cp_warning_level,routineP,error,failure)
       if (.not. standalone) then
          CPAssert(size(pao_env%nr_min_bas)==size(atom_info),cp_warning_level,routineP,error,failure)
       end if
       failureIf2: if (.not.failure .and.(cp_debug.or.test_pos)) then 
! I could use all() and any(), but then the info on the atom nr 
! would be lost...
          do i=1,size(pao_env%nr_min_bas)
             m_bas_att=pao_env%nr_min_bas(i)
             f_bas_att=pao_env%nr_full_bas(i)
             call cp_assert (m_bas_att >= 1, cp_warning_level, &
                  cp_assertion_failed, routineP,&
                  'the atom '//cp_to_string(i)// &
                  ' has a minmal basis that is less than 1', error,failure)
             call cp_assert(m_bas_att<= f_bas_att,cp_warning_level, &
                  cp_assertion_failed, routineP,&
                  ' min bas>number of orbitals for atom'//cp_to_string(i),&
                  error,failure)
             if (.not.standalone) then
                call cp_assert( f_bas_att == &
                     atom_info(i)%last_sgf - atom_info(i)%first_sgf +1, &
                     cp_warning_level,cp_assertion_failed, routineP, &
                     'full bas is different from hamiltonian basis on atom '//&
                     cp_to_string(i),error,failure)
             end if
          end do
       end if failureIf2
       CPAssert(associated(pao_env%projection),cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          do i=1,size(pao_env%projection)
             CPAssert(associated(pao_env%projection(i)%array),cp_warning_level,routineP,error,failure)
             if (.not.failure) then
                call cp_assert(size(pao_env%projection(i)%array)==&
                     pao_env%nr_min_bas(i),cp_warning_level,&
                     cp_assertion_failed,routineP,&
                     "projection has invalid size in "//&
                     CPSourceFileRef,&
                     error,failure)
                CPAssert(all(pao_env%projection(i)%array>0),cp_warning_level,routineP,error,failure)
                do j=2,size(pao_env%projection(i)%array)
                   call cp_assert(pao_env%projection(i)%array(j) > &
                        pao_env%projection(i)%array(j-1),cp_warning_level,&
                        cp_assertion_failed, routineP, &
                        "projection indices are not ordered in "//&
                        CPSourceFileRef,&
                        error,failure)
                end do
             end if
          end do
       end if
! angles_begin_at
       CPPrecondition(associated(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(pao_env%n_atoms+1==size(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_assert(all((pao_env%angles_begin_at(2:pao_env%n_atoms+1)-&
               pao_env%angles_begin_at(1:pao_env%n_atoms))==&
               pao_env%nr_min_bas*(pao_env%nr_full_bas-pao_env%nr_min_bas)),&
               cp_warning_level,cp_assertion_failed,routineP,&
               "angles_begin_at is not compatible with nr_min_bas nr_full_bas in "//&
               CPSourceFileRef,&
               error,failure)
       end if
! pos att
       if (test_pos) then
          if (.not.cp_validate(pao_env%angles_att,&
               error=error)) failure=.true.
          if (.not.cp_validate(pao_env%min_density_m,&
               error=error)) failure=.true. ! also test idempotency ?
       end if
       if (.not.standalone) then
          if (.not.cp_validate(pao_env%full_overlap_m,&
               error=error)) failure=.true. ! should be positive definite
          if (.not.cp_validate(pao_env%full_hamiltonian_m,&
               error=error)) failure=.true. ! should be hermitian
       end if
    end if failureIf
    res=.not.failure
  end function cp_pao_env_validate
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_set [1.0] *
!!
!!   NAME
!!     cp_pao_env_set
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_set(pao_env, wanted_subspace_pos_precision,&
!!         wanted_in_subspace_precision, angles_att, min_density_m,&
!!         full_overlap_m, full_hamiltonian_m, min_d_valid, error)
!!       Type(cp_pao_env_type), Target, Intent (INOUT):: pao_env
!!       Real(Kind=wp), Optional, Intent (IN)::&
!!         wanted_in_subspace_precision, wanted_subspace_pos_precision
!!       Type(cp_pao_glob_angles_type), Optional, Intent (IN):: angles_att
!!       Type(cp_pao_matrix_type), Optional, Intent (IN):: min_density_m,&
!!         full_hamiltonian_m, full_overlap_m
!!       Logical, Optional:: min_d_valid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_set
!!
!!   FUNCTION
!!     changes some of the parameters of the environement
!!
!!   NOTES
!!     arguments, unless otherwise stated are always copied
!!
!!   INPUTS
!!     pao_env: the environement that has to be initialized
!!     global_env: the global environement (for the future), optional,
!!                 it is SHARED
!!     n_atoms: the number of atoms (used to initialize some structures)
!!              optional, it is inferred from nr_min_bas or nr_full_bas
!!              they are present
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     for the description of the other variables see the attributes 
!!     of cp_pao_env_type
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_env_set(pao_env, wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, angles_att, &
       min_density_m,full_overlap_m, full_hamiltonian_m,&
       min_d_valid,error)
    type(cp_pao_env_type), target, intent(inout) :: pao_env
    real(kind=wp), optional, intent(in) :: wanted_in_subspace_precision, &
         wanted_subspace_pos_precision
    type(cp_pao_glob_angles_type),optional, intent(in) :: angles_att
    type(cp_pao_matrix_type),optional, intent(in) :: &
         min_density_m, full_hamiltonian_m, full_overlap_m
    logical, optional :: min_d_valid
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_env_methods:cp_pao_env_set'
    integer :: stat
    type(cp_error_type) :: iError
    failure=.false.
    call cp_init_error(iError,template_error=error)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (present(wanted_subspace_pos_precision)) &
         pao_env%wanted_subspace_pos_precision=wanted_subspace_pos_precision
    if (present(wanted_in_subspace_precision)) &
         pao_env%wanted_in_subspace_precision=wanted_in_subspace_precision
    if (present(angles_att)) then
       call cp_set(pao_env%angles_att, source=angles_att, pao_env=pao_env,&
            error=error)
    end if
    if (present(min_density_m)) then
       pao_env%min_d_m_valid=.true.
       call cp_set(pao_env%min_density_m, source=min_density_m,&
            pao_env=pao_env, error=error)
    end if
    if (present(full_hamiltonian_m)) then
       call cp_set(pao_env%full_hamiltonian_m, source=full_hamiltonian_m,&
            pao_env=pao_env, error=error)
    end if
    if (present(full_overlap_m)) then
       call cp_set(pao_env%full_overlap_m, source=full_overlap_m,&
            pao_env=pao_env, error=error)
    end if
    if (present(min_d_valid)) then
       pao_env%min_d_m_valid=min_d_valid
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_env_set
!***************************************************************************

!!****f* cp_pao_env_methods/cp_pao_env_did_change [1.0] *
!!
!!   NAME
!!     cp_pao_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_did_change(pao_env, h_changed, s_changed,&
!!         angles_changed, full_reset, error)
!!       Type(cp_pao_env_type), Intent (INOUT):: pao_env
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_did_change
!!
!!   FUNCTION
!!     invalidates the cache when some data changes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
subroutine cp_pao_env_did_change(pao_env,h_changed, s_changed,&
     angles_changed, full_reset, error)
  type (cp_pao_env_type), intent(inout) :: pao_env
  logical, intent(in), optional :: h_changed, s_changed, &
       angles_changed, full_reset
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle
  character(len=*), parameter :: routineN='cp_pao_env_did_change',&
       routineP=moduleN//':'//routineN
  logical :: h_c, s_c, a_c, full_r
  !call timeset(routineN,'I',moduleN,handle)
  failure=.false.
  h_c=.true. ; s_c=.true. ; a_c=.true. ; full_r=.false.
  
  CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(h_changed)) h_c=h_changed
     if (present(s_changed)) s_c=s_changed
     if (present(angles_changed)) a_c=angles_changed
     if (present(full_reset)) full_r=full_reset
     call cp_did_change(pao_env%angles_att, h_changed, s_changed,&
          angles_changed, full_reset, error)
     
  end if
  !call timestop(0.0_wp,handle)
end subroutine cp_pao_env_did_change
!***************************************************************************

! ====== special methods ======


! ========= private methods ========

!!****f* cp_pao_env_methods/cp_pao_env_init_projection [1.0] *
!!
!!   NAME
!!     cp_pao_env_init_projection
!!
!!   FUNCTION
!!     checks the pao_env projection and if not valid initializes it
!!
!!   NOTES
!!     private method
!!
!!   INPUTS
!!     pao_env: the pao environement where check the projection
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
  subroutine cp_pao_env_init_projection(pao_env, error)
    type(cp_pao_env_type), intent(inout) :: pao_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_env_methods:cp_pao_env_init_projection'
    integer :: natom, i,stat,j
    failure=.false.
    natom=cp_get_n_atom(pao_env)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(natom>0,cp_warning_level,routineP,error,failure)
    CPPrecondition(associated(pao_env%nr_min_bas),cp_warning_level,routineP,error,failure)
    if (.not.failure) then
       CPPrecondition(size(pao_env%nr_min_bas)==natom,cp_warning_level,routineP,error,failure)
    end if
    failureIf: if (.not. failure) then
       if (.not.associated(pao_env%projection)) then
          allocate(pao_env%projection(natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             do i=1,natom
                nullify(pao_env%projection(i)%array)
             end do
          end if
       end if
       CPAssert(size(pao_env%projection)==cp_get_n_atom(pao_env),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          do i=1,natom
             if (associated(pao_env%projection(i)%array)) then
                if (size(pao_env%projection(i)%array) /= pao_env%nr_min_bas(i)) then
                   call cp_error_message(cp_warning_level,&
                        "minbas size changed on atom "// cp_to_string(i)//&
                        CPSourceFileRef,&
                        routineP,error=error)
                   deallocate(pao_env%projection(i)%array ,stat=stat)
                   CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
                end if
             end if
             if (.not.associated(pao_env%projection(i)%array)) then
                allocate(pao_env%projection(i)%array(pao_env%nr_min_bas(i)),&
                     stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                if (.not.failure) then
                   do j=1,pao_env%nr_min_bas(i)
                      pao_env%projection(i)%array(j)=j
                   end do
                end if
             end if
          end do
       end if
    end if failureIf
  end subroutine cp_pao_env_init_projection
!***************************************************************************

end module cp_pao_env_methods

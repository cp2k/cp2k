!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_env_methods [1.0] *
!!
!!   NAME
!!     pao_env_methods
!!
!!   FUNCTION
!!     methods to access and modify the pao environement
!!
!!   NOTES
!!     the get methods are defined in pao_types, because of the circular
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
module pao_env_methods
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_assert, &
       cp_error_message, cp_assertion_failed, cp_error_init,&
       cp_error_get_print_level, cp_error_get_stop_level, cp_error_get_level,&
       cp_dealloc_ref, cp_error_get_logger, cp_internal_error
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types, only: pao_glob_angles_type, pao_env_type, cp_get, cp_valid
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use cp_sparse_matrix, only: cp_sparse_matrix_type, cp_sparse_matrix_p_type,&
       cp_create, cp_dealloc, cp_get_matrix, cp_retain, cp_release
  use pao_glob_angles_methods, only: pao_g_ang_init, cp_dealloc, cp_set
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type
  use pao_bas_projection_types, only: pao_bas_projection_type,cp_proj_get
  use pao_qs_env_methods, only: cp_get_b_projection
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type, &
       cp_bmstruct_create, cp_bmstruct_release
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  logical, private, parameter :: standalone=.true.
  character(len=*), private, parameter :: moduleN='pao_env_methods'

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set

! core functions
  public :: cp_validate

! underlying routines
  public :: pao_env_init, pao_env_dealloc_ref, pao_env_validate,&
       pao_env_set

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure pao_env_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure pao_env_dealloc_ref
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure pao_env_validate
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure pao_env_set
  end interface

! +++ end of the core methods +++

!!***
!****************************************************************************

contains

! ====== core methods =====

!!****f* pao_env_methods/pao_env_init [1.0] *
!!
!!   NAME
!!     pao_env_init
!!
!!   SYNOPSIS
!!     Subroutine pao_env_init(pao_env, qs_env, global_env,&
!!         wanted_subspace_pos_precision, wanted_in_subspace_precision,&
!!         error)
!!       Type(pao_env_type), Intent (OUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Real(Kind=wp), Optional, Intent (IN)::&
!!         wanted_in_subspace_precision, wanted_subspace_pos_precision
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_init
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
!!     pao_env_set, for their description see pao_env_set.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_env_init(pao_env, qs_env, global_env,&
       wanted_subspace_pos_precision,&
       wanted_in_subspace_precision,&
       error)
    type(pao_env_type), intent(out) :: pao_env
    type(global_environment_type), intent(in), target :: global_env
    type(qs_environment_type), intent(inout), target :: qs_env
    real(kind=wp), optional, intent(in) :: wanted_in_subspace_precision, &
         wanted_subspace_pos_precision
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_env_methods:pao_env_init'
    type(cp_error_type) :: iError
    integer :: stat, i, n_atoms, size_excl_bas, size_min_bas, size_full_bas,j
    type(real_matrix_p_type), dimension(:), pointer :: p
    type(pao_bas_projection_type), pointer :: bproj
    integer, dimension(:), allocatable :: min_bas_sizes, full_bas_sizes
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    integer, dimension(:), pointer :: atom_list
    character(len=60) :: name
    failure=.false.
    call cp_error_init(iError, template_error=error)

    nullify(pao_env%angles_begin_at)
    pao_env%min_d_m_valid=.false.

    ! check projection
    call get_qs_env(qs_env, atomic_kind_set= atomic_kind_set,natoms=n_atoms)
    CPPrecondition(associated(atomic_kind_set),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       allocate(min_bas_sizes(n_atoms), full_bas_sizes(n_atoms),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       do i=1,size(atomic_kind_set)
          call get_atomic_kind(atomic_kind_set(i),pao_bas_proj=bproj,&
               atom_list=atom_list)
          if (.not.associated(bproj)) then
             call get_atomic_kind(atomic_kind_set(i),name=name)
             call cp_assert(.false.,cp_failure_level, cp_internal_error,&
                  routineP,&
                  "projection (PAO_MIN_BAS) not defined in kind "//name &
                  //", "//&
                  CPSourceFileRef,&
                  error=error,failure=failure)
          else
             call cp_proj_get(bproj,size_min_bas=size_min_bas,&
                  size_full_bas=size_full_bas,error=error)
             do j=1,size(atom_list)
                min_bas_sizes(atom_list(i))=size_min_bas
                full_bas_sizes(atom_list(i))=size_full_bas
             end do
          end if
       end do
       
       ! min_bas_m_struct / full_bas_m_struct
       if (.not.failure) then
          call cp_bmstruct_create(pao_env%min_bas_m_struct,&
               symmetric=.true.,&
               global_env=global_env,&
               b_row_sizes=min_bas_sizes,&
               b_col_sizes=min_bas_sizes,&
               error=error)
          call cp_bmstruct_create(pao_env%full_bas_m_struct,&
               symmetric=.true.,&
               global_env=global_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
       end if
       deallocate(min_bas_sizes, full_bas_sizes, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if

    ! others
    call get_qs_env(qs_env,p=p, natoms=n_atoms)
    allocate(pao_env%min_density_m(size(p)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       do i=1,size(p)
          nullify(pao_env%min_density_m(i)%matrix)
       end do
    end if
    allocate(pao_env%angles_att,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    pao_env%wanted_subspace_pos_precision=1.0e-4_wp
    pao_env%wanted_in_subspace_precision=1.0e-4_wp
    pao_env%initialized=.true.
    ! min_d
    if (.not.failure) then
       do i=1,size(pao_env%min_density_m)
          call cp_create(pao_env%min_density_m(i)%matrix,&
               matrix_structure=pao_env%min_bas_m_struct,&
               global_env=global_env, error=error)
       end do
    end if
    ! angles_att
    if (.not.failure) then
       call pao_g_ang_init(pao_env%angles_att, global_env=global_env,&
            qs_env=qs_env, error=error)
    end if
    ! angles_begin_at
    if (.not.failure) then
       allocate(pao_env%angles_begin_at(n_atoms+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       if (n_atoms>0) pao_env%angles_begin_at(1)=0
       do i=1,n_atoms
          bproj => cp_get_b_projection(qs_env, atom_nr=i, error=error)
          call cp_get(bproj, size_excl_bas=size_excl_bas, &
               size_min_bas=size_min_bas,error=error)
          pao_env%angles_begin_at(i+1)=pao_env%angles_begin_at(i)+&
               size_excl_bas*size_min_bas
       end do
    end if

    call pao_env_set(pao_env,qs_env=qs_env, global_env=global_env,&
         wanted_subspace_pos_precision=wanted_subspace_pos_precision,&
         wanted_in_subspace_precision=wanted_in_subspace_precision,&
         error=error)
    if (cp_debug.and.debug_this_module) then
       call cp_assert(cp_validate(pao_env,qs_env=qs_env,&
            global_env=global_env,&
            test_pos_att=.false.,error=error),&
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
  end subroutine pao_env_init
!***************************************************************************

!!****f* pao_env_methods/pao_env_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_env_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_env_dealloc_ref(pao_env, qs_env, global_env, error)
!!       Type(pao_env_type), Intent (INOUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_dealloc_ref
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
  subroutine pao_env_dealloc_ref(pao_env,qs_env,global_env, error)
    type(pao_env_type), intent(inout) :: pao_env
    type(global_environment_type), intent(in), target :: global_env
    type(qs_environment_type), intent(inout), target :: qs_env    
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_env_methods:pao_env_dealloc_ref'
    integer :: stat,i
    failure=.false.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(pao_env%angles_att)) then
          call cp_dealloc(pao_env%angles_att,global_env=global_env,&
               error=error)
       end if
       if (associated(pao_env%angles_begin_at)) then
          deallocate(pao_env%angles_begin_at,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       if (associated(pao_env%min_bas_m_struct)) then
          call cp_bmstruct_release(pao_env%min_bas_m_struct,error=error)
          nullify(pao_env%min_bas_m_struct)
       end if
       if (associated(pao_env%full_bas_m_struct)) then
          call cp_bmstruct_release(pao_env%full_bas_m_struct, error=error)
          nullify(pao_env%full_bas_m_struct)
       end if
       pao_env%initialized=.false.
    end if failureIf
  end subroutine pao_env_dealloc_ref
!***************************************************************************


!!****f* pao_env_methods/pao_env_validate [1.0] *
!!
!!   NAME
!!     pao_env_validate
!!
!!   SYNOPSIS
!!     Function pao_env_validate(pao_env, qs_env, global_env, test_pos_att,&
!!         error) Result(res)
!!       Logical:: res
!!       Type(pao_env_type), Intent (IN), Target:: pao_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Logical, Optional:: test_pos_att
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_env_validate
!!
!!   FUNCTION
!!     returns true if the pao environement is valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     qs_env: the environement to be tested
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
  function pao_env_validate(pao_env,qs_env,global_env, test_pos_att, error)&
       result(res)
    logical ::res
    type(pao_env_type), intent(in), target :: pao_env
    type(qs_environment_type), intent(in),target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    logical ,optional ::test_pos_att
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, test_pos
    character(len=*), parameter :: routineP='pao_env_methods:pao_env_validate'
    type(cp_logger_type), pointer :: logger
    integer :: i,m_bas_att,f_bas_att,j
    failure=.false.
    logger => cp_error_get_logger(error)
    if (present(test_pos_att)) then
       test_pos=test_pos_att
    else
       test_pos=.true.
    end if
    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
! to do
!!FM    failureIf: if (.not. failure) then
!!FM       CPAssert(associated(pao_env%global_env),cp_warning_level,routineP,error,failure)
!!FM       CPAssert(pao_env%wanted_subspace_pos_precision<1.0e-2_wp,cp_warning_level,routineP,error,failure)
!!FM       CPAssert(pao_env%wanted_in_subspace_precision<1.0e-2_wp,cp_warning_level,routineP,error,failure)
!!FM       CPAssert(associated(pao_env%nr_min_bas), cp_warning_level,routineP,error,failure)
!!FM       CPAssert(associated(pao_env%nr_full_bas), cp_warning_level,routineP,error,failure)
!!FM       CPAssert(size(pao_env%nr_min_bas)>0,cp_warning_level,routineP,error,failure)
!!FM       call cp_assert(size(pao_env%nr_min_bas)==size(pao_env%nr_full_bas),&
!!FM            cp_warning_level,cp_assertion_failed,routineP,&
!!FM            "min_bas and full_s must have the same length in"//&
!!FM            CPSourceFileRef,&
!!FM            error,failure)
!!FM       CPAssert(size(pao_env%nr_min_bas)==pao_env%n_atoms,cp_warning_level,routineP,error,failure)
!!FM       if (.not. standalone) then
!!FM          !CPAssert(size(pao_env%nr_min_bas)==size(atom_info),cp_warning_level,routineP,error,failure)
!!FM       end if
!!FM       failureIf2: if (.not.failure .and.(cp_debug.or.test_pos)) then 
!!FM! I could use all() and any(), but then the info on the atom nr 
!!FM! would be lost...
!!FM          do i=1,size(pao_env%nr_min_bas)
!!FM             m_bas_att=pao_env%nr_min_bas(i)
!!FM             f_bas_att=pao_env%nr_full_bas(i)
!!FM             call cp_assert (m_bas_att >= 1, cp_warning_level, &
!!FM                  cp_assertion_failed, routineP,&
!!FM                  'the atom '//cp_to_string(i)// &
!!FM                  ' has a minmal basis that is less than 1', error,failure)
!!FM             call cp_assert(m_bas_att<= f_bas_att,cp_warning_level, &
!!FM                  cp_assertion_failed, routineP,&
!!FM                  ' min bas>number of orbitals for atom'//cp_to_string(i),&
!!FM                  error,failure)
!!FM             if (.not.standalone) then
!!FM                !call cp_assert( f_bas_att == &
!!FM                !     atom_info(i)%last_sgf - atom_info(i)%first_sgf +1, &
!!FM                !     cp_warning_level,cp_assertion_failed, routineP, &
!!FM                !     'full bas is different from hamiltonian basis on atom '//&
!!FM                !     cp_to_string(i),error,failure)
!!FM             end if
!!FM          end do
!!FM       end if failureIf2
!!FM       CPAssert(associated(pao_env%projection),cp_warning_level,routineP,error,failure)
!!FM       if (.not.failure) then
!!FM          do i=1,size(pao_env%projection)
!!FM             CPAssert(associated(pao_env%projection(i)%array),cp_warning_level,routineP,error,failure)
!!FM             if (.not.failure) then
!!FM                call cp_assert(size(pao_env%projection(i)%array)==&
!!FM                     pao_env%nr_min_bas(i),cp_warning_level,&
!!FM                     cp_assertion_failed,routineP,&
!!FM                     "projection has invalid size in "//&
!!FM                     CPSourceFileRef,&
!!FM                     error,failure)
!!FM                CPAssert(all(pao_env%projection(i)%array>0),cp_warning_level,routineP,error,failure)
!!FM                do j=2,size(pao_env%projection(i)%array)
!!FM                   call cp_assert(pao_env%projection(i)%array(j) > &
!!FM                        pao_env%projection(i)%array(j-1),cp_warning_level,&
!!FM                        cp_assertion_failed, routineP, &
!!FM                        "projection indexes are not ordered in "//&
!!FM                        CPSourceFileRef,&
!!FM                        error,failure)
!!FM                end do
!!FM             end if
!!FM          end do
!!FM       end if
!!FM! angles_begin_at
!!FM       CPPrecondition(associated(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
!!FM       if (.not.failure) then
!!FM          CPPrecondition(pao_env%n_atoms+1==size(pao_env%angles_begin_at),cp_warning_level,routineP,error,failure)
!!FM       end if
!!FM       if (.not.failure) then
!!FM          call cp_assert(all((pao_env%angles_begin_at(2:pao_env%n_atoms+1)-&
!!FM               pao_env%angles_begin_at(1:pao_env%n_atoms))==&
!!FM               pao_env%nr_min_bas*(pao_env%nr_full_bas-pao_env%nr_min_bas)),&
!!FM               cp_warning_level,cp_assertion_failed,routineP,&
!!FM               "angles_begin_at is not compatible with nr_min_bas nr_full_bas in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
!!FM       end if
!!FM! pos att
!!FM       if (test_pos) then
!!FM          if (.not.cp_validate(pao_env%angles_att,&
!!FM               error=error)) failure=.true.
!!FM          if (.not.cp_validate(pao_env%min_density_m,&
!!FM               error=error)) failure=.true. ! also test idempotency ?
!!FM       end if
!!FM       if (.not.standalone) then
!!FM          if (.not.cp_validate(pao_env%full_overlap_m,&
!!FM               error=error)) failure=.true. ! should be positive definite
!!FM          if (.not.cp_validate(pao_env%full_hamiltonian_m,&
!!FM               error=error)) failure=.true. ! should be hermitian
!!FM       end if
!!FM    end if failureIf
    res=.not.failure
  end function pao_env_validate
!***************************************************************************

!!****f* pao_env_methods/pao_env_set [1.0] *
!!
!!   NAME
!!     pao_env_set
!!
!!   SYNOPSIS
!!     Subroutine pao_env_set(pao_env, qs_env, global_env,&
!!         wanted_subspace_pos_precision, wanted_in_subspace_precision,&
!!         angles_att, min_density_m, min_d_valid, error)
!!       Type(pao_env_type), Target, Intent (INOUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Real(Kind=wp), Optional, Intent (IN)::&
!!         wanted_in_subspace_precision, wanted_subspace_pos_precision
!!       Type(pao_glob_angles_type), Optional, Intent (IN):: angles_att
!!       Type(cp_sparse_matrix_p_type), Optional, Intent (IN), Dimension(:),&
!!         Target:: min_density_m
!!       Logical, Optional:: min_d_valid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_set
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
!!     min_density_m: the density in the minimal basis (SHARED)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     for the description of the other variables see the attributes 
!!     of pao_env_type
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_env_set(pao_env, qs_env, global_env, &
       wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, angles_att, &
       min_density_m, min_d_valid,error)
    type(pao_env_type), target, intent(inout) :: pao_env
    type(global_environment_type), intent(in), target :: global_env
    type(qs_environment_type), intent(in), target :: qs_env
    real(kind=wp), optional, intent(in) :: wanted_in_subspace_precision, &
         wanted_subspace_pos_precision
    type(pao_glob_angles_type),optional, intent(in) :: angles_att
    type(cp_sparse_matrix_p_type),optional, intent(in), dimension(:),&
         target :: min_density_m
    logical, optional :: min_d_valid
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_env_methods:pao_env_set'
    integer :: stat
    type(cp_error_type) :: iError
    failure=.false.
    call cp_error_init(iError,template_error=error)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    if (present(wanted_subspace_pos_precision)) &
         pao_env%wanted_subspace_pos_precision=wanted_subspace_pos_precision
    if (present(wanted_in_subspace_precision)) &
         pao_env%wanted_in_subspace_precision=wanted_in_subspace_precision
    if (present(angles_att)) then
       call cp_set(pao_env%angles_att, qs_env=qs_env, global_env=global_env,&
            source=angles_att, error=error)
    end if
    if (present(min_density_m)) then
       pao_env%min_d_m_valid=.true.
       ! to do
       CPAssert(.false.,cp_warning_level,routineP,error,failure)
       pao_env%min_density_m => min_density_m
    end if
    if (present(min_d_valid)) then
       pao_env%min_d_m_valid=min_d_valid
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine pao_env_set
!***************************************************************************

end module pao_env_methods

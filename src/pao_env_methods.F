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
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     02.2002 created [fawzi]
!!
!!   SOURCE
!*****************************************************************************
module pao_env_methods
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_assert, &
       cp_error_message, cp_assertion_failed, cp_error_init,&
       cp_error_get_print_level, cp_error_get_stop_level, cp_error_get_level,&
       cp_error_dealloc_ref, cp_error_get_logger, cp_internal_error
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types, only: pao_glob_angles_type, pao_env_type, cp_get, pao_env_valid
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use cp_block_matrix, only: cp_block_matrix_type, cp_block_matrix_p_type,&
       cp_sp_create, cp_get_matrix, cp_retain, cp_release
  use pao_glob_angles_methods, only: pao_g_ang_init, cp_dealloc, cp_set
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type
  use pao_bas_projection_types, only: pao_bas_projection_type,cp_proj_get
  use pao_qs_env_methods, only: cp_get_b_projection
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind, &
       get_atomic_kind_set
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type, &
       cp_bmstruct_create, cp_bmstruct_release
  use cp_para_types, only: cp_para_env_type
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
!!
!!   INPUTS
!!     - pao_env: the environement that has to be initialized
!!     - qs_env: the quickstep environement that contains this pao_env
!!     - global_env: the global environement
!!     - wanted_subspace_pos_precision: the wanted precision to find 
!!       the subspace
!!     - wanted_in_subspace_precision: the wanted precision in the subspace
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!     05.2002 updated doc [fawzi]
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
    type(pao_bas_projection_type), pointer :: bproj
    integer, dimension(:), allocatable :: min_bas_sizes, full_bas_sizes
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    integer, dimension(:), pointer :: atom_list
    character(len=60) :: name
    failure=.false.
    call cp_error_init(iError, template_error=error)

    nullify(pao_env%angles_begin_at, pao_env%atomic_ortho_valid,&
         pao_env%atomic_ortho, pao_env%min_density_m,&
         pao_env%min_bas_m_struct, pao_env%full_bas_m_struct,&
         pao_env%full_min_bas_m_struct, pao_env%min_full_bas_m_struct,&
         pao_env%full_bas_diag_m_struct)
    pao_env%min_d_m_valid=.false.

    ! check projection
    call get_qs_env(qs_env, atomic_kind_set= atomic_kind_set)
    n_atoms=size(qs_env%particle_set)
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
                min_bas_sizes(atom_list(j))=size_min_bas
                full_bas_sizes(atom_list(j))=size_full_bas
             end do
          end if
       end do
       
       ! min_bas_m_struct / full_bas_m_struct
       if (.not.failure) then
          call cp_bmstruct_create(pao_env%min_bas_m_struct,&
               symmetric=.true.,&
               para_env=global_env%para_env,&
               b_row_sizes=min_bas_sizes,&
               b_col_sizes=min_bas_sizes,&
               error=error)
          call cp_bmstruct_create(pao_env%full_bas_m_struct,&
               symmetric=.false.,&
               para_env=global_env%para_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
          call cp_bmstruct_create(pao_env%full_min_bas_m_struct,&
               symmetric=.false.,&
               para_env=global_env%para_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=min_bas_sizes,&
               error=error)
          call cp_bmstruct_create(pao_env%min_full_bas_m_struct,&
               symmetric=.false.,&
               para_env=global_env%para_env,&
               b_row_sizes=min_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
          call cp_bmstruct_create(pao_env%full_bas_diag_m_struct,&
               symmetric=.false.,&
               para_env=global_env%para_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
          pao_env%full_bas_diag_m_struct%sparsity=.false.
          do i=1,size(pao_env%full_bas_diag_m_struct%sparsity,1)
             pao_env%full_bas_diag_m_struct%sparsity(i,i)=.true.
          end do
       end if
       ! angles_begin_at
       if (.not.failure) then
          allocate(pao_env%angles_begin_at(n_atoms+1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          if (n_atoms>0) pao_env%angles_begin_at(1)=0
          do i=1,n_atoms
             pao_env%angles_begin_at(i+1)=pao_env%angles_begin_at(i)+&
                  min_bas_sizes(i)*(full_bas_sizes(i)-min_bas_sizes(i))
          end do
       end if

       deallocate(min_bas_sizes, full_bas_sizes, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if

    ! cache_max vals
    call get_atomic_kind_set(atomic_kind_set,&
         pao_max_min_bas=pao_env%max_min_bas,&
         pao_max_full_bas=pao_env%max_full_bas,&
         pao_max_excl_bas=pao_env%max_excl_bas,&
         pao_tot_min_bas=pao_env%tot_min_bas,&
         pao_tot_full_bas=pao_env%tot_full_bas,&
         pao_tot_excl_bas=pao_env%tot_excl_bas,&
         pao_max_l_angles=pao_env%max_l_angles)
        
    ! others
    n_atoms=size(qs_env%particle_set)
    allocate(pao_env%min_density_m(1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       do i=1,size(pao_env%min_density_m)
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
          call cp_sp_create(pao_env%min_density_m(i)%matrix,&
               matrix_struct=pao_env%min_bas_m_struct,&
               error=error)
       end do
    end if
    ! angles_att
    if (.not.failure) then
       call pao_g_ang_init(pao_env%angles_att, global_env=global_env,&
            qs_env=qs_env, error=error)
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
    end if
    call cp_error_dealloc_ref(iError,error)
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
!!     - pao_env: the environement to be deallocated
!!     - qs_env: the quickstep environement that contains this pao_env
!!     - global_env: the global environement
!!     - error: variable to control error logging, stopping,... 
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
       if (associated(pao_env%full_min_bas_m_struct)) then
          call cp_bmstruct_release(pao_env%full_min_bas_m_struct, error=error)
          nullify(pao_env%full_min_bas_m_struct)
       end if
       if (associated(pao_env%min_full_bas_m_struct)) then
          call cp_bmstruct_release(pao_env%min_full_bas_m_struct, error=error)
          nullify(pao_env%min_full_bas_m_struct)
       end if
       if (associated(pao_env%full_bas_diag_m_struct)) then
          call cp_bmstruct_release(pao_env%full_bas_diag_m_struct, error=error)
          nullify(pao_env%full_bas_diag_m_struct)
       end if
       if (associated(pao_env%atomic_ortho)) then
          call cp_release(pao_env%atomic_ortho,error=error)
       end if
       if (associated(pao_env%atomic_ortho_valid)) then
          deallocate(pao_env%atomic_ortho_valid, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
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
!!     - pao_env: the environement to be tested
!!     - qs_env: the quickstep environement that contains this pao_env
!!     - global_env: the global environement
!!     - test_pos_att: if it is false does not test the validity of the
!!       actual angles and unitary transformation. defaults to true
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
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
    failure=failure.and..not.pao_env_valid(pao_env,error=error)
    ! to do: more tests using qs_env
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
!!       Type(cp_block_matrix_p_type), Optional, Intent (IN), Dimension(:),&
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
!!     - pao_env: the environement that has to be initialized
!!     - qs_env: the quickstep environement where the pao_env lives
!!     - global_env: the global environement
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
    type(cp_block_matrix_p_type),optional, intent(in), dimension(:),&
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
    call cp_error_dealloc_ref(iError,error=error)
  end subroutine pao_env_set
!***************************************************************************

end module pao_env_methods

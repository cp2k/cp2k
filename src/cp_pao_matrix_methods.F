!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_matrix_methods [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_methods
!!
!!   FUNCTION
!!     methods to access and modify the pao matrices
!!
!!   NOTES
!!     type defined in cp_pao_types
!!     to do : renaming and public underlying functions
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
module cp_pao_matrix_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_pao_types
  use cp_pao_utils
  use matrix_types
  use cp_pao_plain_rot, only: cp_left_compose_p_rot
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_matrix_methods'

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, cp_create, cp_dealloc

! core functions
  public :: cp_validate, cp_valid

! special procedures
  !public :: cp_set_to_identity

! special functions
  public :: cp_get_pao_env, cp_get_matrix

! underlying functions
  public :: cp_pao_matrix_init, cp_pao_matrix_dealloc_ref,&
       cp_pao_matrix_create, cp_pao_matrix_dealloc,&
       cp_pao_matrix_validate, cp_pao_matrix_valid,&
       cp_pao_matrix_set, cp_pao_matrix_get,&
       cp_pao_matrix_get_pao_env,cp_pao_matrix_set_to_unity,&
       cp_pao_matrix_get_matrix

! ============ interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure cp_pao_matrix_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure cp_pao_matrix_dealloc_ref
  end interface

!! Allocates and initialize a new object returning a pointer to it in 
!! the first argument
  interface cp_create
     module procedure cp_pao_matrix_create
  end interface

!! Deallocates the memory used by the first argument, and the object
!! itself
  interface cp_dealloc
     module procedure cp_pao_matrix_dealloc
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure cp_pao_matrix_validate
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure cp_pao_matrix_valid
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure cp_pao_matrix_set
  end interface

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure cp_pao_matrix_get
  end interface

! +++ end of the core methods +++

!! returns the pao environement the first object lives in
  interface cp_get_pao_env
     module procedure cp_pao_matrix_get_pao_env
  end interface

!! sets the matrix to a unit matrix (and thus also guarantee that 
!! the diagonal blocks are present)
  interface cp_set_to_identity
     module procedure cp_pao_matrix_set_to_unity
  end interface

!! gets the underlying real matrix
  interface cp_get_matrix
     module procedure cp_pao_matrix_get_matrix
  end interface

!!***
!****************************************************************************

contains

! ====== core methods ======

!!****f* cp_pao_matrix_methods/cp_pao_matrix_init [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_init
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_init(pao_matrix, pao_env, matrix, source,&
!!         min_bas_rows, min_bas_cols, should_dealloc_matrix, error)
!!       Type(cp_pao_matrix_type), Intent (OUT):: pao_matrix
!!       Type(cp_pao_env_type), Target:: pao_env
!!       Type(real_matrix_type), Optional, Target:: matrix
!!       Type(cp_pao_matrix_type), Intent (IN), Optional:: source
!!       Logical, Intent (IN), Optional:: min_bas_cols, min_bas_rows,&
!!         should_dealloc_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_init
!!
!!   FUNCTION
!!     initializes a pao matrix. If no matrix is given allocates a block
!!     diagonal matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     For the others see cp_pao_matrix_set (that is called
!!     by this routine).
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *********************************************************************
  subroutine cp_pao_matrix_init(pao_matrix, pao_env,&
       matrix, source, min_bas_rows, min_bas_cols,&
       should_dealloc_matrix, error)
    type(cp_pao_matrix_type), intent(out) :: pao_matrix
    type(cp_pao_env_type), target :: pao_env
    type(real_matrix_type), optional, target :: matrix
    type(cp_pao_matrix_type), intent(in), optional :: source
    logical, intent(in), optional :: min_bas_cols, min_bas_rows, &
         should_dealloc_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_init'
    failure=.false.

    pao_matrix%initialized=.true.
    pao_matrix%should_dealloc_matrix=.true.
    nullify(pao_matrix%matrix)
! nullify(pao_matrix%pao_env)
    pao_matrix%pao_env => pao_env
    pao_matrix%min_bas_rows=.true.
    pao_matrix%min_bas_rows=.true.
    call cp_pao_matrix_set(pao_matrix,pao_env,&
         matrix, source, min_bas_rows, min_bas_cols,&
         should_dealloc_matrix, error)
  end subroutine cp_pao_matrix_init
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_dealloc_ref(pao_matrix, error)
!!       Type(cp_pao_matrix_type), Intent (INOUT):: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_dealloc_ref
!!
!!   FUNCTION
!!     Dealloc the memory that the matrix has allocated
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that should dealloc its memory
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
!!   SOURCE
!!*** **********************************************************************
  subroutine cp_pao_matrix_dealloc_ref(pao_matrix, error)
    type(cp_pao_matrix_type), intent(inout) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_dealloc_ref'
    failure=.false.

    CPPrecondition(pao_matrix%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       nullify(pao_matrix%pao_env)
       if (pao_matrix%should_dealloc_matrix) then
          if (associated(pao_matrix%matrix)) then
             call deallocate_matrix(pao_matrix%matrix)
          end if
       else
          nullify(pao_matrix%matrix)
       end if
    end if failureIf
  end subroutine cp_pao_matrix_dealloc_ref
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_create [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_create
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_create(pao_matrix, pao_env, matrix, source,&
!!         min_bas_rows, min_bas_cols, should_dealloc_matrix, error)
!!       Type(cp_pao_matrix_type), Pointer:: pao_matrix
!!       Type(cp_pao_env_type), Target:: pao_env
!!       Type(real_matrix_type), Optional, Target:: matrix
!!       Type(cp_pao_matrix_type), Intent (IN), Optional:: source
!!       Logical, Intent (IN), Optional:: min_bas_cols, min_bas_rows,&
!!         should_dealloc_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_create
!!
!!   FUNCTION
!!     allocates the matrix and initializes it
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: will point to the newly created matrix
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     for the others arguments see cp_init_block_pao_matrix
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_matrix_create(pao_matrix, pao_env,&
       matrix, source, min_bas_rows, min_bas_cols,&
       should_dealloc_matrix, error)
    type(cp_pao_matrix_type), pointer :: pao_matrix
    type(cp_pao_env_type), target :: pao_env
    type(real_matrix_type), optional, target :: matrix
    type(cp_pao_matrix_type), intent(in), optional :: source
    logical, intent(in), optional :: min_bas_cols, min_bas_rows, &
         should_dealloc_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_create'
    integer :: stat
    failure=.false.

    allocate(pao_matrix, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_pao_matrix_init(pao_matrix, pao_env,&
            matrix, source, min_bas_rows, min_bas_cols,&
            should_dealloc_matrix, error)
    end if
  end subroutine cp_pao_matrix_create
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_dealloc [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_dealloc(pao_matrix, error)
!!       Type(cp_pao_matrix_type), Pointer:: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_dealloc
!!
!!   FUNCTION
!!     deallocates the memory allocated by the matrix ans the matrix
!!     object itself.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix to be deallocated
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
  subroutine cp_pao_matrix_dealloc(pao_matrix, error)
    type(cp_pao_matrix_type), pointer :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_dealloc'
    integer :: stat
    failure=.false.

    if (associated(pao_matrix)) then
       call cp_pao_matrix_dealloc_ref(pao_matrix,error=error)
       deallocate(pao_matrix,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
  end subroutine cp_pao_matrix_dealloc
!***************************************************************************


!!****f* cp_pao_matrix_methods/cp_pao_matrix_validate [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_validate
!!
!!   SYNOPSIS
!!     Function cp_pao_matrix_validate(pao_matrix, error) Result(res)
!!       Logical:: res
!!       Type(cp_pao_matrix_type), Intent (IN):: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_matrix_validate
!!
!!   FUNCTION
!!     controls if the matrix  has the right size, returns true
!!     if it is valid, false if it is not.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix to be tested
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
  function cp_pao_matrix_validate(pao_matrix, error) result(res)
    logical ::res
    type(cp_pao_matrix_type), intent(in) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_validate'
    integer, dimension(:), pointer :: cols_size, rows_size
    failure=.false.
    nullify(cols_size,rows_size)

    CPPrecondition(pao_matrix%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       CPAssert(associated(pao_matrix%matrix),cp_warning_level,routineP,error,failure)
       CPAssert(associated(pao_matrix%pao_env), cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          CPAssert(associated(pao_matrix%pao_env%nr_full_bas),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_matrix%pao_env%nr_min_bas),cp_warning_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_pao_matrix_get(pao_matrix, cols_size=cols_size,&
               rows_size=rows_size)
          call cp_assert(cp_consistent_matrix(matrix=pao_matrix%matrix,&
               rows_size=rows_size, cols_size=cols_size,error=error),&
               cp_warning_level, cp_assertion_failed, routineP,&
               "size in inconsistent with given basis! in "//&
               CPSourceFileRef,&
               error, failure)
       end if
    end if failureIf
    res = .not.failure
  end function cp_pao_matrix_validate
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_set [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_set
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_set(pao_matrix, pao_env, matrix, source,&
!!         min_bas_rows, min_bas_cols, should_dealloc_matrix, error)
!!       Type(cp_pao_matrix_type), Intent (INOUT):: pao_matrix
!!       Type(cp_pao_env_type), Optional, Target:: pao_env
!!       Type(real_matrix_type), Optional, Target:: matrix
!!       Type(cp_pao_matrix_type), Intent (IN), Optional:: source
!!       Logical, Intent (IN), Optional:: min_bas_cols, min_bas_rows,&
!!         should_dealloc_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_set
!!
!!   FUNCTION
!!     sets the values of the various elements of the stucture
!!
!!   NOTES
!!     you should rarely nedd this (initialize in init)
!!
!!   INPUTS
!!     pao_matrix: the transformation you want to change
!!     pao_env: the pao environement this transformation refers to.
!!           it is SHARED, and it should be already valid.
!!     matrix: normally COPIES the matrix into this object, but if at the
!!             same time you also set should_dealloc_matrix to false
!!             (i.e. the matrix is shared) it is SHARED.
!!     source: a transformation that is copied into this one(other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     min_bas_rows: if it is true (default) the rows are in the 
!!             minimal basis
!!     min_bas_cols: if it is true (default) the cols are in the
!!             minimal basis
!!     should_dealloc_matrix: if true (default) the matrix will be 
!!             deallocated in cp_dealloc_ref
!!     error: variable to control error logging, stopping,... 
!!           see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_matrix_set(pao_matrix, pao_env,&
       matrix, source, min_bas_rows, min_bas_cols,&
       should_dealloc_matrix, error)
    type(cp_pao_matrix_type), intent(inout) :: pao_matrix
    type(cp_pao_env_type), optional, target :: pao_env
    type(real_matrix_type), optional, target :: matrix
    type(cp_pao_matrix_type), intent(in), optional :: source
    logical, intent(in), optional :: min_bas_cols, min_bas_rows,&
         should_dealloc_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, should_share_matrix
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_set'
    integer :: stat
    character (len=60) :: matrix_name
    type(real_matrix_type), pointer :: tmp_ptr
    failure=.false.
    nullify(tmp_ptr)

    CPPrecondition(pao_matrix%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       sourceIf: if (present(source)) then
          CPPrecondition(source%initialized, cp_failure_level,routineP,error,failure)
          if (present(should_dealloc_matrix)) then
             should_share_matrix=.not.should_share_matrix
          else
             should_share_matrix=.false.
          end if
          if (.not.failure) then
             pao_matrix%pao_env => source%pao_env
             pao_matrix%min_bas_rows=source%min_bas_rows
             pao_matrix%min_bas_cols=source%min_bas_cols
             if (associated(source%matrix)) then
                if (.not.present(should_dealloc_matrix)) &
                     should_share_matrix=.not.source%should_dealloc_matrix
                if (.not.should_share_matrix) then
                   if (.not. associated(pao_matrix%matrix)) then
                      call get_matrix_info(source%matrix,matrix_name=matrix_name)
                      call replicate_matrix(source%matrix, &
                           pao_matrix%matrix,&
                           matrix_name)
                   else
                      call copy_matrix(source%matrix, &
                           pao_matrix%matrix)
                   end if
                else ! the matrix is shared
                   pao_matrix%matrix => source%matrix
                end if
             else if (associated(pao_matrix%matrix))then
                if (pao_matrix%should_dealloc_matrix) then
                   deallocate(pao_matrix%matrix, stat=stat)
                   CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
                else
                   nullify(pao_matrix%matrix)
                end if
             end if
          end if
          pao_matrix%should_dealloc_matrix=source%should_dealloc_matrix
       end if sourceIf
       if (present(pao_env)) then
          pao_matrix%pao_env => pao_env
       end if
       matrixIf: if (present(matrix)) then
          if (should_share_matrix) then
             pao_matrix%matrix => matrix
          else
             tmp_ptr => matrix
             if (.not. associated(pao_matrix%matrix)) then
                call get_matrix_info(source%matrix,matrix_name=matrix_name)
                call replicate_matrix(tmp_ptr, &
                     pao_matrix%matrix,&
                     matrix_name)
             else
                call copy_matrix(tmp_ptr, &
                     pao_matrix%matrix)
             end if
          end if
       end if matrixIf
       if (present(min_bas_rows)) &
            pao_matrix%min_bas_rows=min_bas_rows
       if (present(min_bas_cols)) &
            pao_matrix%min_bas_cols=min_bas_cols
       if (present(should_dealloc_matrix)) &
            pao_matrix%should_dealloc_matrix=should_dealloc_matrix
       if (associated(pao_matrix%pao_env).and. &
            .not. associated(pao_matrix%matrix)) then
          call cp_init_block_pao_matrix(pao_matrix, error=error)
       end if
    end if failureIf
  end subroutine cp_pao_matrix_set
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_get [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_get
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_get(pao_matrix, pao_env, matrix,&
!!         min_bas_cols, min_bas_rows, cols_size, rows_size,&
!!         should_dealloc_matrix, error)
!!       Type(cp_pao_matrix_type), Intent (IN):: pao_matrix
!!       Type(cp_pao_env_type), Optional, Pointer:: pao_env
!!       Type(real_matrix_type), Optional, Pointer:: matrix
!!       Logical, Optional, Intent (OUT):: min_bas_rows, min_bas_cols,&
!!         should_dealloc_matrix
!!       Integer, Dimension(:), Pointer, Optional:: cols_size, rows_size
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_get
!!
!!   FUNCTION
!!     Get the attributes of the given matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the pao transformation you want info from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     For the other arguments see the attributes of cp_pao_matrix_type
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_matrix_get(pao_matrix,pao_env,matrix,min_bas_cols,&
       min_bas_rows, cols_size, rows_size, should_dealloc_matrix, error)
    type(cp_pao_matrix_type), intent(in) :: pao_matrix
    type(cp_pao_env_type), optional, pointer :: pao_env
    type(real_matrix_type), optional, pointer :: matrix
    logical, optional, intent(out) :: min_bas_rows, min_bas_cols, &
         should_dealloc_matrix
    integer, dimension(:), pointer, optional :: cols_size, rows_size
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_get'
    failure=.false.

    CPPrecondition(pao_matrix%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(pao_env)) pao_env => pao_matrix%pao_env
       if (present(matrix)) matrix => pao_matrix%matrix
       if (present(min_bas_cols)) min_bas_cols=pao_matrix%min_bas_cols
       if (present(min_bas_rows)) min_bas_rows=pao_matrix%min_bas_rows
       if (present(should_dealloc_matrix)) &
            should_dealloc_matrix=pao_matrix%should_dealloc_matrix
       if (present(cols_size)) then
          CPPrecondition(associated(pao_matrix%pao_env),cp_warning_level,routineP,error,failure)
          if (.not.failure) then
             if (pao_matrix%min_bas_cols) then
                cols_size => pao_matrix%pao_env%nr_min_bas
             else
                cols_size => pao_matrix%pao_env%nr_full_bas
             end if
          else
             nullify(cols_size)
          end if
       end if
       if (present(rows_size)) then
          CPPrecondition(associated(pao_matrix%pao_env),cp_warning_level,routineP,error,failure)
          if (.not.failure) then
             if (pao_matrix%min_bas_rows) then
                rows_size => pao_matrix%pao_env%nr_min_bas
             else
                rows_size => pao_matrix%pao_env%nr_full_bas
             end if
          else
             nullify(rows_size)
          end if
       end if
    end if failureIf
  end subroutine cp_pao_matrix_get
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_valid [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_valid
!!
!!   SYNOPSIS
!!     Function cp_pao_matrix_valid(pao_matrix, error) Result(res)
!!       Logical:: res
!!       Type(cp_pao_matrix_type), Intent (IN):: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_matrix_valid
!!
!!   FUNCTION
!!     performs minimal validadion (no dangling pointers)
!!
!!   NOTES
!!     can be called often
!!
!!   INPUTS
!!     pao_matrix: the object to validate
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
  function cp_pao_matrix_valid(pao_matrix, error) result(res)
    logical ::res
    type(cp_pao_matrix_type), intent(in) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_valid'
    failure=.false.

    CPPrecondition(pao_matrix%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (.not.(associated(pao_matrix%matrix).and.&
            associated(pao_matrix%pao_env)) ) then
          CPAssert(associated(pao_matrix%matrix),cp_warning_level,routineP,error,failure)
          CPAssert(associated(pao_matrix%pao_env),cp_warning_level,routineP,error,failure)
       end if
    end if
    res=.not.failure
  end function cp_pao_matrix_valid
!***************************************************************************

! ===== special methods =====

!!****f* cp_pao_matrix_methods/cp_pao_matrix_set_to_unity [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_set_to_unity
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_matrix_set_to_unity(pao_matrix, error)
!!       Type(cp_pao_matrix_type), Intent (INOUT):: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_matrix_set_to_unity
!!
!!   FUNCTION
!!     sets the matrix to the unit matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix to be initialized
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
  subroutine cp_pao_matrix_set_to_unity(pao_matrix, error)
    type(cp_pao_matrix_type), intent(inout) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_set_to_unity'
    type(cp_error_type) :: iError
    failure=.false.
    call cp_init(iError,template_error=error)

    CPPrecondition(pao_matrix%initialized,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       if (.not.associated(pao_matrix%matrix)) then
          call cp_init_block_pao_matrix(pao_matrix, error=iError)
          call cp_propagate_error(iError,routineP,&
               CPSourceFileRef,&
               error=error,failure=failure)
       end if
    end if
    if (.not. failure) then
       CPAssert(.false.,cp_failure_level,routineP,error,failure)
       ! to do
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_matrix_set_to_unity
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_get_pao_env [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_get_pao_env
!!
!!   SYNOPSIS
!!     Function cp_pao_matrix_get_pao_env(pao_matrix, error) Result(res)
!!       Type(cp_pao_env_type), Pointer:: res
!!       Type(cp_pao_matrix_type), Intent (IN):: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_matrix_get_pao_env
!!
!!   FUNCTION
!!     returns the global environement the matrix object lives in
!!
!!   NOTES
!!     this might go away (to optimize), so do not rely on it too heavily
!!
!!   INPUTS
!!     pao_matrix: the matrix object you get the info from
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
  function cp_pao_matrix_get_pao_env(pao_matrix, error) result(res)
    type(cp_pao_env_type), pointer ::res
    type(cp_pao_matrix_type), intent(in) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_get_pao_env'
    failure=.false.
    CPPrecondition(pao_matrix%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       res => pao_matrix%pao_env
    else failureIf
       nullify(res)
    end if failureIf
  end function cp_pao_matrix_get_pao_env
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_pao_matrix_get_matrix [1.0] *
!!
!!   NAME
!!     cp_pao_matrix_get_matrix
!!
!!   SYNOPSIS
!!     Function cp_pao_matrix_get_matrix(pao_matrix, error) Result(res)
!!       Type(real_matrix_type), Pointer:: res
!!       Type(cp_pao_matrix_type), Intent (IN):: pao_matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_matrix_get_matrix
!!
!!   FUNCTION
!!     returns the undelying matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_matrix: the matrix you get the info from
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
  function cp_pao_matrix_get_matrix(pao_matrix, error) result(res)
    type(real_matrix_type), pointer ::res
    type(cp_pao_matrix_type), intent(in) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_pao_matrix_get_matrix'
    failure=.false.
    CPPrecondition(pao_matrix%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       res => pao_matrix%matrix
    else
       nullify(res)
    end if
  end function cp_pao_matrix_get_matrix
!***************************************************************************


! ===== private utilities =====

!!****f* cp_pao_matrix_methods/cp_init_block_pao_matrix [1.0] *
!!
!!   NAME
!!     cp_init_block_pao_matrix
!!
!!   FUNCTION
!!     if the matrix is not already right
!!     allocate a blocks sizes of the matrix using the pao_env info
!!     but no blocks are allocated (use cp_set_to_identity to have a 
!!     block diagonal matrix).
!!
!!   NOTES
!!     pao_env should be valid.
!!     private method
!!
!!   INPUTS
!!     pao_matrix: the trasformation that should be initialized
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
  subroutine cp_init_block_pao_matrix(pao_matrix, error)
    type(cp_pao_matrix_type), intent(inout) :: pao_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_init_block_pao_matrix'
    integer :: n_atom, i, stat
    integer, dimension(:), allocatable :: firstRow, lastRow, firstCol, lastCol
    integer, dimension(:), pointer :: cols_size, rows_size
    failure=.false.
    n_atom=-1

    CPPrecondition(pao_matrix%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(pao_matrix%pao_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_assert(cp_valid(pao_matrix%pao_env, &
            error=error), cp_failure_level, &
            cp_precondition_failed, routineP,& 
            'pao_env must be valid to initialize a pao matrix in '//&
            CPSourceFileRef,&
            error,failure) !remove this test?
    end if
    if (.not.failure) then
       call cp_get(pao_matrix, cols_size=cols_size, rows_size=rows_size)
       CPPrecondition(associated(cols_size),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(rows_size),cp_failure_level,routineP,error,failure)
    end if
    failureIf: if (.not. failure) then
       if (associated(pao_matrix%matrix)) then
          if (.not.cp_consistent_matrix(pao_matrix%matrix,&
               cols_size=cols_size,&
               rows_size=rows_size, error=error)) then
             CPErrorMessage(cp_warning_level,routineP,"the matrix was inconsistent",error)
             call deallocate_matrix(pao_matrix%matrix)
          end if
       end if
       if (.not.associated(pao_matrix%matrix)) then
          n_atom=cp_get_n_atom(pao_matrix%pao_env)
          CPPrecondition(n_atom==size(cols_size),cp_failure_level,routineP,error,failure)
          CPPrecondition(n_atom==size(rows_size),cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             allocate(firstRow(n_atom), lastRow(n_atom),&
                  firstCol(n_atom), lastCol(n_atom), stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
          if (.not.failure) then
             call cp_first_last_from_size(cols_size,first_index=firstCol,&
                  last_index=lastCol,error=error)
             call cp_first_last_from_size(rows_size,first_index=firstRow,&
                  last_index=lastRow,error=error)
             call allocate_matrix(matrix=pao_matrix%matrix, &
                  nblock_row=n_atom, nblock_col=n_atom, &
                  nrow=sum(rows_size), ncol=sum(cols_size), &
                  first_row=firstRow, last_row=lastRow, &
                  first_col=firstCol, last_col=lastCol, &
                  matrix_name="pao matrix", &
                  matrix_symmetry="block diagonal")
             if (debug_this_module .and. cp_debug) then
                call cp_assert(&
                     cp_consistent_matrix(pao_matrix%matrix,&
                     cols_size=cols_size, rows_size=rows_size,error=error),&
                     cp_failure_level,cp_assertion_failed,&
                     routineP,"internal error: after its init the pao "//&
                     "transformation is not consistent in "//&
                     CPSourceFileRef,&
                     error, failure)
             end if
             deallocate(firstRow,lastRow, firstCol, lastCol, stat=stat) 
! just to be sure...
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
    end if failureIf
  end subroutine cp_init_block_pao_matrix
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_first_last_from_size [1.0] *
!!
!!   NAME
!!     cp_first_last_from_size
!!
!!   FUNCTION
!!     generates the first index and last index from the sizes of the blocks
!!
!!   NOTES
!!     private method
!!
!!   INPUTS
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
  subroutine cp_first_last_from_size(block_size, first_index, last_index, error)
    type(cp_error_type), optional, intent(inout) :: error
    integer, dimension(:), intent(out) :: first_index, last_index
    integer, dimension(:), intent(in) :: block_size
    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_first_last_from_size'
    integer :: n_blocks, i
    failure=.false.
    n_blocks=size(block_size)

    CPPrecondition(size(first_index)==n_blocks,cp_failure_level,routineP,error,failure)
    CPPrecondition(size(last_index)==n_blocks,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (n_blocks>0) then
          first_index(1)=1
       end if
       do i=1,n_blocks-1
          last_index(i)=first_index(i)+block_size(i)-1
          first_index(i+1)=last_index(i)+1
       end do
       if (n_blocks>0) then
          last_index(n_blocks)=first_index(n_blocks)+block_size(n_blocks)
       end if
    end if
  end subroutine cp_first_last_from_size
!***************************************************************************

!!****f* cp_pao_matrix_methods/cp_consistent_matrix [1.0] *
!!
!!   NAME
!!     cp_consistent_matrix
!!
!!   FUNCTION
!!     Returns true if the matrix size is consistent
!!     with the given sizes
!!
!!   NOTES
!!     private method
!!
!!   INPUTS
!!     matrix: the matrix to check
!!     cols_size: the size of the columns of the blocks
!!     rows_size the size of the rows of the blocks
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
  function cp_consistent_matrix(matrix,cols_size,rows_size, error) result(res)
    logical ::res
    type(real_matrix_type), intent(in), target :: matrix
    integer, dimension(:), intent(in) :: cols_size, rows_size
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_matrix_methods:cp_consistent_matrix'
    integer :: n_atom, i, stat
    integer, dimension(size(cols_size)) :: firstRow, lastRow,&
         firstCol,lastCol
    integer, dimension(:), pointer :: nr_full_bas, first_col, last_col, &
         first_row, last_row
    type(real_matrix_type), pointer :: tstPtr
    failure=.false.
    nullify(nr_full_bas, first_col, last_col, first_row, last_row)
    tstPtr => matrix

    CPPrecondition(associated(tstPtr),cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       n_atom=size(cols_size)
       CPPrecondition(n_atom==size(rows_size),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          call cp_first_last_from_size(cols_size,first_index=firstCol,&
               last_index=lastCol)
          call cp_first_last_from_size(rows_size,first_index=firstRow,&
               last_index=lastRow)
          call get_matrix_info(tstPtr, &
               first_col=first_col , last_col=last_col, &
               first_row=first_row , last_row=last_row)
          CPAssert(associated(first_row).and.associated(last_row),cp_warning_level,routineP,error,failure)
          CPAssert(associated(first_col).and.associated(last_row),cp_warning_level,routineP,error,failure)
          if (.not.failure) then
             CPAssert(size(first_row)==n_atom,cp_warning_level,routineP,error,failure)
             CPAssert(size(last_row)==n_atom,cp_warning_level,routineP,error,failure)
             CPAssert(size(first_col)==n_atom,cp_warning_level,routineP,error,failure)
             CPAssert(size(last_col)==n_atom,cp_warning_level,routineP,error,failure)
             if (.not.failure) then
                CPAssert(all(first_row==firstRow),cp_warning_level,routineP,error,failure)
                CPAssert(all(first_col==firstCol),cp_warning_level,routineP,error,failure)
                CPAssert(all(last_row==lastRow),cp_warning_level,routineP,error,failure)
                CPAssert(all(last_col==lastCol),cp_warning_level,routineP,error,failure)
             end if
          end if
          nullify(first_col, last_col, first_row, last_row)
       end if
    end if failureIf
    res= .not. failure
  end function cp_consistent_matrix
!***************************************************************************
  
end module cp_pao_matrix_methods

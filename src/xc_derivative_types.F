!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/xc_derivative_types *
!!
!!   NAME
!!     xc_derivative_types
!!
!!   FUNCTION
!!     Provides types for the management of the xc-functionals and
!!     their derivatives.
!!
!!   NOTES
!!     
!!     
!!****
module xc_derivative_types

  use kinds,              only: dp
  use cp_error_handling,  only: cp_error_type, cp_a_l, cp_assert, cp_error_get_logger
  use cp_log_handling,    only: cp_logger_type, cp_failure_level, cp_warning_level, &
                                cp_to_string, cp_logger_get_default_unit_nr
  use xc_derivative_desc, only: MAX_DERIVATIVE_DESC_LENGTH, &
                                standardize_dervitative_desc

  implicit none

  private

  integer, save :: derivative_id_nr = 1
  integer, save :: derivative_set_id_nr  = 1

  character(len=*), parameter :: moduleN='xc_derivative_types'

  !!****t* cp2k/xc_derivative_type  *
  !!
  !!   NAME
  !!     xc_derivative_type
  !!
  !!   FUNCTION
  !!     A node of a linked list contains the actual derivative together
  !!     with a string describing it.
  !!
  !!****
  type xc_derivative_type
     integer                                   :: ref_count, id_nr
     character(len=MAX_DERIVATIVE_DESC_LENGTH) :: description
     real(KIND = dp), dimension(:), pointer    :: data
     type(xc_derivative_type), pointer         :: next
  end type xc_derivative_type

  !!****t* cp2k/xc_derivative_set_type  *
  !!
  !!   NAME
  !!     xc_derivative_set_type
  !!
  !!   FUNCTION
  !!     A derivative set contains the different derivatives of a xc-functional
  !!     in form of a linked list
  !!
  !!   NOTE
  !!     To all USERS of this type: never reference data directly. data is an
  !!     internal variable, which due to fortrans limitations cannot be made
  !!     private
  !!
  !!****
  type xc_derivative_set_type
     integer                           :: ref_count, id_nr, npoints
     type(xc_derivative_type), pointer :: first
     real(KIND = dp), dimension(:,:), pointer :: data
  end type xc_derivative_set_type

  !!****t* cp2k/xc_rho_set_type  *
  !!
  !!   NAME
  !!     xc_rho_set_type
  !!
  !!   FUNCTION
  !!     A commodity type to hold the various needed (spin) densities
  !!     and gradients of densities
  !!
  !!****
  type xc_rho_set_type
     integer                           :: ref_count, id_nr, npoints
     ! for spin restricted systems
     real(KIND = dp), dimension(:,:,:), pointer :: rho
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drho
     real(KIND = dp), dimension(:,:,:), pointer :: norm_drho
     ! for UNrestricted systems
     real(KIND = dp), dimension(:,:,:), pointer :: rhoa, rhob
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drhoa, drhob
     real(KIND = dp), dimension(:,:,:), pointer :: norm_drhoa, norm_drhob
     real(KIND = dp), dimension(:,:,:), pointer :: norm_drhoa_drhob
  end type xc_rho_set_type


  public :: xc_derivative_type, &
            xc_derivative_set_type, &
            xc_derivative_set_create, &
            xc_derivative_set_release, &
            xc_dset_get_derivative, &
            xc_rho_set_type

contains

  ! *****************************************************************************

  !!****f* xc_derivative_types/xc_dset_get_derivative [1.0] *
  !!
  !!   NAME
  !!     xc_dset_get_derivative
  !!
  !!   FUNCTION
  !!     returns the requested xc_derivative
  !!
  !!   NOTES
  !!     -
  !!     
  !!   INPUTS
  !!     - derivative_set: the set where to search for the derivative
  !!     - description   : the description of the derivative you want to have
  !!     - allocate_deriv: if the derivative should be allocated when not present
  !!                       Defaults to false.
  !!     - error         : variable to control error logging, stopping,... 
  !!                       see module cp_error_handling 
  !!
  !!*** **********************************************************************
  function xc_dset_get_derivative(derivative_set, description, allocate_deriv, error) &
       result(res)

    ! arguments
    type(xc_derivative_set_type), pointer        :: derivative_set
    character(len=*), intent(in)                 :: description
    logical, intent(in), optional                :: allocate_deriv
    type(cp_error_type), optional, intent(inout) :: error
    type(xc_derivative_type), pointer            :: res

    ! locals
    logical                                   :: failure, my_allocate_deriv
    character(len=MAX_DERIVATIVE_DESC_LENGTH) :: std_deriv_desc
    type(xc_derivative_type), pointer         :: pos

    ! parameters
    character(len=*), parameter :: routineN = 'xc_dset_get_derivative',&
                                   routineP = moduleN//':'//routineN

    failure=.false.
    CPPrecondition(associated(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative_set%ref_count>0,cp_failure_level,routineP,error,failure)
    my_allocate_deriv=.false.
    if (present(allocate_deriv)) my_allocate_deriv=allocate_deriv
    nullify(res)
    std_deriv_desc = standardize_dervitative_desc(description)
    if (.not. failure) then
       pos => derivative_set%first
       do while(associated(pos))
          if (pos%description == std_deriv_desc) then
             res => pos
             exit
          end if
          pos => pos%next
       end do
       if (.not.associated(res).and.my_allocate_deriv) then
          call xc_derivative_create(res, std_deriv_desc, & 
                                    derivative_set%npoints, error=error)
          res%next => derivative_set%first
          derivative_set%first => res
       end if
    end if
  end function xc_dset_get_derivative
  !***************************************************************************

  !!****f* xc_derivative_types/xc_derivative_set_create [1.0] *
  !!
  !!   NAME
  !!     xc_derivative_set_create
  !!
  !!   FUNCTION
  !!     creates a derivative set object
  !!
  !!   INPUTS
  !!     - derivative_set : the set where to search for the derivative
  !!     - npoints : the number of points in the derivatives
  !!     - error   : variable to control error logging, stopping,... 
  !!                 see module cp_error_handling 
  !!
  !!*** **********************************************************************
  subroutine xc_derivative_set_create(derivative_set, npoints, error)

    ! arguments
    type(xc_derivative_set_type), pointer        :: derivative_set
    integer, intent(in)                          :: npoints
    type(cp_error_type), optional, intent(inout) :: error

    ! locals
    logical :: failure
    integer :: stat

    ! parameters
    character(len=*), parameter :: routineN = 'xc_derivative_set_create',&
                                   routineP = moduleN//':'//routineN

    failure=.false.
    
    CPPrecondition(.not.associated(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(npoints>=0,cp_failure_level,routineP,error,failure)

    allocate(derivative_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    if (.not. failure) then
       nullify(derivative_set%first)
       nullify(derivative_set%data)
       derivative_set%ref_count  = 1
       derivative_set%id_nr      = derivative_set_id_nr
       derivative_set_id_nr      = derivative_set_id_nr + 1
       derivative_set%npoints    = npoints
    end if

  end subroutine xc_derivative_set_create
  !***************************************************************************

  !!****f* xc_derivative_types/xc_derivative_set_release [1.0] *
  !!
  !!   NAME
  !!     xc_derivative_set_release
  !!
  !!   FUNCTION
  !!     releases a derivative set
  !!
  !!   INPUTS
  !!     - set     : the set to release
  !!     - error   : variable to control error logging, stopping,... 
  !!                 see module cp_error_handling 
  !!
  !!*** **********************************************************************
  subroutine xc_derivative_set_release(derivative_set, error)

    ! arguments
    type(xc_derivative_set_type), pointer        :: derivative_set
    type(cp_error_type), optional, intent(inout) :: error

    ! locals
    logical :: failure
    integer :: stat
    type(xc_derivative_type), pointer :: pos, tmp

    ! parameters
    character(len=*), parameter :: routineN = 'xc_derivative_set_release',&
                                   routineP = moduleN//':'//routineN

    failure=.false.
    
    CPPrecondition(associated(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative_set%ref_count>0,cp_failure_level,routineP,error,failure)

    if (.not. failure) then
       derivative_set%ref_count = derivative_set%ref_count - 1
       if (derivative_set%ref_count == 0) then
          if (associated(derivative_set%data)) then
             deallocate(derivative_set%data, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          end if
          pos => derivative_set%first
          do while (associated(pos))
             tmp => pos
             pos => pos%next
             call xc_derivative_release(tmp, error=error)
          end do
          deallocate(derivative_set, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       end if
    end if
    NULLIFY(derivative_set)

  end subroutine xc_derivative_set_release
  !***************************************************************************

  !!****f* xc_derivative_types/xc_derivative_create [1.0] *
  !!
  !!   NAME
  !!     xc_derivative_create
  !!
  !!   FUNCTION
  !!     allocates and initializes a derivative type
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - derivative  : the object to create
  !!     - description : the derivative description
  !!     - npoints     : size of the data array
  !!     - error       : variable to control error logging, stopping,... 
  !!                     see module cp_error_handling 
  !!
  !!*** **********************************************************************
  subroutine xc_derivative_create(derivative, description, npoints, error)

    ! arguments
    type(xc_derivative_type), pointer              :: derivative
    character(len=*), intent(in)                   :: description
    integer, intent(in)                            :: npoints
    type(cp_error_type), optional, intent(inout)   :: error

    ! locals
    logical                                   :: failure
    integer                                   :: stat
    character(len=MAX_DERIVATIVE_DESC_LENGTH) :: my_description

    ! parameters
    character(len=*), parameter :: routineN = 'xc_derivative_create',&
         routineP = moduleN//':'//routineN

    failure=.false.

    allocate(derivative, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    if (.not. failure) then
       nullify(derivative%data)
       nullify(derivative%next)
       derivative%ref_count = 1
       derivative%id_nr     = derivative_id_nr
       derivative_id_nr     = derivative_id_nr + 1
       my_description = standardize_dervitative_desc(description)
       derivative%description = my_description
       allocate(derivative%data(npoints), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if

  end subroutine xc_derivative_create
  !***************************************************************************

  !!****f* xc_derivative_types/xc_derivative_release [1.0] *
  !!
  !!   NAME
  !!     xc_derivative_release
  !!
  !!   FUNCTION
  !!     allocates and initializes a derivative type
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - derivative  : the object to create
  !!     - error       : variable to control error logging, stopping,... 
  !!                     see module cp_error_handling 
  !!
  !!*** **********************************************************************
  subroutine xc_derivative_release(derivative, error)

    ! arguments
    type(xc_derivative_type), pointer              :: derivative
    type(cp_error_type), optional, intent(inout)   :: error

    ! locals
    logical                                   :: failure
    integer                                   :: stat
    character(len=MAX_DERIVATIVE_DESC_LENGTH) :: my_description

    ! parameters
    character(len=*), parameter :: routineN = 'xc_derivative_release',&
                                   routineP = moduleN//':'//routineN

    failure=.false.

    CPPrecondition(associated(derivative),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative%ref_count>=1,cp_failure_level,routineP,error,failure)

    if (.not.failure) then
       derivative%ref_count = derivative%ref_count - 1
       if (derivative%ref_count == 0) then
          deallocate(derivative%data, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          deallocate(derivative, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       end if
    end if
    NULLIFY(derivative)

  end subroutine xc_derivative_release
  !***************************************************************************

  !!****f* xc_derivative_types/xc_rho_set_create [1.0] *
  !!
  !!   NAME
  !!     xc_rho_set_create
  !!
  !!   FUNCTION
  !!     creates an rho_set object
  !!
  !!   INPUTS
  !!     - rho_set : the rho set to create
  !!     - error   : variable to control error logging, stopping,... 
  !!                 see module cp_error_handling 
  !!
  !!*** **********************************************************************
  subroutine xc_rho_set_create(rho_set, error)

    ! arguments
    type(xc_derivative_set_type), pointer        :: rho_set
    type(cp_error_type), optional, intent(inout) :: error

    ! locals
    logical :: failure
    integer :: stat

    ! parameters
    character(len=*), parameter :: routineN = 'xc_rho_set_create',&
                                   routineP = moduleN//':'//routineN

    failure=.false.
    
    CPPrecondition(.not.associated(rho_set),cp_failure_level,routineP,error,failure)

    if (.not. failure) then

    end if

  end subroutine xc_rho_set_create
  !***************************************************************************

end module xc_derivative_types


!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_derivatives [1.0] *
!!
!!   NAME
!!     xc_derivatives
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!****
module xc_derivatives

  use cp_log_handling,     only: cp_fatal_level, cp_failure_level,&
                                 cp_warning_level, cp_note_level, &
                                 cp_to_string, cp_logger_type
  use cp_error_handling,   only: cp_debug, cp_error_type, cp_error_init, &
                                 cp_error_dealloc_ref, cp_error_message, &
                                 cp_assert, cp_assertion_failed, &
                                 cp_internal_error, cp_a_l, cp_error_check
  use dft_types,           only: xc_control_type
  use kinds,               only: dp
  use global_types,        only: global_environment_type
  use timings,             only: timeset, timestop
  use xc_derivative_types, only: xc_derivative_set_type, &
                                 xc_derivative_type
  use xc_functionals,      only: xc_set, xc_lda, xc_lsd, &
                                 xc_get_pot_size

  implicit none

  private

  logical, parameter          :: debug_this_module=.false.
  character(len=*), parameter :: moduleN = 'xc_derivatives'

  public :: xc_calc_derivatives

!****************************************************************************

contains

  !!****f* cp2k/xc_derivatives/xc_calc_derivatives [1.0] *
  !!
  !!   NAME
  !!     xc_calc_derivatives
  !!
  !!   FUNCTION
  !!     -
  !!
  !!   ARGUMENTS
  !!     - derivative_set : will store the derivatives
  !!     - functional     : strings describing the functionals
  !!     - gradient_functionals:
  !!     - crossterms     :
  !!     - order          : up to which order to derive
  !!     - rhoa           : in spin-restricted calc the density (on a grid)
  !!                        in spin-unrestricted calc the alpha spin-density
  !!     - rhob, optional : in spin-unrestricted calc the beta spin-density
  !!     - drhoa, optional: in spin-restricted calc the grdient density (on a grid)
  !!                        in spin-unrestricted calc the gradient alpha spin-density
  !!     - drhob, optional: in spin-unrestricted calc the gradient beta spin-density
  !!     - error          : variable to control error logging, stopping,... 
  !!                        see module cp_error_handling 
  !!
  !!****
  subroutine xc_calc_derivatives(&
       derivative_set, &
       xc_control, &
       order, &
       rhoa, rhob, &
       drhoa, drhob, &
       error)

    implicit none

    ! arguments
    type(xc_derivative_set_type), pointer        :: derivative_set
    type(xc_control_type), intent(in)            :: xc_control
    integer, intent(in)                          :: order
    real(KIND = dp), dimension(:), intent(in)    :: rhoa
    real(KIND = dp), dimension(:), intent(in), optional :: rhob, drhoa, drhob
    type(cp_error_type), optional, intent(inout) :: error
    
    ! locals
    logical :: do_lsd
    logical :: failure
    integer :: npot, stat
    
    ! parameters
    character(len=*), parameter :: routineN = 'xc_calc_derivatives', &
                                   routineP = moduleN//'/'//routineN

    failure = .false.
    do_lsd = .false.
    if (present(rhob)) do_lsd = .true.

    !-------------!
    ! some checks !
    !-------------!
    CPPrecondition(associated(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(abs(order)<=3,cp_failure_level,routineP,error,failure)
    if (any(xc_control%gradient_functionals)) then
       CPPrecondition(present(drhoa),cp_failure_level,routineP,error,failure)
       if (do_lsd) then
          CPPrecondition(present(drhob),cp_failure_level,routineP,error,failure)
       end if
    end if
    
    if (.not.failure) then

       !-------------------------!
       ! allocate the data array !
       !-------------------------!
       if (.not.failure) then
          npot = xc_get_pot_size(do_lsd, any(xc_control%gradient_functionals), &
                                 do_lsd .and. any(xc_control%crossterms), order)
          ! for compatability 0:npot-1 later simply npot
          allocate(derivative_set%data(derivative_set%npoints,0:npot-1), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       end if
       
       !-------------------------!
       ! call the (old) routines !
       !-------------------------!
       call xc_set(xc_control%density_cut, xc_control%gradient_cut)
       if (any(xc_control%gradient_functionals)) then
          if (do_lsd) then
             call xc_lsd ( xc_control%functionals, &
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:), rhob(:), &
                           derivative_set, order, &
                           drho_a=drhoa(:), drho_b=drhob(:), &
                           error=error )
          else
             call xc_lda ( xc_control%functionals,&
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:),  &
                           derivative_set, order=order, &
                           drho=drhoa(:), &
                           error=error )
          end if
       else
          if (do_lsd) then
             call xc_lsd ( xc_control%functionals, &
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:), rhob(:), &
                           derivative_set, order, &
                           error=error )
          else
             call xc_lda ( xc_control%functionals,&
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:), &
                           derivative_set, order, &
                           error=error )
          end if
       end if
       
    end if

  end subroutine xc_calc_derivatives

end module xc_derivatives

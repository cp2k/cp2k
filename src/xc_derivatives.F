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
  use kinds,               only: wp=>dp
  use global_types,        only: global_environment_type
  use timings,             only: timeset, timestop
  use xc_derivative_types, only: xc_derivative_set_type, &
                                 xc_derivative_type
  use xc_functionals,      only: xc_set, xc_lda, xc_lsd

  implicit none

  private

  logical, parameter          :: debug_this_module=.false.
  character(len=*), parameter :: moduleN = 'xc_derivatives'

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
  !!     - density_cutoff :
  !!     - gradient_cutoff:
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
       functional, gradient_functionals, crossterms, &
       order, &
       density_cutoff, gradient_cutoff, &
       rhoa, rhob, &
       drhoa, drhob, &
       error)

    implicit none

    ! arguments
    type(xc_derivative_set_type), pointer        :: derivative_set
    character(len=40), dimension(3), intent(in)  :: functional
    logical, dimension(3), intent(in)            :: gradient_functionals, crossterms
    integer, intent(in)                          :: order
    real(kind=wp), intent(in)                    :: density_cutoff, gradient_cutoff
    real(kind=wp), dimension(:,:), intent(in)    :: rhoa
    real(kind=wp), dimension(:,:), intent(in), optional :: rhob, drhoa, drhob
    type(cp_error_type), optional, intent(inout) :: error
    
    ! locals
    logical :: do_lsd, gradient_functional
    logical :: failure
    integer :: npot
    real(kind=wp), dimension(:,:), pointer :: pot
    
    ! parameters
    character(len=*), parameter :: routineN = 'xc_calc_derivatives', &
                                   routineP = moduleN//'/'//routineN

    failure = .false.

    CPPrecondition(associated(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(abs(order)<3,cp_failure_level,routineP,error,failure)

    if (present(rhob)) do_lsd = .true.
    gradient_functional = ANY(gradient_functionals)

    if (gradient_functional) then
       CPPrecondition(present(drhoa),cp_failure_level,routineP,error,failure)
       if (do_lsd) then
          CPPrecondition(present(drhob),cp_failure_level,routineP,error,failure)
       end if
    end if
    
    if (.not.failure) then
       
       call xc_set(density_cutoff, gradient_cutoff)
       if (gradient_functional) then
          if (do_lsd) then
             CALL xc_lsd ( functional, gradient_functionals, crossterms, &
                           rhoa(:,1), rhob(:,1), &
                           pot, order, &
                           drho_a=drhoa(:,1), drho_b=drhob(:,1), &
                           error=error )
          else
             CALL xc_lsd ( functional, gradient_functionals, crossterms, &
                           rhoa(:,1), rhob(:,1), &
                           pot(:,1:npot), order, &
                           error=error )
          end if
       else
          if (gradient_functional) then
             call xc_lda ( functional,&
                           gradient_functionals, crossterms, &
                           rhoa(:,1),  &
                           pot, order=order, &
                           drho=drhoa(:,1), error=error )
          else
             call xc_lda ( functional,&
                           gradient_functionals, crossterms, &
                           rhoa(:,1), pot, &
                           order, error=error )
          end if
       end if

    end if

  end subroutine xc_calc_derivatives

end module xc_derivatives

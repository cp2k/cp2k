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
MODULE xc_derivatives
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE dft_types,                       ONLY: xc_control_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals,                  ONLY: xc_get_pot_size,&
                                             xc_lda,&
                                             xc_lsd,&
                                             xc_set
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_type
  USE xc_xbecke88,                     ONLY: xb88_lda_eval,&
                                             xb88_lda_info,&
                                             xb88_lsd_eval,&
                                             xb88_lsd_info
  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER          :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER :: moduleN = 'xc_derivatives'

  PUBLIC :: xc_calc_derivatives, xc_functional_get_info, xc_functional_eval

!****************************************************************************

CONTAINS

!!****f* cp2k/xc_derivatives/xc_calc_derivatives [1.0] *
!!
!!   NAME
!!     xc_calc_derivatives
!!
!!   SYNOPSIS
!!     Subroutine xc_calc_derivatives(derivative_set, xc_control, order,&
!!         rhoa, rhob, drhoa, drhob, error)
!!       Implicit None
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Type(xc_control_type), Intent (IN):: xc_control
!!       Integer, Intent (IN):: order
!!       Real(Kind=dp), Dimension(:), Intent (IN):: rhoa
!!       Real(Kind=dp), Dimension(:), Intent (IN), Optional:: rhob, drhoa,&
!!         drhob
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_calc_derivatives
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
  SUBROUTINE xc_calc_derivatives(&
       derivative_set, &
       xc_control, &
       order, &
       rhoa, rhob, &
       drhoa, drhob, &
       error)

    IMPLICIT NONE

    ! arguments
    TYPE(xc_derivative_set_type), POINTER        :: derivative_set
    TYPE(xc_control_type), INTENT(in)            :: xc_control
    INTEGER, INTENT(in)                          :: order
    REAL(KIND = dp), DIMENSION(:), INTENT(in)    :: rhoa
    REAL(KIND = dp), DIMENSION(:), INTENT(in), OPTIONAL :: rhob, drhoa, drhob
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    ! locals
    LOGICAL :: do_lsd
    LOGICAL :: failure
    INTEGER :: npot, stat
    
    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_derivatives', &
                                   routineP = moduleN//'/'//routineN

    failure = .FALSE.
    do_lsd = .FALSE.
    IF (PRESENT(rhob)) do_lsd = .TRUE.

    !-------------!
    ! some checks !
    !-------------!
    CPPrecondition(ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ABS(order)<=3,cp_failure_level,routineP,error,failure)
    IF (ANY(xc_control%gradient_functionals)) THEN
       CPPrecondition(PRESENT(drhoa),cp_failure_level,routineP,error,failure)
       IF (do_lsd) THEN
          CPPrecondition(PRESENT(drhob),cp_failure_level,routineP,error,failure)
       END IF
    END IF
    
    IF (.not.failure) THEN

       !-------------------------!
       ! allocate the data array !
       !-------------------------!
       IF (.not.failure) THEN
          npot = xc_get_pot_size(do_lsd, ANY(xc_control%gradient_functionals), &
                                 do_lsd .AND. ANY(xc_control%crossterms), order)
          ! for compatability 0:npot-1 later simply npot
          ALLOCATE(derivative_set%data(SIZE(rhoa),0:npot-1), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       
       !-------------------------!
       ! call the (old) routines !
       !-------------------------!
       CALL xc_set(xc_control%density_cut, xc_control%gradient_cut)
       IF (ANY(xc_control%gradient_functionals)) THEN
          IF (do_lsd) THEN
             CALL xc_lsd ( derivative_set, &
                           xc_control%functionals, &
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:), rhob(:), &
                           order, &
                           drho_a=drhoa(:), drho_b=drhob(:), &
                           error=error )
          ELSE
             CALL xc_lda ( derivative_set, &
                           xc_control%functionals,&
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:),  &
                           order=order, &
                           drho=drhoa(:), &
                           error=error )
          END IF
       ELSE
          IF (do_lsd) THEN
             CALL xc_lsd ( derivative_set, &
                           xc_control%functionals, &
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:), rhob(:), &
                           order, &
                           error=error )
          ELSE
             CALL xc_lda ( derivative_set, &
                           xc_control%functionals,&
                           xc_control%gradient_functionals, &
                           xc_control%crossterms, &
                           rhoa(:), &
                           order, &
                           error=error )
          END IF
       END IF
       
    END IF

  END SUBROUTINE xc_calc_derivatives

!!****f* xc_derivatives/xc_functional_get_info [1.0] *
!!
!!   NAME
!!     xc_functional_get_info
!!
!!   SYNOPSIS
!!     Subroutine xc_functional_get_info(functional, lsd, reference,&
!!         shortform, needs, error)
!!       Character(Len=*), Intent (IN):: functional
!!       Logical, Intent (IN):: lsd
!!       Character(Len=*), Intent (OUT), Optional:: reference
!!       Character(Len=*), Intent (OUT), Optional:: shortform
!!       Type(xc_rho_cflags_type), Intent (INOUT), Optional:: needs
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_functional_get_info
!!
!!   FUNCTION
!!     get the information about the given functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - functional: the functional you want info about
!!     - lsd: if you are using lsd or lda
!!     - reference: the reference to the acticle where the functional is
!!       explained
!!     - shortform: the short definition of the functional
!!     - needs: the flags corresponding to the inputs needed by this 
!!       functional are set to true (the flags not needed aren't touched)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xc_functional_get_info(functional,lsd,reference,shortform,&
     needs, error)
  CHARACTER (len=*), INTENT(in) :: functional
  LOGICAL, INTENT(in) :: lsd
  CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
  CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform
  TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL :: needs
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xc_functional_get_info',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     SELECT CASE(functional)
     CASE("BECKE88")
        IF (lsd) THEN
           CALL xb88_lsd_info(reference=reference,shortform=shortform,&
                needs=needs, error=error)
        ELSE
           CALL xb88_lda_info(reference=reference,shortform=shortform,&
                needs=needs, error=error)
        END IF
     CASE("BECKE88_LDA")
        CALL xb88_lda_info(reference=reference,shortform=shortform,&
             needs=needs, error=error)
     CASE("BECKE88_LSD")
        CALL xb88_lsd_info(reference=reference,shortform=shortform,&
             needs=needs, error=error)
     CASE("NONE")
        IF (PRESENT(reference)) reference=""
        IF (PRESENT(shortform)) shortform=""
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(functional)//"' in "//&
             CPSourceFileRef,&
             error=error)
     END SELECT
  END IF
END SUBROUTINE xc_functional_get_info
!***************************************************************************

!!****f* xc_derivatives/xc_functional_eval [1.0] *
!!
!!   NAME
!!     xc_functional_eval
!!
!!   SYNOPSIS
!!     Subroutine xc_functional_eval(functional, lsd, rho_set, deriv_set,&
!!         grad_deriv, error)
!!       Character(Len=*), Intent (IN):: functional
!!       Logical, Intent (IN):: lsd
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: grad_deriv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_functional_eval
!!
!!   FUNCTION
!!     evaluate a functional (and its derivatives)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - functional: a string that describes which functional should be used
!!     - lsd: if a local spin desnity is performed
!!     - rho_set: a rho set where all the arguments needed by this functional
!!       should be valid (which argument are needed acn be found with
!!       xc_functional_get_info)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xc_functional_eval(functional, lsd, rho_set,deriv_set,&
     grad_deriv,error)
  CHARACTER(len=*),INTENT(in) :: functional
  LOGICAL, INTENT(in) :: lsd
  TYPE(xc_rho_set_type), POINTER :: rho_set
  TYPE(xc_derivative_set_type), POINTER :: deriv_set
  INTEGER, INTENT(in) :: grad_deriv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xc_functional_eval',&
       routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     SELECT CASE(functional)
     CASE("BECKE88")
        IF (lsd) THEN
           CALL xb88_lsd_eval(rho_set=rho_set,deriv_set=deriv_set,&
                grad_deriv=grad_deriv,error=error)
        ELSE
           CALL xb88_lda_eval(rho_set=rho_set,deriv_set=deriv_set,&
                grad_deriv=grad_deriv,error=error)
        END IF
     CASE("BECKE88_LDA")
        CALL xb88_lda_eval(rho_set=rho_set,deriv_set=deriv_set,&
                grad_deriv=grad_deriv,error=error)
     CASE("BECKE88_LSD")
        CALL xb88_lsd_eval(rho_set=rho_set,deriv_set=deriv_set,&
                grad_deriv=grad_deriv,error=error)
     CASE("NONE")
        ;
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(functional)//"' in "//&
             CPSourceFileRef,&
             error=error)
     END SELECT
  END IF
END SUBROUTINE xc_functional_eval
!***************************************************************************

END MODULE xc_derivatives

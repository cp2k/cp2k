!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_optx [1.0] *
!!
!!   NAME
!!     xc_optx
!!
!!   FUNCTION
!!     calculate optx
!!
!!   NOTES
!!     will need proper testing / review
!!
!!   AUTHOR
!!     Joost VandeVondele [03.2004]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_optx
  USE cp_error_handling,               ONLY: cp_a_l,&
       cp_assert,&
       cp_debug,&
       cp_error_get_logger,&
       cp_error_message,&
       cp_error_type,&
       cp_unimplemented_error, cp_assertion_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
       cp_fatal_level,&
       cp_note_level,&
       cp_to_string,&
       cp_warning_level
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE timings,                         ONLY: timeset,&
       timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
       xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
       xc_derivative_type
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
       xc_rho_set_get,&
       xc_rho_set_type

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_optx'

  PUBLIC :: optx_lda_info, optx_lda_eval
  !!***
  !****************************************************************************
CONTAINS
  !!****f* xc_optx/optx_lda_info [1.0] *
  !!
  !!   NAME
  !!     optx_lda_info
  !!
  !!   FUNCTION
  !!     info about the optx functional
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - reference: string with the reference of the actual functional
  !!     - shortform: string with the shortform of the functional name
  !!     - needs: the components needed by this functional are set to
  !!       true (does not set the unneeded components to false)
  !!     - maximum implemented derivative of the xc functional
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Joost
  !!
  !!*** **********************************************************************
  SUBROUTINE optx_lda_info(reference,shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
         INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'optx_lda_info', &
         routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "OPTX, Handy NC and Cohen AJ,  JCP 116, p. 5411 (2002), btw April 1st"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
          shortform = "OPTX exchange (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1
  END SUBROUTINE optx_lda_info
  !***************************************************************************
  !!****f* xc_optx/optx_lda_eval [1.0] *
  !!
  !!   NAME
  !!     optx_lda_eval
  !!
  !!   FUNCTION
  !!     evaluates the optx functional for lda
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - rho_set: the density where you want to evaluate the functional
  !!     - deriv_set: place where to store the functional derivatives (they are
  !!       added to the derivatives)
  !!     - grad_deriv: degree of the derivative that should be evalated,
  !!       if positive all the derivatives up to the given degree are evaluated,
  !!       if negative only the given degree is calculated
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Joost 
  !!
  !!*** **********************************************************************
  SUBROUTINE optx_lda_eval(rho_set,deriv_set,grad_deriv,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'optx_lda_eval', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, &
         e_rho, norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    failure=.FALSE.
    NULLIFY(bo,e_0, e_ndrho, e_rho, norm_drho, rho)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
            drho_cutoff=epsilon_drho,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       IF (grad_deriv>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
       IF (grad_deriv>1.OR.grad_deriv<-1) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 1 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL optx_lda_calc(rho=rho, norm_drho=norm_drho,&
               e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
               npoints=npoints,epsilon_rho=epsilon_rho,&
               epsilon_drho=epsilon_drho, error=error)
    END IF
  END SUBROUTINE optx_lda_eval
  !***************************************************************************

  !!****f* xc_optx/optx_lda_calc *
  !!
  !!   NAME
  !!     optx_lda_calc
  !!
  !!   FUNCTION
  !!     optx exchange functional
  !!
  !!   NOTES
  !!
  !!   ARGUMENTS
  !!     - rho: the density
  !!     - norm_drho: the norm of the gradient of the density
  !!     - e_0: the value of the functional in that point
  !!     - e_rho: the derivative of the functional wrt. rho
  !!     - e_ndrho: the derivative of the functional wrt. norm_drho
  !!     - epsilon_rho: the cutoff on rho
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Joost VandeVondele
  !!
  !!   MODIFICATION HISTORY
  !!
  !!*** *********************************************************************
  SUBROUTINE optx_lda_calc(rho,norm_drho,e_0,e_rho,e_ndrho,&
       epsilon_rho,epsilon_drho,npoints,error)
    REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho, norm_drho
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rho, e_ndrho
    REAL(kind=dp), INTENT(in) :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(in) :: npoints
    TYPE(cp_error_type), INTENT(inout), optional :: error
 
    REAL(KIND=dp), PARAMETER :: o43=4.0_dp/3.0_dp ,gam=0.006_dp, &
                               a1cx=0.9784571170284421_dp,a2=1.43169_dp
    REAL(KIND=dp) :: myrho,myndrho,xs,denom,rho43,ex,gamxsxs,tmp
    INTEGER :: ii
    DO ii=1,npoints
       ! we get the full density and need spin parts -> 0.5
       myrho   = 0.5_dp * rho(ii)
       myndrho = 0.5_dp * MAX(norm_drho(ii),epsilon_drho)
       IF (myrho>epsilon_rho) THEN
          rho43   = myrho**o43
          xs      = (myndrho / rho43)
          gamxsxs = gam * xs * xs
          denom   = 1.0_dp / (1.0_dp + gamxsxs )
          ex      = rho43*(a1cx+a2*(gamxsxs*denom)**2)
          ! 2.0 for both spins
          e_0(ii)    = e_0(ii)     - 2.0_dp * ex
          tmp = rho43 * 2.0_dp * a2 * gamxsxs * denom**2 * ( 1.0_dp - gamxsxs * denom) 
          ! derive e_0 wrt to rho (full) and ndrho (also full)
          e_rho(ii)  = e_rho(ii)   - ( o43 * ex / myrho + tmp * gamxsxs * (-2.0_dp * o43 ) / myrho  )
          e_ndrho(ii)= e_ndrho(ii) - ( tmp * gam * 2.0_dp * myndrho / rho43**2 )
       END IF
    END DO

  END SUBROUTINE optx_lda_calc
  !***************************************************************************
END MODULE xc_optx

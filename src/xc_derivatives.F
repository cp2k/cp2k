!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
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
  USE cp_control_types,                ONLY: xc_control_type
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
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_cs1_new,                      ONLY: cs1_lda_eval,&
                                             cs1_lda_info,&
                                             cs1_lsd_eval,&
                                             cs1_lsd_info
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_exchange_gga_new,             ONLY: xgga_b88,&
                                             xgga_eval,&
                                             xgga_info,&
                                             xgga_opt,&
                                             xgga_pw86,&
                                             xgga_pw91,&
                                             xgga_revpbe
  USE xc_functionals,                  ONLY: xc_lda,&
                                             xc_lsd,&
                                             xc_set
  USE xc_hcth,                         ONLY: hcth_lda_eval,&
                                             hcth_lda_info
  USE xc_ke_gga_new,                   ONLY: ke_gga_info,&
                                             ke_gga_lda_eval,&
                                             ke_gga_lsd_eval,&
                                             ke_lc,&
                                             ke_llp,&
                                             ke_ol1,&
                                             ke_ol2,&
                                             ke_pw86,&
                                             ke_pw91,&
                                             ke_t92
  USE xc_lyp_new,                      ONLY: lyp_lda_eval,&
                                             lyp_lda_info,&
                                             lyp_lsd_eval,&
                                             lyp_lsd_info
  USE xc_optx,                         ONLY: optx_lda_eval,&
                                             optx_lda_info,&
                                             optx_lsd_eval,&
                                             optx_lsd_info
  USE xc_pade,                         ONLY: pade,&
                                             pade_lda_pw_eval,&
                                             pade_lsd_pw_eval
  USE xc_pbe,                          ONLY: pbe_lda_eval,&
                                             pbe_lda_info,&
                                             pbe_lsd_eval,&
                                             pbe_lsd_info
  USE xc_perdew86_new,                 ONLY: p86_lda_eval,&
                                             p86_lda_info
  USE xc_perdew_wang_new,              ONLY: perdew_wang_info,&
                                             perdew_wang_lda_eval,&
                                             perdew_wang_lsd_eval
  USE xc_perdew_zunger_new,            ONLY: pz_info,&
                                             pz_lda_eval,&
                                             pz_lsd_eval
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type
  USE xc_tfw_new,                      ONLY: tfw_lda_eval,&
                                             tfw_lda_info,&
                                             tfw_lsd_eval,&
                                             tfw_lsd_info
  USE xc_thomas_fermi_new,             ONLY: thomas_fermi_info,&
                                             thomas_fermi_lda_eval,&
                                             thomas_fermi_lsd_eval
  USE xc_tpss,                         ONLY: my_tpss_lda_eval => tpss_lda_eval,&
                                             my_tpss_lda_info => tpss_lda_info,&
                                             my_tpss_lsd_eval => tpss_lsd_eval,&
                                             my_tpss_lsd_info => tpss_lsd_info
  USE xc_tpss_p,                       ONLY: tpss_lsd_eval,&
                                             tpss_lsd_info
  USE xc_tpss_u,                       ONLY: tpss_lda_eval,&
                                             tpss_lda_info
  USE xc_vwn_new,                      ONLY: vwn_lda_eval,&
                                             vwn_lda_info
  USE xc_xalpha_new,                   ONLY: xalpha_info,&
                                             xalpha_lda_eval,&
                                             xalpha_lsd_eval
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


    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    TYPE(xc_control_type), INTENT(in)        :: xc_control
    INTEGER, INTENT(in)                      :: order
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: rhoa
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: rhob, drhoa, drhob
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_derivatives', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: do_lsd, failure

    CALL timeset(routineN,"I","",handle)

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

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xc_calc_derivatives

!!****f* xc_derivatives/xc_functional_get_info [1.0] *
!!
!!   NAME
!!     xc_functional_get_info
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
!!     - max_deriv: the maximal derivative available
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE xc_functional_get_info(functional,lsd,reference,shortform,&
     needs, max_deriv, error)
    CHARACTER(len=*), INTENT(in)             :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_functional_get_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     SELECT CASE(functional)
     CASE("BECKE88")
        IF (lsd) THEN
           CALL xb88_lsd_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        ELSE
           CALL xb88_lda_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        END IF
     CASE("LYP")
        IF (lsd) THEN
           CALL lyp_lsd_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        ELSE
           CALL lyp_lda_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        END IF
     CASE("PADE")
        CALL pade(reference, shortform, lsd=lsd, needs=needs, error=error)
     CASE("HCTH93")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(93,reference,shortform,needs,max_deriv,error=error)
     CASE("HCTH120")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(120,reference,shortform,needs,max_deriv,error=error)
     CASE("HCTH147")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(147,reference,shortform,needs,max_deriv,error=error)
     CASE("HCTH407")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(407,reference,shortform,needs,max_deriv,error=error)
     CASE("OPTX")
        IF (lsd) THEN
           CALL optx_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL optx_lda_info(reference,shortform,needs,max_deriv,error=error)
        ENDIF
     CASE("CS1")
        IF (lsd) THEN
           CALL cs1_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL cs1_lda_info(reference, shortform, needs=needs, max_deriv=max_deriv,&
                error=error)
        ENDIF
     CASE("BECKE88X")
        CALL xgga_info(xgga_b88,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PW86X")
        CALL xgga_info(xgga_pw86,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PW91X")
        CALL xgga_info(xgga_pw91,lsd,reference,shortform,needs,max_deriv,error=error)
!     CASE("PBEX")
!        CALL xgga_info(xgga_pbe,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("REV_PBEX")
        CALL xgga_info(xgga_revpbe,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("OPTX_2")
        CALL xgga_info(xgga_opt,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("OL1","KE_OL1")
        CALL ke_gga_info(ke_ol1,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("OL2","KE_OL2")
        CALL ke_gga_info(ke_ol2,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("LLP","KE_LLP")
        CALL ke_gga_info(ke_llp,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("KE_PW86")
        CALL ke_gga_info(ke_pw86,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("KE_PW91")
        CALL ke_gga_info(ke_pw91,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("LC","KE_LC")
        CALL ke_gga_info(ke_LC,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("T92","KE_T92")
        CALL ke_gga_info(ke_t92,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("P86C")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL p86_lda_info(reference,shortform,needs,max_deriv,error=error)
     CASE("PW","PW92","PW92_PWO")
        CALL perdew_wang_info("PWO",lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PW92_DMC")
        CALL perdew_wang_info("DMC",lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PW92_VMC")
        CALL perdew_wang_info("VMC",lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PZ81_ORIG")
        CALL pz_info("ORIG",lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PZ","PZ81_DMC")
        CALL pz_info("DMC",lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PZ81_VMC")
        CALL pz_info("VMC",lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("TFW")
        IF (lsd) THEN
           CALL tfw_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL tfw_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE("TF")
        CALL thomas_fermi_info(lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("VWN")
        CPPrecondition(.NOT.lsd,cp_failure_level,routineP,error,failure)
        CALL vwn_lda_info(reference,shortform,needs,max_deriv,error=error)
     CASE("XALPHA","SLATER")
        !FM to use an input defined xa_parameter we should probably switch
        !FM to a derived type instead of a simple string as mean to select
        !FM the functional, which is probably a good idea, but for the moment
        !FM left for the future.
        CALL xalpha_info(lsd,reference,shortform,needs,max_deriv,&
             xa_parameter=2._dp/3._dp,error=error)
     CASE("TPSS")
        IF (lsd) THEN
           CALL tpss_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL tpss_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE("MY_TPSS")
        IF (lsd) THEN
           CALL my_tpss_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL my_tpss_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE("PBEX")
        IF (PRESENT(reference)) reference=""
        IF (PRESENT(shortform)) shortform=""
        IF (PRESENT(max_deriv)) max_deriv=HUGE(0)
     CASE("PBEC")
        IF (lsd) THEN
           CALL pbe_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL pbe_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE("NONE")
        IF (PRESENT(reference)) reference=""
        IF (PRESENT(shortform)) shortform=""
        IF (PRESENT(max_deriv)) max_deriv=HUGE(0)
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
!!       should be valid (which argument are needed can be found with
!!       xc_functional_get_info)
!!     - deriv_set: place where to store the functional derivatives (they are
!!       added to the derivatives)
!!     - deriv_order: degree of the derivative that should be evalated,
!!       if positive all the derivatives up to the given degree are evaluated,
!!       if negative only the given degree is requested (but to simplify 
!!       the code all the derivatives might be calculated, you should ignore
!!       them: when adding derivatives of various functionals they might contain
!!       the derivative of just one functional)
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
SUBROUTINE xc_functional_eval(functional, lsd, rho_set, deriv_set, &
     deriv_order, error)

    CHARACTER(len=*), INTENT(in)             :: functional
    LOGICAL, INTENT(in)                      :: lsd
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: deriv_order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_functional_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: density_cut, gradient_cut

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.

  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,rho_cutoff=density_cut,&
          drho_cutoff=gradient_cut,error=error)

     SELECT CASE(functional)
     CASE("BECKE88")
        IF (lsd) THEN
           CALL xb88_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL xb88_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("LYP")
        IF (lsd) THEN
           CALL lyp_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL lyp_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PADE")
        CALL pade(density_cut)
        IF (lsd) THEN
           CALL pade_lsd_pw_eval(deriv_set, rho_set, deriv_order, error=error)
        ELSE
           CALL pade_lda_pw_eval(deriv_set, rho_set, deriv_order, error=error)
        END IF
     CASE("HCTH93")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(93,rho_set,deriv_set,deriv_order,error=error)
     CASE("HCTH120")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(120,rho_set,deriv_set,deriv_order,error=error)
     CASE("HCTH147")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(147,rho_set,deriv_set,deriv_order,error=error)
     CASE("HCTH407")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(407,rho_set,deriv_set,deriv_order,error=error)
     CASE("OPTX")
        IF (lsd) THEN
           CALL optx_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL optx_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        ENDIF
     CASE("CS1")
        IF (lsd) THEN
           CALL cs1_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL cs1_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        ENDIF
     CASE("BECKE88X")
        CALL xgga_eval(xgga_b88,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE("PW86X")
        CALL xgga_eval(xgga_pw86,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE("PW91X")
        CALL xgga_eval(xgga_pw91,lsd,rho_set,deriv_set,deriv_order,error=error)
!     CASE("PBEX")
!        CALL xgga_eval(xgga_pbe,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE("REV_PBEX")
        CALL xgga_eval(xgga_revpbe,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE("OPTX_2")
        CALL xgga_eval(xgga_opt,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE("OL1","KE_OL1")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_ol1,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_ol1,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("OL2","KE_OL2")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_ol2,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_ol2,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("LLP","KE_LLP")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_llp,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_llp,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("KE_PW86")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_pw86,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_pw86,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("KE_PW91")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_pw91,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_pw91,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("LC","KE_LC")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_LC,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_LC,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("T92","KE_T92")
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_t92,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_t92,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("P86","P86C")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL p86_lda_eval(rho_set,deriv_set,deriv_order,error=error)
     CASE("PW","PW92","PW92_PWO")
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval("PWO",rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval("PWO",rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PW92_DMC")
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval("DMC",rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval("DMC",rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PW92_VMC")
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval("VMC",rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval("VMC",rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PZ81_ORIG")
        IF (lsd) THEN
           CALL pz_lsd_eval("ORIG",rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval("ORIG",rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PZ","PZ81_DMC")
        IF (lsd) THEN
           CALL pz_lsd_eval("DMC",rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval("DMC",rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PZ81_VMC")
        IF (lsd) THEN
           CALL pz_lsd_eval("VMC",rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval("VMC",rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("TFW")
        IF (lsd) THEN
           CALL tfw_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL tfw_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE ("TF")
        IF (lsd) THEN
           CALL thomas_fermi_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL thomas_fermi_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF        
     CASE("VWN")
        CPPrecondition(.NOT.lsd,cp_failure_level,routineP,error,failure)
        CALL vwn_lda_eval(rho_set,deriv_set,deriv_order,error=error)
     CASE("XALPHA","SLATER")
        !FM to use an input defined xa_parameter we should probably switch
        !FM to a derived type instead of a simple string as mean to select
        !FM the functional, which is probably a good idea, but for the moment
        !FM left for the future.
        IF (lsd) THEN
           CALL xalpha_lsd_eval(rho_set,deriv_set,deriv_order,&
                xa_parameter=2._dp/3._dp,error=error)
        ELSE
           CALL xalpha_lda_eval(rho_set,deriv_set,deriv_order,&
                xa_parameter=2._dp/3._dp,error=error)
        END IF
     CASE("TPSS")
        IF (lsd) THEN
           CALL tpss_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL tpss_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("MY_TPSS")
        IF (lsd) THEN
           CALL my_tpss_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL my_tpss_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PBEX")
     CASE("PBEC")
        IF (lsd) THEN
           CALL pbe_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pbe_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("NONE")
        !
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(functional)//"' in "//&
             CPSourceFileRef,&
             error=error)
     END SELECT
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE xc_functional_eval
!***************************************************************************

END MODULE xc_derivatives

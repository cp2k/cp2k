!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_cp2k_dft, only: xc_functionals_expand
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_cs1_new,                      ONLY: cs1_lda_eval,&
                                             cs1_lda_info,&
                                             cs1_lsd_eval,&
                                             cs1_lsd_info
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_exchange_gga_new,             ONLY: xgga_eval,&
                                             xgga_info
  USE xc_hcth,                         ONLY: hcth_lda_eval,&
                                             hcth_lda_info
  USE xc_ke_gga_new,                   ONLY: ke_gga_info,&
                                             ke_gga_lda_eval,&
                                             ke_gga_lsd_eval
  USE input_constants,                 ONLY: xgga_b88,&
                                             xgga_b88x,&  
                                             xgga_opt,&
                                             xgga_pbex,&
                                             xgga_pw86,&
                                             xgga_pw91,&
                                             xgga_revpbe,&
                                             ke_lc,&
                                             ke_llp,&
                                             ke_ol1,&
                                             ke_ol2,&
                                             ke_pw86,&
                                             ke_pw91,&
                                             ke_pbe,&
                                             ke_revpbe,&
                                             ke_t92,&
                                             xc_funct_no_shortcut,&
                                             xc_funct_blyp,&
                                             xc_funct_pade,&
                                             xc_funct_pbe,&
                                             xc_funct_tpss,&
                                             xc_funct_hcth120
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
  USE input_constants,                 ONLY: pw_dmc,&
                                             pw_orig,&
                                             pw_vmc
  USE input_constants,                 ONLY: pz_dmc,&
                                             pz_orig,&
                                             pz_vmc
  USE xc_perdew_zunger_new,            ONLY: pz_info,&
                                             pz_lda_eval,&
                                             pz_lsd_eval
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type,&
       xc_rho_cflags_setall
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
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
  USE xc_vwn_new,                      ONLY: vwn_lda_eval,&
                                             vwn_lda_info
  USE xc_xalpha_new,                   ONLY: xalpha_info,&
                                             xalpha_lda_eval,&
                                             xalpha_lsd_eval
  USE xc_xbecke88,                     ONLY: xb88_lda_eval,&
                                             xb88_lda_info,&
                                             xb88_lsd_eval,&
                                             xb88_lsd_info
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER          :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER :: moduleN = 'xc_derivatives'

  PUBLIC :: xc_functional_get_info, xc_functionals_eval, xc_functionals_get_needs

!****************************************************************************

CONTAINS

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
    TYPE(section_vals_type), POINTER         :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_functional_get_info_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_param
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: r_param

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(functional),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     SELECT CASE(functional%section%name)
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
     CASE("HCTH")
        CALL section_vals_val_get(functional,"PARAMETER_SET",i_val=i_param,&
             error=error)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(i_param,reference,shortform,needs,max_deriv,error=error)
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
     CASE("XGGA")
        CALL section_vals_val_get(functional,"FUNCTIONAL",i_val=i_param,&
             error=error)
        CALL xgga_info(i_param,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("KE_GGA")
        CALL section_vals_val_get(functional,"FUNCTIONAL",i_val=i_param,&
             error=error)
        CALL ke_gga_info(i_param,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("P86C")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        IF(failure) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "BP functional not implemented with LSD")
        END IF
        CALL p86_lda_info(reference,shortform,needs,max_deriv,error=error)
     CASE("PW92")
        CALL section_vals_val_get(functional,"PARAMETRIZATION",i_val=i_param,&
             error=error)
        CALL perdew_wang_info(i_param,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE("PZ81")
        CALL section_vals_val_get(functional,"PARAMETRIZATION",i_val=i_param,&
             error=error)
        CALL pz_info(i_param,lsd,reference,shortform,needs,max_deriv,error=error)
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
     CASE("XALPHA")
        CALL section_vals_val_get(functional,"XA",r_val=r_param,&
             error=error)
        CALL xalpha_info(lsd,reference,shortform,needs,max_deriv,&
             xa_parameter=r_param,error=error)
     CASE("TPSS")
        IF (lsd) THEN
           CALL my_tpss_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL my_tpss_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE("PBE")
        IF (lsd) THEN
           CALL pbe_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL pbe_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(functional%section%name)//&
             "' in "//&
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
!!     - functional: a section that describes the functional to be added
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

    TYPE(section_vals_type), POINTER         :: functional
    LOGICAL, INTENT(in)                      :: lsd
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: deriv_order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_functional_eval_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_param
    LOGICAL                                  :: failure, fun_active
    REAL(KIND=dp)                            :: density_cut, gradient_cut, &
                                                r_param

  CALL timeset(routineN,"I"," ",handle)
  failure=.FALSE.

  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,rho_cutoff=density_cut,&
          drho_cutoff=gradient_cut,error=error)
     CALL section_vals_val_get(functional,"_SECTION_PARAMETERS_",&
          l_val=fun_active,error=error)
     IF (.NOT.fun_active) THEN
        go to 123
     END IF

     SELECT CASE(functional%section%name)
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
     CASE("HCTH")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL section_vals_val_get(functional,"PARAMETER_SET",i_val=i_param,&
             error=error)
        CALL hcth_lda_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
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
     CASE("XGGA")
        CALL section_vals_val_get(functional,"FUNCTIONAL",i_val=i_param,&
             error=error)
        CALL xgga_eval(i_param,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE("KE_GGA")
        CALL section_vals_val_get(functional,"FUNCTIONAL",i_val=i_param,&
             error=error)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("P86C")
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL p86_lda_eval(rho_set,deriv_set,deriv_order,error=error)
     CASE("PW92")
        CALL section_vals_val_get(functional,"PARAMETRIZATION",i_val=i_param,&
             error=error)
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PZ81")
        CALL section_vals_val_get(functional,"PARAMETRIZATION",i_val=i_param,&
             error=error)
        IF (lsd) THEN
           CALL pz_lsd_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval(i_param,rho_set,deriv_set,deriv_order,error=error)
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
     CASE("XALPHA")
        CALL section_vals_val_get(functional,"XA",r_val=r_param,&
             error=error)
        IF (lsd) THEN
           CALL xalpha_lsd_eval(rho_set,deriv_set,deriv_order,&
                xa_parameter=r_param,error=error)
        ELSE
           CALL xalpha_lda_eval(rho_set,deriv_set,deriv_order,&
                xa_parameter=r_param,error=error)
        END IF
     CASE("TPSS")
        IF (lsd) THEN
           CALL my_tpss_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL my_tpss_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE("PBE")
        IF (lsd) THEN
           CALL pbe_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pbe_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(functional%section%name)//&
             "' in "//&
             CPSourceFileRef,&
             error=error)
     END SELECT
  END IF
123 continue
  CALL timestop(0.0_dp,handle)
END SUBROUTINE xc_functional_eval
!***************************************************************************

!!****f* xc_derivatives/xc_functionals_eval *
!!
!!   NAME
!!     xc_functionals_eval
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - functionals: a section containing the functional combination to be
!!       applied
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
!!*** **********************************************************************
SUBROUTINE xc_functionals_eval(functionals, lsd, rho_set, deriv_set, &
     deriv_order, error)
  TYPE(section_vals_type), pointer :: functionals
  LOGICAL, INTENT(in) :: lsd
  TYPE(xc_rho_set_type), POINTER :: rho_set
  TYPE(xc_derivative_set_type), pointer :: deriv_set
  INTEGER, intent(in) :: deriv_order
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='xc_functionals_eval',&
        routineP=moduleN//':'//routineN
  integer :: ifun
  TYPE(section_vals_type), pointer :: xc_fun

  failure=.false.
  
  CPPrecondition(ASSOCIATED(functionals),cp_failure_level,routineP,error,failure)
  CALL xc_functionals_expand(functionals,error=error)
  IF (.NOT. failure) THEN
     ifun=0
     DO 
        ifun=ifun+1
        xc_fun => section_vals_get_subs_vals2(functionals,i_section=ifun,&
             error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        CALL xc_functional_eval(xc_fun, &
             lsd=lsd,&
             rho_set=rho_set, &
             deriv_set=deriv_set,&
             deriv_order=deriv_order, &
             error=error)
     END DO
  END IF
END SUBROUTINE xc_functionals_eval
!***************************************************************************

!!****f* xc_derivatives/xc_functionals_get_needs *
!!
!!   NAME
!!     xc_functionals_get_needs
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - functionals: a section containing the functional combination to be
!!       applied
!!     - lsd: if a local spin desnity is performed
!!     - add_basic_components: makes sure that if some gradient combination
!!       is needed also the gradient components are requested, and if lsd
!!       is true rho_spin is requested. Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION xc_functionals_get_needs(functionals, lsd, add_basic_components,&
     error) RESULT(needs)
  TYPE(section_vals_type), pointer :: functionals
  LOGICAL, INTENT(in) :: lsd
  LOGICAL, INTENT(in), OPTIONAL :: add_basic_components
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(xc_rho_cflags_type) :: needs
  
  LOGICAL :: failure, my_add_basic_components
  character(len=*), parameter :: routineN='xc_functionals_get_needs',&
        routineP=moduleN//':'//routineN
  integer :: ifun
  TYPE(section_vals_type), POINTER :: xc_fun

  failure=.false.
  my_add_basic_components=.false.
  if (present(add_basic_components)) my_add_basic_components=add_basic_components
  
  CPPrecondition(ASSOCIATED(functionals),cp_failure_level,routineP,error,failure)
  CALL xc_functionals_expand(functionals,error=error)
  IF (.NOT. failure) THEN
     CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
     ifun=0
     DO
        ifun=ifun+1
        xc_fun => section_vals_get_subs_vals2(functionals,i_section=ifun,&
             error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        CALL xc_functional_get_info(xc_fun,lsd=lsd,&
             needs=needs, error=error)
     END DO
     IF (my_add_basic_components) THEN
        IF (lsd) THEN
           needs%rho_spin=.true.
           needs%tau_spin=needs%tau_spin.OR.needs%tau
        ELSE
           needs%rho=.true.
        END IF
        IF (needs%drho .OR. needs%drho_spin .OR. needs%norm_drho .OR.&
             needs%norm_drho_spin .OR. needs%drhoa_drhob) THEN
           IF (lsd) THEN
              needs%drho_spin=.TRUE.
           ELSE
              needs%drho=.TRUE.
           END IF
        END IF
     END IF
  END IF
END FUNCTION xc_functionals_get_needs
!***************************************************************************

END MODULE xc_derivatives

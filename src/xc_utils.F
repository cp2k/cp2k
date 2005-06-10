!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/xc_utils [1.0] *
!!
!!   NAME
!!     xc_utils
!!
!!   FUNCTION
!!     Utilities to set up the xc functional
!!     read xc_control from input, get the info about the needs of the
!!     required functionals, and the driver to the actual calculation of
!!     energy and potentials once rho and its derivatives are given
!!
!!   SOURCE
!******************************************************************************

MODULE xc_utils

  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_reset,&
                                             cp_error_type
  USE cp_keywords,                     ONLY: find_keyword,&
                                             keyword_type,&
                                             list_keywords
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: cgga_lyp, c_pade, &
                                             cgga_hcth93, cgga_hcth120, &
                                             cgga_hcth147, cgga_hcth407,&
                                             cgga_cs1, cgga_pw86, c_pw92,&
                                             c_pw92dmc, c_pw92vmc,c_pz, &
                                             c_pzdmc, c_pzvmc, c_vwn,&
                                             cgga_pbe
  USE input_constants,                 ONLY: ke_ol1, ke_ol2,ke_llp,&
                                             ke_pw86, ke_pw91, ke_lc,&
                                             ke_t92, ke_pbe,ke_revpbe, &
                                             ke_tf, ke_tfw, xgga_b88,&
                                             xgga_b88x, xgga_pw86,&
                                             xgga_pw91,xgga_pbex,&
                                             xgga_revpbe, xgga_optx,&
                                             xgga_hcth, xgga_tpss, &
                                             xalpha, xgga_pbe, xgga_opt
  USE input_constants,                 ONLY: xc_new_f_routine,&
                                             xc_old_f_routine,&
                                             xc_test_lsd_f_routine,&
                                             xc_debug_new_routine
  USE input_cp2k,                      ONLY: parsed_cp2k_input
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type
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
   
                                             
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "xc_utils"

! *** Public subroutines ***

  PUBLIC :: read_xc_control, xc_get_info_yat!, xc_eval_yat
  
!!***
! *****************************************************************************

CONTAINS

  SUBROUTINE read_xc_control(xc_control,nspins,kenergy_den,globenv,error)

    TYPE(xc_control_type)                    :: xc_control
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(INOUT)                   :: kenergy_den
    TYPE(global_environment_type)            :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER ::  routineN = "read_xc_control", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: i,i_param,istat
    LOGICAL                                  :: failure, was_present, xc_fun_set
    TYPE(section_vals_type), POINTER         :: xc_section, &
                                                func_section,&
                                                grid_section,&
                                                tmp_section
    LOGICAL, PARAMETER                       :: CROSS_TERM = .TRUE., &
                                                GRADIENT_FUNCTIONAL = .TRUE.

    failure=.FALSE.
    was_present = .FALSE.

    NULLIFY(xc_section,func_section,grid_section)
    xc_section => section_vals_get_subs_vals(globenv%input_file,&
                 "FORCE_EVAL",error=error)
    xc_section => section_vals_get_subs_vals(xc_section,"DFT",error=error)
    xc_section => section_vals_get_subs_vals(xc_section,"XC",error=error)

    
    func_section => section_vals_get_subs_vals(xc_section,&
                    "XC_FUNCTIONAL",error=error)
    grid_section => section_vals_get_subs_vals(xc_section,&
                    "XC_GRID",error=error)

    xc_control%functionals = "NONE" 
    xc_control%functionals_new = 0

    CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",&
         r_val=xc_control%density_cut,error=error)
    CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",&
         r_val=xc_control%gradient_cut,error=error)
    xc_control%density_smooth_cut_range = 0.0_dp
    CALL section_vals_val_get(xc_section,"TAU_CUTOFF",&
         r_val=xc_control%tau_cut,error=error)
    xc_control%gradient_functionals = .FALSE.
    xc_control%crossterms = .FALSE.
    CALL section_vals_val_get(xc_section,"FUNCTIONAL_ROUTINE",&
         i_val=xc_control%functional_routine,error=error)

    CALL xc_rho_cflags_setall(xc_control%needs,.FALSE.,error=error)    

    CALL section_vals_get(func_section,explicit=was_present,error=error)
    if (was_present) then
      tmp_section => section_vals_get_subs_vals(func_section,"BECKE88",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(1)=tmp_section%section%name
         xc_control%gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
         xc_control%functionals_new(1)= xgga_b88
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"LYP",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(2)=tmp_section%section%name
         xc_control%gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
         xc_control%crossterms ( 2 ) = CROSS_TERM
         xc_control%functionals_new(2)=cgga_lyp
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"PADE",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(1)="NONE"
         xc_control%functionals(2)=tmp_section%section%name
         xc_control%crossterms ( 2 ) = CROSS_TERM
         xc_control%functionals_new(2)=c_pade
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"HCTH",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         CALL section_vals_val_get(tmp_section,"PARAMETER_SET",i_val=i_param,error=error) 
         xc_control%functionals(1)=tmp_section%section%name//"X"
         xc_control%functionals(2)=tmp_section%section%name//TRIM(ADJUSTL(cp_to_string(i_param)))
         xc_control%gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
         xc_control%crossterms ( 2 ) = CROSS_TERM
         xc_control%functionals_new(1)=xgga_hcth
         SELECT CASE (i_param)
         CASE(93)
            xc_control%functionals_new(2) = cgga_hcth93
         CASE(120)
            xc_control%functionals_new(2) = cgga_hcth120
         CASE(147)
            xc_control%functionals_new(2) = cgga_hcth147
         CASE(407)
            xc_control%functionals_new(2) = cgga_hcth407
         CASE DEFAULT
           CALL stop_program ( 'read_xc_functional', 'Not known HCTH parameter' )
         END SELECT
         
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"OPTX",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(1)=tmp_section%section%name
         xc_control%gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
         xc_control%functionals_new(1)=xgga_opt
      END IF
      
      tmp_section => section_vals_get_subs_vals(func_section,"CS1",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(2)=tmp_section%section%name
         xc_control%gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
         xc_control%crossterms ( 2 ) = CROSS_TERM
         xc_control%functionals_new(2)=cgga_cs1
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"XGGA",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         CALL section_vals_val_get(tmp_section,"FUNCTIONAL",i_val=i_param,error=error)
        xc_control%gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
        xc_control%functionals_new(1) = i_param
        SELECT CASE(i_param)
        CASE(xgga_b88x)
           xc_control%functionals(1)="BECKE88X"
        CASE(xgga_pw86)
           xc_control%functionals(1)="PW86X"
        CASE(xgga_pw91)
           xc_control%functionals(1)="PW91X"
        CASE(xgga_pbex)
           xc_control%functionals(1)="PBEXX"
        CASE(xgga_revpbe)
           xc_control%functionals(1)="REV_PBEX"
        CASE(xgga_opt)
           xc_control%functionals(1)="OPTX"
        CASE DEFAULT
           CALL stop_program ( 'read_xc_functional', 'Not known XGGA' )
        END SELECT
      END IF
      
      tmp_section => section_vals_get_subs_vals(func_section,"KE_GGA",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         CALL section_vals_val_get(tmp_section,"FUNCTIONAL",i_val=i_param,error=error)
        xc_control%gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
        xc_control%functionals_new(3) = i_param
        SELECT CASE (i_param)
        CASE(ke_ol1)
           xc_control%functionals(3)="KE_OL1"
        CASE(ke_ol2)
           xc_control%functionals(3)="KE_OL2"
        CASE(ke_llp)
           xc_control%functionals(3)="KE_LLP"
        CASE(ke_pw86)
           xc_control%functionals(3)="KE_PW86"
        CASE(ke_pw91)
           xc_control%functionals(3)="KE_PW91"
        CASE(ke_lc)
           xc_control%functionals(3)="KE_LC"
        CASE(ke_t92)
           xc_control%functionals(3)="KE_T92"
        CASE(ke_pbe)
           xc_control%functionals(3)="KE_PBE"
        CASE(ke_revpbe)
           xc_control%functionals(3)="KE_REVPBE"
        CASE DEFAULT
           CALL stop_program ( 'read_xc_functional', 'Not known KE_GGA' )
        END SELECT
      END IF 
      
      tmp_section => section_vals_get_subs_vals(func_section,"P86C",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(2)=tmp_section%section%name
         xc_control%gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
         xc_control%crossterms ( 2 ) = CROSS_TERM
          xc_control%functionals_new(2)=cgga_pw86
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"PW92",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         CALL section_vals_val_get(tmp_section,"PARAMETRIZATION",i_val=i_param,error=error)
         xc_control%functionals_new(2)=i_param
         SELECT CASE (i_param)
         CASE(c_pw92)
            xc_control%functionals(2)= "PW92"
         CASE(c_pw92dmc)
            xc_control%functionals(2)= "PW92_DMC"    
         CASE(c_pw92vmc)
            xc_control%functionals(2)= "PW92_VMC"
         CASE DEFAULT
             CALL stop_program ( 'read_xc_functional', &
             'Not known PW parametrization' )
         END SELECT
         xc_control%crossterms ( 2 ) = CROSS_TERM
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"PZ",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         CALL section_vals_val_get(tmp_section,"PARAMETRIZATION",i_val=i_param,error=error)
         xc_control%functionals_new(2)=i_param
         SELECT CASE (i_param)
         CASE(c_pz)
            xc_control%functionals(2)= "PZ81_ORIG"
         CASE(c_pzdmc)
            xc_control%functionals(2)= "PZ81_DMC"    
         CASE(c_pzvmc)
            xc_control%functionals(2)= "PZ81_VMC"
         CASE DEFAULT
             CALL stop_program ( 'read_xc_functional', &
             'Not known PZ parametrization' )
         END SELECT
         xc_control%crossterms ( 2 ) = CROSS_TERM
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"TFW",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(3)=tmp_section%section%name
         xc_control%gradient_functionals ( 3 ) = GRADIENT_FUNCTIONAL
          xc_control%functionals_new(2)=ke_tfw
      END IF
      
      tmp_section => section_vals_get_subs_vals(func_section,"TF",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(3)=tmp_section%section%name
          xc_control%functionals_new(2)=ke_tfw
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"VWN",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(2)=tmp_section%section%name
         xc_control%crossterms ( 2 ) = CROSS_TERM
          xc_control%functionals_new(2)=c_vwn
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"XALPHA",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(1)="SLATER"
         CALL section_vals_val_get(tmp_section,"XA",r_val=xc_control%xa_param,error=error)
          xc_control%functionals_new(1)=xalpha
      END IF   
      
      tmp_section => section_vals_get_subs_vals(func_section,"TPSS",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(1)="MY_TPSS" !tmp_section%section%name
         xc_control%gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
         xc_control%crossterms ( 1 ) = CROSS_TERM
         xc_control%functionals_new(1)=xgga_tpss
      END IF
      
      tmp_section => section_vals_get_subs_vals(func_section,"PBEX",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(1)=tmp_section%section%name
         xc_control%gradient_functionals ( 1 ) = GRADIENT_FUNCTIONAL
         xc_control%functionals_new(1)=xgga_pbe
      END IF

      tmp_section => section_vals_get_subs_vals(func_section,"PBEC",error=error)
      CALL section_vals_get(tmp_section,explicit=was_present,error=error)
      IF(was_present)THEN
         xc_control%functionals(2)=tmp_section%section%name
         xc_control%gradient_functionals ( 2 ) = GRADIENT_FUNCTIONAL
         xc_control%crossterms ( 2 ) = CROSS_TERM
         xc_control%functionals_new(2)=cgga_pbe
      END IF

    end if

    ! If the functional is not defined from input use defaults
    xc_fun_set=.FALSE.
    DO i=1,SIZE( xc_control%functionals_new )
       IF (xc_control%functionals_new(i) == 0 ) THEN
          xc_control % functionals(i)="NONE"
       ELSE
          xc_fun_set=.TRUE.
       END IF
    END DO
    IF (.NOT.xc_fun_set) THEN
       xc_control % functionals ( 1 ) = "XALPHA"
       xc_control%functionals_new(1) = xalpha
       xc_control%xa_param = 2.0_dp/3.0_dp
       xc_control % functionals ( 2 ) = "PZ81_DMC"
       xc_control%functionals_new(2) = c_pzdmc
       IF (nspins /= 1) THEN
          xc_control %needs%rho_spin=.TRUE.
          xc_control %needs%tau_spin= .FALSE.
       ELSE
          xc_control % needs%rho=.TRUE.
       END IF
    END IF    

    DO i = 1,SIZE( xc_control%functionals_new )
      CALL xc_get_info_yat(functional=xc_control%functionals_new(i),&
           lsd=(nspins/=1),needs=xc_control%needs,&
           max_deriv=xc_control%max_deriv,xparam=xc_control%xa_param,&
           error=error)
    END DO

    ! sets the flag that the kinetic energy needs to be calculated
    ! if some functional needs it
    kenergy_den= kenergy_den .OR.&
         xc_control%needs%tau .OR. xc_control%needs%tau_spin
    ! fills in the old structure (I would like to remove it)
    xc_control%gradient_functionals=&
         xc_control%gradient_functionals .OR. xc_control%needs%norm_drho&
        .OR. xc_control%needs%norm_drho_spin .OR.&
         xc_control%needs%drho .OR. xc_control%needs%drho_spin

    CALL section_vals_get(grid_section,explicit=was_present,error=error)
    if (was_present) then
       CALL section_vals_val_get(grid_section,"XC_SMOOTH_RHO",&
            i_val=xc_control%xc_rho_smooth_id,error=error)
       CALL section_vals_val_get(grid_section,"XC_DERIV",&
            i_val=xc_control%xc_deriv_method_id,error=error)
    end if

  END SUBROUTINE read_xc_control

!***************************************************************************

!!****f* xc_utils/xc_get_info_yat [1.0] *
!!  
!!   NAME
!!     xc_get_info_yat
!!
!!   FUNCTION
!!     get the information about the given functional
!!
!!   NOTES
!!     - yet another time a get_info
!!       I rewrote it to make the transition cleaner
!!       Possibly we should get rid of the functionals names as strings
!!       This get_info should be used only when the functional are
!!       defined from the input and not at every iteration
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
!!     - xparam: input parameter (possibly more than one)  
!!       used to define the functional
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************  
  SUBROUTINE xc_get_info_yat(functional, lsd, reference, shortform,&
             needs, max_deriv, xparam, error)

    INTEGER, INTENT(IN)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    REAL(dp), INTENT(IN)                     :: xparam
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_info_yat', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    SELECT CASE (functional)
      CASE(xgga_b88)
        IF (lsd) THEN
           CALL xb88_lsd_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        ELSE
           CALL xb88_lda_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        END IF
     CASE(cgga_lyp)
        IF (lsd) THEN
           CALL lyp_lsd_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        ELSE
           CALL lyp_lda_info(reference=reference,shortform=shortform,&
                needs=needs, max_deriv=max_deriv, error=error)
        END IF
     CASE(c_pade)
        CALL pade(reference, shortform, lsd=lsd, needs=needs, error=error)
     CASE(cgga_hcth93)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(93,reference,shortform,needs,max_deriv,error=error)
     CASE(cgga_hcth120)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(120,reference,shortform,needs,max_deriv,error=error)
     CASE(cgga_hcth147)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(147,reference,shortform,needs,max_deriv,error=error)
     CASE(cgga_hcth407)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_info(407,reference,shortform,needs,max_deriv,error=error)
     CASE(xgga_opt)
        IF (lsd) THEN
           CALL optx_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE 
           CALL optx_lda_info(reference,shortform,needs,max_deriv,error=error)
        ENDIF
     CASE(cgga_cs1)
        IF (lsd) THEN
           CALL cs1_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL cs1_lda_info(reference, shortform, needs=needs, max_deriv=max_deriv,&
                error=error)
        ENDIF
     CASE(xgga_b88x)
        CALL xgga_info(xgga_b88x,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(xgga_pw86)
        CALL xgga_info(xgga_pw86,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(xgga_pw91)
        CALL xgga_info(xgga_pw91,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(xgga_pbex)
        CALL xgga_info(xgga_pbex,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(xgga_revpbe)
        CALL xgga_info(xgga_revpbe,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(xgga_optx)
        CALL xgga_info(xgga_optx,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_ol1)
        CALL ke_gga_info(ke_ol1,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_ol2)
        CALL ke_gga_info(ke_ol2,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_llp)
        CALL ke_gga_info(ke_llp,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_pw86)
        CALL ke_gga_info(ke_pw86,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_pw91)
        CALL ke_gga_info(ke_pw91,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_pbe)
        CALL ke_gga_info(ke_pbe,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_revpbe)
        CALL ke_gga_info(ke_revpbe,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_lc)
        CALL ke_gga_info(ke_LC,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_t92)
        CALL ke_gga_info(ke_t92,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(cgga_pw86)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        IF(failure) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "BP functional not implemented with LSD")
        END IF
        CALL p86_lda_info(reference,shortform,needs,max_deriv,error=error)
     CASE(c_pw92)
        CALL perdew_wang_info(c_pw92,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(c_pw92dmc)
        CALL perdew_wang_info(c_pw92dmc,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(c_pw92vmc)
        CALL perdew_wang_info(c_pw92vmc,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(c_pz)
        CALL pz_info(c_pz,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(c_pzdmc)
        CALL pz_info(c_pzdmc,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(c_pzvmc)
        CALL pz_info(c_pzvmc,lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(ke_tfw)
        IF (lsd) THEN
           CALL tfw_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL tfw_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE(ke_tf)
        CALL thomas_fermi_info(lsd,reference,shortform,needs,max_deriv,error=error)
     CASE(c_vwn)
        CPPrecondition(.NOT.lsd,cp_failure_level,routineP,error,failure)
        CALL vwn_lda_info(reference,shortform,needs,max_deriv,error=error)
     CASE(xalpha)
        CALL xalpha_info(lsd,reference,shortform,needs,max_deriv,&
             xa_parameter=xparam,error=error)
!     CASE(xgga_tpss)
!        IF (lsd) THEN
!           CALL tpss_lsd_info(reference,shortform,needs,max_deriv,error=error)
!        ELSE
!           CALL tpss_lda_info(reference,shortform,needs,max_deriv,error=error)
!        END IF
     CASE(xgga_tpss)
        IF (lsd) THEN
           CALL my_tpss_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL my_tpss_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE(xgga_pbe)
        IF (PRESENT(reference)) reference=" "
        IF (PRESENT(shortform)) shortform=" "
        IF (PRESENT(max_deriv)) max_deriv=HUGE(0)
     CASE(cgga_pbe)
        IF (lsd) THEN
           CALL pbe_lsd_info(reference,shortform,needs,max_deriv,error=error)
        ELSE
           CALL pbe_lda_info(reference,shortform,needs,max_deriv,error=error)
        END IF
     CASE(0)
        IF (PRESENT(reference)) reference=" "
        IF (PRESENT(shortform)) shortform=" "
        IF (PRESENT(max_deriv)) max_deriv=HUGE(0)
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(cp_to_string(functional))//"' in "//&
             CPSourceFileRef,&
             error=error)
     END SELECT
  
  END SUBROUTINE
  
!!****f* xc_derivatives/xc_eval_yat [1.0] *
!!
!!   NAME
!!     xc_eval_yat
!!
!!   FUNCTION
!!     evaluate a functional (and its derivatives)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - functional: an  that describes which functional should be used
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
     deriv_order, xparam, error)

    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: deriv_order
    REAL(dp), INTENT(IN)                     :: xparam
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_eval_yat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: density_cut, gradient_cut

 !     CALL timeset(routineP,"I"," ",handle)

      failure=.FALSE.

    CALL xc_rho_set_get(rho_set,rho_cutoff=density_cut,&
          drho_cutoff=gradient_cut,error=error)

    SELECT CASE (functional)
      CASE(xgga_b88)
        IF (lsd) THEN
           CALL xb88_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL xb88_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(cgga_lyp)
        IF (lsd) THEN
           CALL lyp_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL lyp_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(c_pade)
        CALL pade(density_cut)
        IF (lsd) THEN
           CALL pade_lsd_pw_eval(deriv_set, rho_set, deriv_order, error=error)
        ELSE
           CALL pade_lda_pw_eval(deriv_set, rho_set, deriv_order, error=error)
        END IF
     CASE(cgga_hcth93)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(93,rho_set,deriv_set,deriv_order,error=error)
     CASE(cgga_hcth120)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(120,rho_set,deriv_set,deriv_order,error=error)
     CASE(cgga_hcth147)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(147,rho_set,deriv_set,deriv_order,error=error)
     CASE(cgga_hcth407)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        CALL hcth_lda_eval(407,rho_set,deriv_set,deriv_order,error=error)
     CASE(xgga_opt)
        IF (lsd) THEN
           CALL optx_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE 
           CALL optx_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        ENDIF
     CASE(cgga_cs1)
        IF (lsd) THEN
           CALL cs1_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL cs1_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        ENDIF
     CASE(xgga_b88x)
        CALL xgga_eval(xgga_b88x,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE(xgga_pw86)
        CALL xgga_eval(xgga_pw86,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE(xgga_pw91)
        CALL xgga_eval(xgga_pw91,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE(xgga_pbex)
        CALL xgga_eval(xgga_pbex,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE(xgga_revpbe)
        CALL xgga_eval(xgga_revpbe,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE(xgga_optx)
        CALL xgga_eval(xgga_opt,lsd,rho_set,deriv_set,deriv_order,error=error)
     CASE(ke_ol1)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_ol1,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_ol1,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_ol2)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_ol2,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_ol2,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_llp)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_llp,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_llp,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_pw86)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_pw86,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_pw86,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_pw91)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_pw91,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_pw91,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_pbe)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_pbe,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_pbe,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_revpbe)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_revpbe,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_revpbe,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_lc)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_LC,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_LC,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_t92)
        IF (lsd) THEN
           CALL ke_gga_lsd_eval(ke_t92,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL ke_gga_lda_eval(ke_t92,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(cgga_pw86)
        CPPrecondition(.NOT.lsd,cp_warning_level,routineP,error,failure)
        IF(failure) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "BP functional not implemented with LSD")
        END IF
        CALL p86_lda_eval(rho_set,deriv_set,deriv_order,error=error)
     CASE(c_pw92)
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval(c_pw92,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval(c_pw92,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(c_pw92dmc)
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval(c_pw92dmc,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval(c_pw92dmc,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(c_pw92vmc)
        IF (lsd) THEN
           CALL perdew_wang_lsd_eval(c_pw92vmc,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL perdew_wang_lda_eval(c_pw92vmc,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(c_pz)
        IF (lsd)  THEN
           CALL pz_lsd_eval(c_pz,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval(c_pz,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(c_pzdmc)
        IF (lsd) THEN 
           CALL pz_lsd_eval(c_pzdmc,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval(c_pzdmc,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(c_pzvmc)
        IF (lsd) THEN
           CALL pz_lsd_eval(c_pzvmc,rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pz_lda_eval(c_pzvmc,rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_tfw)
        IF (lsd) THEN
           CALL tfw_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL tfw_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(ke_tf)
        IF (lsd) THEN
           CALL thomas_fermi_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL thomas_fermi_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF        
     CASE(c_vwn)
        CPPrecondition(.NOT.lsd,cp_failure_level,routineP,error,failure)
        CALL vwn_lda_eval(rho_set,deriv_set,deriv_order,error=error)
     CASE(xalpha)
        IF (lsd) THEN
           CALL xalpha_lsd_eval(rho_set,deriv_set,deriv_order,&
                xa_parameter=xparam,error=error)
        ELSE
           CALL xalpha_lda_eval(rho_set,deriv_set,deriv_order,&
                xa_parameter=xparam,error=error)
        END IF
!     CASE(xgga_tpss)
!        IF (lsd) THEN
!           CALL tpss_lsd_info(reference,shortform,needs,max_deriv,error=error)
!        ELSE
!           CALL tpss_lda_info(reference,shortform,needs,max_deriv,error=error)
!        END IF
     CASE(xgga_tpss)
        IF (lsd) THEN
           CALL my_tpss_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL my_tpss_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(xgga_pbe)
     CASE(cgga_pbe)
        IF (lsd) THEN
           CALL pbe_lsd_eval(rho_set,deriv_set,deriv_order,error=error)
        ELSE
           CALL pbe_lda_eval(rho_set,deriv_set,deriv_order,error=error)
        END IF
     CASE(0)
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
             message="unknown functional '"//TRIM(cp_to_string(functional))//"' in "//&
             CPSourceFileRef,&
             error=error)
     END SELECT
  
  END SUBROUTINE xc_functional_eval
!!*** **********************************************************************
END MODULE xc_utils

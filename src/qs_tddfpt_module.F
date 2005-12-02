!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_tddfpt_module *
!!
!!   NAME
!!     qs_tddfpt_module
!!
!!   FUNCTION
!!     Performs density functional perturbation theory (tddfpt) calculations.
!!     Uses the self consistent approach. The tddfpt calculation uses the ground
!!     state of the unperturbed system as the initial state.
!!
!!****

MODULE qs_tddfpt_module

  USE input_constants,                 ONLY: oe_none,&
                                             oe_lb,&
                                             oe_gllb,&
                                             oe_saop,&
                                             oe_sic, tddfpt_excitations
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_error_get_logger,&
                                             cp_error_type,cp_assert
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,cp_to_string
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_tddfpt_eigensolver,           ONLY: eigensolver
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_type
  USE qs_tddfpt_utils,                 ONLY: find_contributions,&
                                             tddfpt_cleanup,&
                                             tddfpt_init,&
                                             tddfpt_write_banner
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_pot_saop,                     ONLY: add_saop_pot
  USE input_section_types, ONLY:section_vals_type, section_vals_val_get, &
       section_vals_get_subs_vals,  section_vals_retain, section_vals_create, &
       section_vals_set_subs_vals, section_vals_release
  IMPLICIT NONE


  PRIVATE
  PUBLIC :: tddfpt_calculation

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_tddfpt_module'

CONTAINS

! *****************************************************************************

!!****f* cp2k/tddfpt_module/qs_tddfpt_calculation
!!
!!   NAME
!!     tddfpt_calculation
!!
!!   FUNCTION
!!     Performs the perturbation calculation
!!
!!****
  SUBROUTINE tddfpt_calculation(qs_env, glob_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: glob_env
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_calculation', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER                                  :: handle, excitations, oe_corr
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(tddfpt_env_type)                    :: t_env
    TYPE(section_vals_type),POINTER :: input,xc_section,tddfpt_section

    NULLIFY(input)
    CALL get_qs_env(qs_env, input=input,error=error)
    
    CALL section_vals_val_get(input,"dft%excitations",i_val=excitations,error=error)
    IF (.NOT.excitations==tddfpt_excitations) RETURN

    CALL timeset(routineN,"I","",handle)

    tddfpt_section => section_vals_get_subs_vals(input,"DFT%TDDFPT",error=error)
    CALL section_vals_val_get(tddfpt_section,"oe_corr",i_val=oe_corr,error=error)
    IF (oe_corr /= oe_none) THEN
       CALL orbital_eigenvalue_correction(qs_env, error)
    END IF
    
    NULLIFY(p_env)

    CALL tddfpt_write_banner(glob_env)

    !---------------------------------------!
    ! we don't want to update the KS matrix !
    !---------------------------------------!
    qs_env%ks_env%rho_changed = .FALSE.

    CALL tddfpt_init(p_env, t_env, qs_env)

    CALL eigensolver(p_env, qs_env, t_env, error)

    CALL find_contributions(qs_env, t_env, error)

    CALL tddfpt_cleanup(t_env, p_env)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tddfpt_calculation
  
! *****************************************************************************

  !!****f* qs_scf_post/orbital_eigenvalue_correction *
  !!
  !!   NAME
  !!      orbital_eigenvalue_correction
  !!
  !!   FUNCTION
  !!      Apply a special potential to obtain better
  !!      orbital eigenvalues.
  !!
  !!****
  SUBROUTINE orbital_eigenvalue_correction(qs_env, error)
    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL             :: error

    INTEGER                                  :: output_unit,oe_corr
    LOGICAL                                  :: failure, ionode
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(section_vals_type), POINTER :: input, xc_fun_orig, xc_fun_tmp

    CHARACTER(LEN=*), PARAMETER              :: routineN = "orbital_eigenvalue_correction", &
                                                routineP = moduleN//"/"//routineN

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)    
    
    failure = .FALSE.
    NULLIFY(logger, scf_env, input)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    IF(ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    IF (ionode) THEN
       WRITE (output_unit, "(2X,A)") "",&
       "-----------------------------------------------------------------------------",&
       "-                    orbital eigenvalue correction started                  -",&
       "-----------------------------------------------------------------------------",&
       ""    
    END IF

    CALL get_qs_env(qs_env, scf_env=scf_env, input=input, error=error)

    !----------------------!
    ! KS matrix without XC !
    !----------------------!
    xc_fun_orig => section_vals_get_subs_vals(input,"DFT%XC%XC_FUNCTIONAL",&
         error=error)
    CALL section_vals_retain(xc_fun_orig,error=error)
    NULLIFY(xc_fun_tmp)
    CALL section_vals_create(xc_fun_tmp,xc_fun_orig%section,error=error)
    CALL section_vals_set_subs_vals(input,"DFT%XC%XC_FUNCTIONAL",xc_fun_tmp,&
         error=error)
    CALL section_vals_release(xc_fun_tmp,error=error)

    call qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env, qs_env=qs_env, &
                                      ks_matrix=qs_env%matrix_ks, rho=qs_env%rho, &
                                      energy=qs_env%energy,&
                                      calculate_forces=.false.,&
                                      just_energy=.false.,error=error)

    CALL section_vals_set_subs_vals(input,"DFT%XC%XC_FUNCTIONAL",xc_fun_orig,&
         error=error)
    CALL section_vals_release(xc_fun_orig,error=error)

    CALL section_vals_val_get(input,"DFT%TDDFPT%OE_CORR",i_val=oe_corr,error=error)
    IF (oe_corr == oe_saop .OR. &
        oe_corr == oe_lb .OR. &
        oe_corr == oe_gllb) THEN
       CALL add_saop_pot(qs_env%matrix_ks, qs_env, oe_corr, error)
    ELSE IF (oe_corr == oe_sic) THEN
    END IF

  END SUBROUTINE orbital_eigenvalue_correction

END MODULE qs_tddfpt_module

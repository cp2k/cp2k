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

  USE cp_control_utils,                ONLY: read_tddfpt_control
  USE cp_control_types,                ONLY: dft_control_type,&
                                             oe_none,&
                                             oe_saop,&
                                             oe_sic,&
                                             sic_none,&
                                             tddfpt_control_type,&
                                             tddfpt_singlet,&
                                             tddfpt_lanczos,&
                                             xc_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_scf,                          ONLY: scf_env_do_scf
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_tddfpt_eigensolver,           ONLY: eigensolver
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_type
  USE qs_tddfpt_utils,                 ONLY: tddfpt_cleanup,&
                                             tddfpt_init,&
                                             tddfpt_write_banner
  USE scf_control_types,               ONLY: scf_control_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_pot_saop,                     ONLY: add_saop_pot

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
!!   ARGUMENTS
!!     - qs_env      - in    - qs environment with ground state
!!     - glob_env    - in    - global environment
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
 
    INTEGER                                  :: handle
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control
    TYPE(tddfpt_env_type)                    :: t_env

    IF (.NOT.qs_env%dft_control%do_tddfpt_calculation) RETURN

    CALL timeset(routineN,"I","",handle)

    NULLIFY(tddfpt_control)
    CALL get_qs_env(qs_env, tddfpt_control=tddfpt_control)

    !---------------------!
    ! some default values !
    !---------------------!
    tddfpt_control%n_ev = 1
    tddfpt_control%n_reortho = 2
    tddfpt_control%do_kernel = .TRUE.
    tddfpt_control%invert_S = .TRUE.
    tddfpt_control%res_etype = tddfpt_singlet
    tddfpt_control%diag_method = tddfpt_lanczos
    tddfpt_control%oe_corr = oe_none
    tddfpt_control%sic_method_id = sic_none
    tddfpt_control%xc_control = qs_env%dft_control%xc_control

    CALL read_tddfpt_control(tddfpt_control, glob_env)

    IF (tddfpt_control%oe_corr /= oe_none) THEN
       CALL orbital_eigenvalue_correction(qs_env, glob_env, error)
    END IF
    
    NULLIFY(p_env)

    CALL tddfpt_write_banner(glob_env)

    !---------------------------------------!
    ! we don't want to update the KS matrix !
    !---------------------------------------!
    qs_env%ks_env%rho_changed = .FALSE.

    CALL tddfpt_init(p_env, t_env, qs_env, glob_env)

    CALL eigensolver(t_env, p_env, qs_env, glob_env)

    CALL tddfpt_cleanup(t_env, p_env, glob_env)

    CALL timestop(0.0_dbl, handle)

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
  SUBROUTINE orbital_eigenvalue_correction(qs_env, glob_env, error)
    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
         INTENT(IN), TARGET                  :: glob_env
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL             :: error

    INTEGER                                  :: output_unit, i
    LOGICAL                                  :: failure, ionode
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(global_environment_type), POINTER   :: my_glob_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control
    TYPE(xc_control_type), POINTER           :: xc_control

    ! variables to remember from the environments
    INTEGER                                  :: scf_control_max_scf
    TYPE(xc_control_type)                    :: my_xc_control

    CHARACTER(LEN=*), PARAMETER              :: routineN = "orbital_eigenvalue_correction", &
                                                routineP = moduleN//"/"//routineN

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)    
    
    failure = .FALSE.
    NULLIFY(logger, my_glob_env, scf_env, &
            dft_control, scf_control, xc_control, tddfpt_control)
    my_glob_env => glob_env
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    output_unit = cp_logger_get_default_unit_nr(logger)

    IF (ionode) THEN
       WRITE (output_unit, "(2X,A)") "",&
       "-----------------------------------------------------------------------------",&
       "-                    orbital eigenvalue correction started                  -",&
       "-----------------------------------------------------------------------------",&
       ""    
    END IF

    CALL get_qs_env(qs_env, scf_env=scf_env, dft_control=dft_control, &
                    tddfpt_control=tddfpt_control)
    xc_control => dft_control%xc_control

    my_xc_control%functionals = xc_control%functionals

    !----------------------!
    ! KS matrix without XC !
    !----------------------!
    xc_control%functionals = "NONE"
    call qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env, qs_env=qs_env, &
                                      ks_matrix=qs_env%matrix_ks, rho=qs_env%rho, &
                                      energy=qs_env%energy,&
                                      calculate_forces=.false.,&
                                      just_energy=.false.,error=error)

    xc_control%functionals = my_xc_control%functionals

    IF (tddfpt_control%oe_corr == oe_saop) THEN
       CALL add_saop_pot(qs_env%matrix_ks, qs_env, error)
    ELSE IF (tddfpt_control%oe_corr == oe_sic) THEN
    END IF

    IF (tddfpt_control%oe_corr == oe_saop) THEN
    ELSE IF (tddfpt_control%oe_corr == oe_sic) THEN
    END IF

  END SUBROUTINE orbital_eigenvalue_correction

END MODULE qs_tddfpt_module

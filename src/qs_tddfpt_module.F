!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

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
  USE cp_control_types,                ONLY: tddfpt_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_scf,                          ONLY: scf_env_do_scf
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_tddfpt_eigensolver,           ONLY: eigensolver
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_type
  USE qs_tddfpt_utils,                 ONLY: tddfpt_cleanup,&
                                             tddfpt_init,&
                                             tddfpt_write_banner
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

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
!!   SYNOPSIS
!!     Subroutine tddfpt_calculation(qs_env, glob_env)
!!       Implicit None
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (IN):: glob_env
!!     End Subroutine tddfpt_calculation
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
    tddfpt_control%res_etype = "SINGLET"
    tddfpt_control%diag_method = "DAVIDSON"
    tddfpt_control%oe_corr = "NONE"
    tddfpt_control%xc_control = qs_env%dft_control%xc_control

    CALL read_tddfpt_control(tddfpt_control, glob_env)

    IF (tddfpt_control%oe_corr /= "NONE") THEN
       CALL orbital_eigenvalue_correction(qs_env, glob_env, error)
    END IF
    
    NULLIFY(p_env)

    CALL tddfpt_write_banner(glob_env)

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

    TYPE(global_environment_type), POINTER   :: my_glob_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env

    CHARACTER(LEN=*), PARAMETER              :: routineN = "orbital_eigenvalue_correction", &
                                                routineP = moduleN//"/"//routineN

    ! call qs_scf with just 1 scf cycle, the special pot and diagonalization
    
    NULLIFY(my_glob_env)
    my_glob_env => glob_env

    

    CALL get_qs_env(qs_env, scf_env=scf_env)
    CALL scf_env_do_scf(scf_env, qs_env, my_glob_env, error)

  END SUBROUTINE orbital_eigenvalue_correction

END MODULE qs_tddfpt_module

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
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_tddfpt_eigensolver,           ONLY: tddfpt_es
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_type
  USE qs_tddfpt_utils,                 ONLY: tddfpt_cleanup, &
                                             tddfpt_init, &
                                             tddfpt_write_banner
  USE timings,                         ONLY: timeset, timestop

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
  SUBROUTINE tddfpt_calculation(qs_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_environment_type), INTENT(inout)         :: qs_env
    TYPE(global_environment_type), INTENT(in)        :: glob_env

    ! locals
    TYPE(qs_p_env_type), POINTER           :: p_env         ! perturbation environment
    TYPE(tddfpt_env_type)                  :: t_env         ! tddfpt environment

    INTEGER                                :: handle

    CHARACTER(len=*), PARAMETER            :: routineN = 'tddfpt_calculation', &
                                              routineP = moduleN//'/'//routineN

    ! abort if this calculation is not wanted
    IF (.NOT.qs_env%dft_control%do_tddfpt_calculation) RETURN

    CALL timeset(routineN,"I","",handle)

    !-----------------------------------------------------------!
    ! some default values for qs_env%dft_control%tddfpt_control !
    !-----------------------------------------------------------!
    qs_env%dft_control%tddfpt_control%n_ev = 1
    qs_env%dft_control%tddfpt_control%n_reortho = 2
    qs_env%dft_control%tddfpt_control%do_kernel = .TRUE.
    qs_env%dft_control%tddfpt_control%restricted_excitations_type = "SINGLET"
    qs_env%dft_control%tddfpt_control%xc_control = qs_env%dft_control%xc_control

    CALL read_tddfpt_control(qs_env%dft_control%tddfpt_control, glob_env)

    NULLIFY(p_env)

    CALL tddfpt_write_banner(glob_env)

    CALL tddfpt_init(p_env, t_env, qs_env, glob_env)

    CALL tddfpt_es(t_env, p_env, qs_env, glob_env)

    CALL tddfpt_cleanup(t_env, p_env, glob_env)

    CALL timestop(0.0_dbl, handle)

  END SUBROUTINE tddfpt_calculation
  
! *****************************************************************************

END MODULE qs_tddfpt_module

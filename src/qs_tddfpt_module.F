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
!!     state of the unperturbed system as the initial state. So the ground state
!!     density has either to be provided to the subroutine or the be calculated
!!     before the tddfpt calculation can be started. This is why there are two 
!!     subroutines.
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     09-JUL-2002, TCH, development started
!!
!!   SEE ALSO
!!     
!!
!!****

MODULE qs_tddfpt_module

  USE global_types,         ONLY : global_environment_type
  USE kinds,                ONLY : dbl
  USE qs_build_ks_matrix,   ONLY : qs_ks_env_type
  USE qs_environment_types, ONLY : qs_environment_type
  USE qs_p_types,           ONLY : qs_p_env_type
  USE qs_tddfpt_utils
  USE qs_tddfpt_types,      ONLY : tddfpt_env_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: tddfpt_calculation

  CHARACTER(len=*), PARAMETER :: moduleN='qs_tddfpt_module'

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
!!   ARGUMENTSp
!!     - qs_env      - in    - qs environment with ground state
!!     - ks_env      - in    - KS matrix environment
!!     - glob_env    - in    - global environment
!!
!!   MODIFICATION HISTORY
!!     , TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE tddfpt_calculation(qs_env, ks_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_ks_env_type), INTENT(in)          :: ks_env
    TYPE(qs_environment_type), INTENT(in)     :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    TYPE(qs_p_env_type)         :: p_env
    TYPE(tddfpt_env_type)       :: t_env
    CHARACTER(len=*), PARAMETER :: routineN = 'qs_tddfpt_calculation', &
                                   routineP = moduleN//':'//routineN

    ! abort if this calculation is not wanted
    IF (.NOT.qs_env%dft_control%do_tddfpt_calculation) RETURN

    CALL tddfpt_write_banner(glob_env)

    CALL tddfpt_init(p_env, t_env, qs_env, ks_env, glob_env)

    CALL tddfpt_testing(p_env, t_env)

  END SUBROUTINE tddfpt_calculation

! *****************************************************************************

END MODULE qs_tddfpt_module

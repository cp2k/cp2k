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
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE dft_types,                       ONLY: read_tddfpt_control
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_tddfpt_eigensolver
  USE qs_tddfpt_types
  USE qs_tddfpt_utils
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in)        :: glob_env

    ! locals
    TYPE(qs_p_env_type), POINTER           :: p_env         ! perturbation environment
    TYPE(tddfpt_env_type)                  :: t_env         ! tddfpt environment
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: v ! temporary storage

    INTEGER                                :: stat, handle
    INTEGER                                :: n_ev, n_spins ! cache values
    INTEGER                                :: ev, spin      ! loop indices
    LOGICAL                                :: evstat        ! status of the eigenvalues
    REAL(DBL), DIMENSION(:,:), POINTER     :: val_x, val_y  ! 

    CHARACTER(len=*), PARAMETER            :: routineN = 'tddfpt_calculation', &
                                              routineP = moduleN//'/'//routineN

    ! abort if this calculation is not wanted
    IF (.NOT.qs_env%dft_control%do_tddfpt_calculation) RETURN

    CALL timeset(routineN,"I","",handle)

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

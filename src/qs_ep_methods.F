!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_ep_methods [1.0] *
!!
!!   NAME
!!     qs_ep_methods
!!
!!   FUNCTION
!!     in this module there are the routines for the calculation of 
!!     psi and energy without scf, but perturbatively as in
!!     PRL 82,22 26Nov 2001,"Accurate Total Energies without Self-Consistency",
!!     David M. Benoit, Daniel Sebastiani, and Michele Parrinello 
!!
!!   NOTES
!!     At the moment much is hard-coded
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ep_methods
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE qs_environment_types,            ONLY: qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ep_methods'

  PUBLIC :: qs_ep_calculation

!!***
!****************************************************************************
CONTAINS

!!****f* qs_ep_methods/qs_ep_calculation [1.0] *
!!
!!   NAME
!!     qs_ep_calculation
!!
!!   SYNOPSIS
!!     Subroutine qs_ep_calculation(qs_env, globenv, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (INOUT), Target:: globenv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ep_calculation
!!
!!   FUNCTION
!!     Perform an ep calculation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ep_calculation(qs_env,globenv,error)
  TYPE(qs_environment_type), INTENT(inout),TARGET :: qs_env
  TYPE(global_environment_type), INTENT(inout), TARGET :: globenv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ep_calculation',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  ! to do
  
  ! finish pre initialization of qs_env

  ! split up qs_env

  ! energy in subsistem with diag

  ! perturbative approach for complete system

END SUBROUTINE qs_ep_calculation
!***************************************************************************

END MODULE qs_ep_methods

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"

MODULE farming_types
  USE kinds,                           ONLY: default_string_length
  USE cp_log_handling,                 ONLY: cp_logger_type, &
                                             cp_failure_level, &
                                             cp_logger_create, &
                                             cp_logger_release, &
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_error_type, &
                                             cp_a_l, &
                                             cp_assert, &
                                             cp_error_init, &
                                             cp_error_dealloc_ref, &
                                             cp_assertion_failed
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: farming_env_type

  TYPE farming_env_type
     CHARACTER(LEN=default_string_length)                        :: cwd        ! directory we started from
     INTEGER                                                     :: Njobs      ! how many jobs to run
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: Job_cwd    ! list of directories to go to
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: Job_input  ! name of the input file to use
  END TYPE farming_env_type

END MODULE farming_types
!!*****

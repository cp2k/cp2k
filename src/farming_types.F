!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"

MODULE farming_types
  USE kinds,                           ONLY: default_string_length
  USE cp_log_handling,                 ONLY: cp_logger_type, &
                                             cp_failure_level, &
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_error_type, &
                                             cp_a_l, &
                                             cp_assert, &
                                             cp_error_init, &
                                             cp_error_dealloc_ref, &
                                             cp_assertion_failed
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: farming_env_type, deallocate_farming_env

  TYPE farming_env_type
     INTEGER :: group_size_wish      = 0
     LOGICAL :: group_size_wish_set  = .FALSE.
     INTEGER :: ngroup_wish          = 0
     LOGICAL :: ngroup_wish_set      = .FALSE.
     INTEGER, DIMENSION(:), POINTER                              :: group_partition => NULL() ! user preference for partitioning the cpus
     CHARACTER(LEN=default_string_length)                        :: cwd        ="."           ! directory we started from
     INTEGER                                                     :: Njobs      = 0            ! how many jobs to run
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: Job_cwd    => NULL()      ! list of directories to go to
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: Job_input  => NULL()      ! name of the input file to use
     INTEGER, DIMENSION(:), POINTER                              :: Job_group  => NULL()      ! number of the target group
  END TYPE farming_env_type

CONTAINS

!!****f* farming_types/deallocate_farming_env *
!!
!!   NAME
!!     deallocate_farming_env
!!
!!   FUNCTION
!!     deallocates all associated fields of the farming_env type
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2004 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE deallocate_farming_env(farming_env)
     TYPE(farming_env_type), POINTER :: farming_env
     IF (ASSOCIATED(farming_env)) THEN
        IF (ASSOCIATED(farming_env%Job_cwd)) DEALLOCATE(farming_env%Job_cwd)
        IF (ASSOCIATED(farming_env%Job_input)) DEALLOCATE(farming_env%Job_input)
        IF (ASSOCIATED(farming_env%Job_group)) DEALLOCATE(farming_env%Job_group)
        IF (ASSOCIATED(farming_env%group_partition)) DEALLOCATE(farming_env%group_partition)
     ENDIF 
  END SUBROUTINE deallocate_farming_env
END MODULE farming_types
!!*****

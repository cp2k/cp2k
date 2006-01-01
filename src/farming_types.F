!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!

MODULE farming_types
  USE kinds,                           ONLY: default_path_length
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: farming_env_type, deallocate_farming_env, init_farming_env

  TYPE farming_env_type
     INTEGER :: group_size_wish      
     LOGICAL :: group_size_wish_set  
     INTEGER :: ngroup_wish          
     LOGICAL :: ngroup_wish_set      
     LOGICAL :: restart
     LOGICAL :: cycle
     INTEGER, DIMENSION(:), POINTER                              :: group_partition ! user preference for partitioning the cpus
     CHARACTER(LEN=default_path_length)                          :: cwd        ! directory we started from
     INTEGER                                                     :: Njobs      ! how many jobs to run
     INTEGER                                                     :: restart_n  ! where to start
     INTEGER                                                     :: max_steps  ! max number of steps, 
                                                                               ! results in max_steps*Ngroup jobs being run
     CHARACTER(LEN=default_path_length), DIMENSION(:), POINTER :: Job_cwd      ! list of directories to go to
     CHARACTER(LEN=default_path_length), DIMENSION(:), POINTER :: Job_input    ! name of the input file to use
     INTEGER, DIMENSION(:), POINTER                              :: Job_group  ! number of the target group
  END TYPE farming_env_type

CONTAINS

!!****f* farming_types/init_farming_env *
!!
!!   NAME
!!      init_farming_env
!!   SYNOPSIS
!!     Subroutine init_farming_env(farming_env)
!!       Type(farming_env_type), Pointer:: farming_env
!!     End Subroutine init_farming_env
!!
!!   FUNCTION
!!      help poor compilers do their job
!!      i.e. provide a default initialization
!!   NOTES
!!
!!   INPUTS
!!      an associated farming env pointer
!!
!!   MODIFICATION HISTORY
!!     03.2004 created [Joost VandeVondele ]
!!
!!*** **********************************************************************
SUBROUTINE init_farming_env(farming_env)
    TYPE(farming_env_type), POINTER          :: farming_env

    IF (ASSOCIATED(farming_env)) THEN
       farming_env%group_size_wish     = 0 
       farming_env%group_size_wish_set = .FALSE.
       farming_env%ngroup_wish         = 0
       farming_env%ngroup_wish_set     = .FALSE.
       farming_env%restart             = .FALSE.
       farming_env%restart_n           = 1
       farming_env%cycle               = .FALSE.
       NULLIFY(farming_env%group_partition)
       farming_env%cwd                 = "."
       farming_env%Njobs               = 0
          ! so that maxsteps*ngroup is not overflowing
       farming_env%max_steps           = HUGE(0)/32768 
       NULLIFY(farming_env%Job_cwd)
       NULLIFY(farming_env%Job_input)
       NULLIFY(farming_env%Job_group)
    ENDIF
END SUBROUTINE

!!****f* farming_types/deallocate_farming_env *
!!
!!   NAME
!!     deallocate_farming_env
!!
!!   SYNOPSIS
!!     Subroutine deallocate_farming_env(farming_env)
!!       Type(farming_env_type), Pointer:: farming_env
!!     End Subroutine deallocate_farming_env
!!
!!   FUNCTION
!!     deallocates all associated fields of the farming_env type
!!     and the type itself
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
    TYPE(farming_env_type), POINTER          :: farming_env

     IF (ASSOCIATED(farming_env)) THEN
        IF (ASSOCIATED(farming_env%Job_cwd)) DEALLOCATE(farming_env%Job_cwd)
        IF (ASSOCIATED(farming_env%Job_input)) DEALLOCATE(farming_env%Job_input)
        IF (ASSOCIATED(farming_env%Job_group)) DEALLOCATE(farming_env%Job_group)
        IF (ASSOCIATED(farming_env%group_partition)) DEALLOCATE(farming_env%group_partition)
        DEALLOCATE(farming_env) ! and the type itself
     ENDIF 
  END SUBROUTINE deallocate_farming_env
END MODULE farming_types
!!*****

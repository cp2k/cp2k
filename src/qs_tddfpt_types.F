!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_tddfpt_types *
!!
!!   NAME
!!     qs_tddfpt_types
!!
!!   FUNCTION
!!     
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     05-09-2002, TCH, created
!!
!!   SEE ALSO
!!     
!!
!!****
MODULE qs_tddfpt_types

  USE global_types, ONLY : global_environment_type
  USE qs_blacs,     ONLY : blacs_matrix_p_type, replicate_blacs_matrix
  USE qs_p_utils,   ONLY : p_preortho
  USE qs_p_types,   ONLY : qs_p_env_type
  USE termination,  ONLY : stop_memory

  IMPLICIT NONE

  TYPE tddfpt_env_type
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: c_plus
     TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: c_minus
  END TYPE tddfpt_env_type
  
  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_types"

  PRIVATE
  PUBLIC :: tddfpt_env_type, tddfpt_env_init
  
CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/s *
!!
!!   NAME
!!     tddfpt_env_init
!!
!!   FUNCTION
!!     To initialize a tddfpt_env
!!
!!   ARGUMENTS
!!     -
!!
!!   MODIFICATION HISTORY
!!     , TCH, created
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE tddfpt_env_init(t_env, p_env, glob_env)

    IMPLICIT NONE
  
  	! arguments
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(qs_p_env_type), INTENT(in)           :: p_env
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! locals
    INTEGER                            :: stat
    INTEGER                            :: n_spins, spin
    CHARACTER(len=*), PARAMETER        :: routineN = 'tddfpt_env_init', &
                                          routineP = moduleN//'/'//routineN


  	! code
    n_spins = SIZE(p_env%qs_env%k)

    ALLOCATE(t_env%c_plus(n_spins), t_env%c_minus(n_spins), stat=stat)
    IF (stat /= 0) THEN
       CALL stop_memory(routineP, "t_env%c_plus / t_env%c_minus", n_spins)
    END IF
  
    DO spin=1, n_spins
       CALL replicate_blacs_matrix(p_env%psi0d(spin)%blacs_matrix, &
                                   t_env%c_plus(spin)%blacs_matrix, &
                                   name = "C_plus_matrix")
       CALL replicate_blacs_matrix(p_env%psi0d(spin)%blacs_matrix, &
                                   t_env%c_minus(spin)%blacs_matrix, &
                                   name = "C_minus_matrix")
    END DO

    ! now fill the c_plus and c_minus with something meaningfull

    ! finally make them orthogonal to the p_env%psi0d
    CALL p_preortho(p_env, t_env%c_plus, glob_env)
    CALL p_preortho(p_env, t_env%c_minus, glob_env)

  END SUBROUTINE tddfpt_env_init

! *****************************************************************************



END MODULE qs_tddfpt_types

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
!!***
MODULE qs_tddfpt_types
  USE cp_error_handling,               ONLY: cp_simple_assert
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type, &
                                             cp_fm_struct_create
  USE cp_fm_types,                     ONLY: cp_fm_type, &
                                             cp_fm_p_type, &
                                             cp_fm_create, &
                                             cp_fm_release,&
                                             cp_fm_set_all
  USE cp_log_handling,                 ONLY: cp_to_string
  USE cp_control_types,                       ONLY: dft_control_type, &
                                             tddfpt_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

! TYPE DEFINITIONS
! *****************************************************************************

!!****s* qs_tddfpt_types/qs_tddfpt_types
!!
!! NAME
!!   tddfpt_env_type
!!
!! FUNCTION
!!   A type to hold data necessary for a tddfpt calculation
!!
!! ATTRIBUTES
!!   - control : object with control information (see tddfpt_control_type)
!!   - x       : right eigenvector     
!!   - y       : left eigenvector
!!   - eigenvalues : array containing the eigenvalues
!!   - invS    : A full matrix with the inverse of S (inverse of the metric)
!!
!! AUTHOR
!!   T. Chassaing (tchassai)
!!
!!***
  TYPE tddfpt_env_type
     TYPE(tddfpt_control_type)                          :: control
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER     :: X      ! right eigenvector     
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER     :: Y      ! left eigenvector
     REAL(DBL), DIMENSION(:), POINTER                   :: eigenvalues
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: invS   ! the inverse of the metric
  END TYPE tddfpt_env_type
    
  !----------------------------------------------!
  ! The type that will hold the X and Y matrices !
  !----------------------------------------------!
  TYPE fm_p_spins_p_type
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: spins
  END TYPE fm_p_spins_p_type

! *****************************************************************************

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_types"
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .FALSE.
  
  PRIVATE

  ! TYPES
  PUBLIC :: tddfpt_control_type, &
            tddfpt_env_type, &
            fm_p_spins_p_type

  ! METHODS
  PUBLIC :: tddfpt_env_allocate, &
            tddfpt_env_deallocate
  
CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_allocate *
!!
!!   FUNCTION
!!     Allocate enough memory for the tddfpt_environment_type structure
!!
!!***
  SUBROUTINE tddfpt_env_allocate(t_env, p_env, qs_env, glob_env)
  
    IMPLICIT NONE

  ! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  ! LOCALS
    INTEGER                            :: stat          ! 
    INTEGER                            :: n_ev, n_spins ! (convenience values)
    INTEGER                            :: ev, spin      ! loop indices
    TYPE(cp_fm_struct_type), POINTER   :: fm_struct

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_allocate', &
                                   routineP = moduleN//'/'//routineN
  
    ! -------------------------------------------------------------------

    n_ev    = qs_env%dft_control%tddfpt_control%n_ev
    n_spins = qs_env%dft_control%nspins

    ALLOCATE(t_env%eigenvalues(n_ev), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, "t_env%eigenvalues", n_ev)

    !-----------------------------------------!
    ! allocate memory for t_env%x and t_env%y !
    !-----------------------------------------!
    ALLOCATE(t_env%x(n_ev), t_env%y(n_ev), t_env%eigenvalues(n_ev), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, "t_env%x/y", n_ev)
    DO ev=1, n_ev
       ALLOCATE(t_env%x(ev)%spins(n_spins), t_env%y(ev)%spins(n_spins), stat=stat)
       IF (stat/=0) CALL stop_memory(routineP, "t_env : x%spins y%spins", n_spins)
    END DO    
    DO spin=1, n_spins
       DO ev=1, n_ev
          NULLIFY(t_env%x(ev)%spins(spin)%matrix, t_env%y(ev)%spins(spin)%matrix)
          CALL cp_fm_struct_create(fm_struct, qs_env%para_env, qs_env%blacs_env, &
                                   p_env%n_ao(spin), p_env%n_mo(spin))
          CALL cp_fm_create(t_env%x(ev)%spins(spin)%matrix, fm_struct, &
               routineP//"t_env%x-"//TRIM(cp_to_string(ev))//"-"//TRIM(cp_to_string(spin)))
          CALL cp_fm_set_all(t_env%x(ev)%spins(spin)%matrix, 0.0_dbl) ! fill with zeros
       END DO
    END DO
  
    !-------------------------------------------------!
    ! allocate memory for the inverse of the S matrix !
    !-------------------------------------------------!
    ALLOCATE(t_env%invS(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%invS", n_spins)
    DO spin=1, n_spins
       NULLIFY(fm_struct)
       CALL cp_fm_struct_create(fm_struct, qs_env%para_env, qs_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_ao(spin))
                                
       CALL cp_fm_create(t_env%invS(spin)%matrix, fm_struct, routineP//"invS")
    END DO

  END SUBROUTINE tddfpt_env_allocate

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_deallocate *
!!
!!   FUNCTION
!!     To free memory held by a tddftpt_env_type structure
!!
!!***
  SUBROUTINE tddfpt_env_deallocate(t_env, glob_env)
  
    IMPLICIT NONE

  ! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(inout)      :: t_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

  ! LOCALS
    INTEGER                                   :: stat
    INTEGER                                   :: n_ev, n_spins
    INTEGER                                   :: ev, spin

    CHARACTER(len=*), PARAMETER               :: routineN = 'tddfpt_env_deallocate', &
                                                 routineP = moduleN//'/'//routineN
  
    ! -------------------------------------------------------------------

    n_ev    = size(t_env%x)
    n_spins = size(t_env%x(1)%spins)
    
    DO ev=1, n_ev
       DO spin=1, n_spins
          IF (ASSOCIATED(t_env%x(ev)%spins(spin)%matrix)) &
               CALL cp_fm_release(t_env%x(ev)%spins(spin)%matrix)
          IF (ASSOCIATED(t_env%Y(ev)%spins(spin)%matrix)) &
               CALL cp_fm_release(t_env%y(ev)%spins(spin)%matrix)
       END DO
       DEALLOCATE(t_env%x(ev)%spins, t_env%y(ev)%spins, stat=stat)
       IF (stat /= 0) &
            CALL stop_memory(routineP, "t_env%x%spins, t_env%y%spins")
    END DO
    DEALLOCATE(t_env%x, t_env%y, stat=stat)
    IF (stat /= 0) &
         CALL stop_memory(routineP, "t_env%x, t_env%y")
    NULLIFY(t_env%x, t_env%y)

    DEALLOCATE(t_env%eigenvalues, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%ev")

    DO spin=1, n_spins
       IF (ASSOCIATED(t_env%invS(spin)%matrix)) &
            CALL cp_fm_release(t_env%invS(spin)%matrix)
    END DO
    DEALLOCATE (t_env%invS, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%invS")
    NULLIFY(t_env%invS)
 
  END SUBROUTINE tddfpt_env_deallocate

! *****************************************************************************

END MODULE qs_tddfpt_types

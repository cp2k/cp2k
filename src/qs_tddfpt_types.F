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
  USE cp_control_types,                ONLY: tddfpt_control_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all
  USE cp_log_handling,                 ONLY: cp_to_string
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE hartree_local_types,             ONLY: hartree_local_create,&
                                             hartree_local_release,&
                                             hartree_local_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

! *****************************************************************************

!!****s* qs_tddfpt_types/qs_tddfpt_types
!!
!! NAME
!!   tddfpt_env_type
!!
!! FUNCTION
!!   A type to hold data necessary for a tddfpt calculation
!!
!!***
  TYPE tddfpt_env_type
     TYPE(tddfpt_control_type)                          :: control
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER     :: X      ! right eigenvector     
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER     :: Y      ! left eigenvector
     REAL(DBL), DIMENSION(:), POINTER                   :: eigenvalues
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: invS   ! the inverse of the metric
  END TYPE tddfpt_env_type
    
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
  SUBROUTINE tddfpt_env_allocate(t_env, p_env, qs_env)
  

    TYPE(tddfpt_env_type), INTENT(out)       :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_allocate', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: ev, n_ev, n_spins, spin, stat
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct

    n_ev    = qs_env%dft_control%tddfpt_control%n_ev
    n_spins = qs_env%dft_control%nspins

    ALLOCATE(t_env%eigenvalues(n_ev), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, "t_env%eigenvalues", n_ev)

    !-----------------------------------------!
    ! allocate memory for t_env%x and t_env%y !
    !-----------------------------------------!
    ALLOCATE(t_env%x(n_ev), t_env%y(n_ev), stat=stat)
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
          CALL cp_fm_struct_release(fm_struct)
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
       CALL cp_fm_struct_release(fm_struct)
    END DO

  END SUBROUTINE tddfpt_env_allocate

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_deallocate *
!!
!!   FUNCTION
!!     To free memory held by a tddftpt_env_type structure
!!
!!***
  SUBROUTINE tddfpt_env_deallocate(t_env)
  

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_deallocate', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: ev, n_ev, n_spins, spin, stat

! -------------------------------------------------------------------

    n_ev    = size(t_env%x)
    n_spins = size(t_env%x(1)%spins)
    
    DO ev=1, n_ev
       DO spin=1, n_spins
          IF (ASSOCIATED(t_env%x(ev)%spins(spin)%matrix)) &
               CALL cp_fm_release(t_env%x(ev)%spins(spin)%matrix)
          IF (ASSOCIATED(t_env%y(ev)%spins(spin)%matrix)) &
               CALL cp_fm_release(t_env%y(ev)%spins(spin)%matrix)
       END DO
       DEALLOCATE(t_env%x(ev)%spins, t_env%y(ev)%spins, stat=stat)
       IF (stat /= 0) &
            CALL stop_memory(routineP, "t_env%x%spins, t_env%y%spins")
    END DO
    DEALLOCATE(t_env%x, t_env%y, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%x, t_env%y")

    DEALLOCATE(t_env%eigenvalues, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%ev")

    DO spin=1, n_spins
       IF (ASSOCIATED(t_env%invS(spin)%matrix)) &
            CALL cp_fm_release(t_env%invS(spin)%matrix)
    END DO
    DEALLOCATE (t_env%invS, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%invS")

    NULLIFY(t_env%x, t_env%y, t_env%invS)
  
  END SUBROUTINE tddfpt_env_deallocate

! *****************************************************************************

END MODULE qs_tddfpt_types

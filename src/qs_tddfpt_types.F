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
!!****
MODULE qs_tddfpt_types

  USE cp_error_handling, ONLY : cp_simple_assert
  USE cp_fm_types,    ONLY : cp_fm_set_all, cp_fm_release
  USE cp_log_handling,   ONLY : cp_to_string
  USE global_types,      ONLY : global_environment_type
  USE kinds,             ONLY : dbl, dp_size
  USE qs_blacs,          ONLY : cp_full_matrix_p_type, &
                                cp_fm_vect_dealloc, &
                                replicate_blacs_matrix
  USE qs_p_utils,        ONLY : p_preortho
  USE qs_p_types,        ONLY : qs_p_env_type
  USE termination,       ONLY : stop_memory

  IMPLICIT NONE

  TYPE real_2d_array_p_type
     REAL(dbl), DIMENSION(:,:), POINTER :: values
  END TYPE real_2d_array_p_type

  TYPE fm_p_spins_p_type
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: spins
  END TYPE fm_p_spins_p_type

!!****s* qs_tddfpt_types
!!
!! NAME
!!   tddfpt_env_type
!!
!! FUNCTION
!!   A type to hold data necessary for a tddfpt calculation
!!
!! ATTRIBUTES
!!   - n_ev   : number of eigenvalues 
!!   - x      : right eigenvector     
!!   - y      : left eigenvector      
!!   - eval_x : the eigenvalues       
!!   - eval_y : the eigenvalues       
!!
!! NOTES
!!
!!****
  TYPE tddfpt_env_type
     INTEGER                                           :: n_ev   ! number of eigenvalues 
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER    :: x      ! right eigenvector     
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER    :: y      ! left eigenvector
     TYPE(real_2d_array_p_type), DIMENSION(:), POINTER :: eval_x ! the eigenvalues       
     TYPE(real_2d_array_p_type), DIMENSION(:), POINTER :: eval_y ! the eigenvalues       
  END TYPE tddfpt_env_type
    
  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_types"
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.
  
  PRIVATE
  PUBLIC :: real_2d_array_p_type, fm_p_spins_p_type
  PUBLIC :: tddfpt_env_allocate, tddfpt_env_deallocate, &
            tddfpt_env_type, tddfpt_env_init
  
CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_init *
!!
!!   NAME
!!     tddfpt_env_init
!!
!!   FUNCTION
!!     To initialize a tddfpt_env_type which is to give it
!!     meaningfull initial values
!!
!!   ARGUMENTS
!!     - t_env
!!     - p_env
!!     - glob_env
!!
!!   NOTES
!!     The p_env variable must be fully initialized, since it is used
!!     to determine quantities necessary for the initialization.
!!
!!   MODIFICATION HISTORY
!!     08-2002, TCH, created
!!
!!****
  SUBROUTINE tddfpt_env_init(t_env, p_env, n_ev, glob_env)

    IMPLICIT NONE
  
  	! arguments
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(qs_p_env_type), INTENT(in)           :: p_env
    INTEGER, INTENT(IN)                       :: n_ev
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! locals
    INTEGER                            :: stat
    INTEGER                            :: ev, n_spins, spin
    CHARACTER(len=*), PARAMETER        :: routineN = 'tddfpt_env_init', &
                                          routineP = moduleN//'/'//routineN


  	! code
    t_env%n_ev = n_ev

    ! fill the c_plus and c_minus (x and y) with something meaningfull
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, '(1X,A)') routineP//": Initialization of x and y not yet implemented"
    END IF

    ! finally make them orthogonal to the p_env%psi0d
    DO ev=1, n_ev
       CALL p_preortho(p_env, t_env%x(ev)%spins)
       CALL p_preortho(p_env, t_env%y(ev)%spins)
    END DO

  END SUBROUTINE tddfpt_env_init

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_allocate *
!!
!!   NAME
!!     tddfpt_env_allocate
!!
!!   FUNCTION
!!     Allocate the memory of the fields in an tddftp_env_type.
!!
!!   ARGUMENTS
!!     - t_env : the environment to allocate memory for (inout)
!!     - p_env : the perturbation environment (in)
!!     - glob_env : global environment (in)
!!
!!   NOTES
!!     The p_env variable must be fully initialized, since it is used
!!     to determine quantities necessary for the allocation.
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, created
!!
!!****
  SUBROUTINE tddfpt_env_allocate(t_env, p_env, n_ev, glob_env)
  
    IMPLICIT NONE

  	! arguments
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(qs_p_env_type), INTENT(in)           :: p_env
    INTEGER, INTENT(IN)                       :: n_ev
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! locals
    INTEGER                            :: stat
    INTEGER                            :: ev, n_spins, spin
    CHARACTER(len=*), PARAMETER        :: routineN = 'tddfpt_env_allocate', &
                                          routineP = moduleN//'/'//routineN
  
  
  	! code
    n_spins = p_env%qs_env%dft_control%nspins

    ALLOCATE(t_env%x(n_ev), t_env%y(n_ev), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, "t_env%x, t_env%y", n_ev)

    ALLOCATE(t_env%eval_x(n_spins), t_env%eval_y(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%eval_x, t_env%eval_y", n_spins)

    DO ev=1, n_ev
       ALLOCATE(t_env%x(ev)%spins(n_spins), t_env%y(ev)%spins(n_spins), stat=stat)
       IF (stat/=0) CALL stop_memory(routineP, "t_env%x%spins, t_env%y%spins", n_spins)
    END DO
    
    DO spin=1, n_spins

       ! arpack need one field more, hence the +1
       ALLOCATE(t_env%eval_x(spin)%values(n_ev+1, 2), &
                t_env%eval_y(spin)%values(n_ev+1, 2), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, "t_env%eval_x%values, t_env%eval_y%values", &
                                       2*n_ev*dp_size)
       t_env%eval_x(spin)%values = 0.0_dbl
       t_env%eval_y(spin)%values = 0.0_dbl

       DO ev=1, n_ev
          CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, &
               t_env%x(ev)%spins(spin)%matrix, &
               name = "X ev = "//cp_to_string(ev)//" spin = "//cp_to_string(spin))
          CALL cp_fm_set_all(t_env%x(ev)%spins(spin)%matrix, 0.0_dbl) ! fill with zeros
          CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, &
               t_env%y(ev)%spins(spin)%matrix, &
               name = "Y ev = "//cp_to_string(ev)//" spin = "//cp_to_string(spin))
          CALL cp_fm_set_all(t_env%y(ev)%spins(spin)%matrix, 0.0_dbl) ! fill with zeros
       END DO
    END DO
  
  END SUBROUTINE tddfpt_env_allocate

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_deallocate *
!!
!!   NAME
!!     tddfpt_env_deallocate
!!
!!   FUNCTION
!!     To free memory held by a tddftpt_env_type.
!!
!!   ARGUMENTS
!!     - t_env
!!
!!   NOTES
!!     -
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, created
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE tddfpt_env_deallocate(t_env, glob_env)
  
    IMPLICIT NONE

  	! arguments
    TYPE(tddfpt_env_type), INTENT(inout)      :: t_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

  	! locals
    INTEGER                                   :: ev, n_spins, spin, stat
    CHARACTER(len=*), PARAMETER               :: routineN = 'tddfpt_env_deallocate', &
                                                 routineP = moduleN//'/'//routineN
  
  	! code
    n_spins = SIZE(t_env%x(1)%spins)

    CPAAssert(t_env%n_ev==SIZE(t_env%x))
    CPAAssert(SIZE(t_env%y(1)%spins)==n_spins)
    CPAAssert(SIZE(t_env%eval_x)==SIZE(t_env%x(1)%spins))
    CPAAssert(SIZE(t_env%eval_y)==SIZE(t_env%y(1)%spins))

    DO spin=1, n_spins
       DEALLOCATE(t_env%eval_x(spin)%values, t_env%eval_y(spin)%values, stat=stat)
       IF (stat /= 0) &
            CALL stop_memory(routineP, "t_env%eval_x%values, t_env%eval_y%values")
       NULLIFY(t_env%eval_x(spin)%values)
       NULLIFY(t_env%eval_y(spin)%values)       
    END DO
    DEALLOCATE(t_env%eval_x, t_env%eval_y, stat=stat)
    IF (stat /= 0) &
         CALL stop_memory(routineP, "t_env%eval_x, t_env%eval_y")
    NULLIFY(t_env%eval_x); NULLIFY(t_env%eval_y)
    NULLIFY(t_env%x); NULLIFY(t_env%y)

    DO ev=1, t_env%n_ev
       DO spin=1, n_spins
          CALL cp_fm_release(t_env%x(ev)%spins(spin)%matrix)
          NULLIFY(t_env%x(ev)%spins(spin)%matrix)
          CALL cp_fm_release(t_env%y(ev)%spins(spin)%matrix)
          NULLIFY(t_env%y(ev)%spins(spin)%matrix)
       END DO
       DEALLOCATE(t_env%x(ev)%spins, t_env%y(ev)%spins, stat=stat)
       IF (stat /= 0) &
            CALL stop_memory(routineP, "t_env%x%spins, t_env%y%spins")
       NULLIFY(t_env%x(ev)%spins)
       NULLIFY(t_env%y(ev)%spins)
    END DO
    DEALLOCATE(t_env%x, t_env%y, stat=stat)
    IF (stat /= 0) &
         CALL stop_memory(routineP, "t_env%x, t_env%y")
    NULLIFY(t_env%x)
    NULLIFY(t_env%y)

  END SUBROUTINE tddfpt_env_deallocate

! *****************************************************************************

END MODULE qs_tddfpt_types

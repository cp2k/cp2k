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
  USE cp_fm_types,                     ONLY: cp_fm_release,&
                                             cp_fm_set_all
  USE cp_log_handling,                 ONLY: cp_to_string
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp_size
  USE qs_blacs,                        ONLY: cp_fm_vect_dealloc,&
                                             cp_full_matrix_p_type,&
                                             replicate_blacs_matrix
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_types,                      ONLY: qs_p_env_type
  USE qs_p_utils,                      ONLY: p_preortho
  USE termination,                     ONLY: stop_memory, stop_program
  IMPLICIT NONE

  TYPE real_2d_array_p_type
     REAL(dbl), DIMENSION(:,:), POINTER :: values
  END TYPE real_2d_array_p_type

  TYPE fm_p_spins_p_type
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: spins
  END TYPE fm_p_spins_p_type

  TYPE tddfpt_control_type
     INTEGER :: n_ev    ! number of eigenvalues
     INTEGER :: n_spins ! number of spins used in the calculation (cached value)
  END TYPE tddfpt_control_type

!!****s* qs_tddfpt_types/qs_tddfpt_types
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
!!***
  TYPE tddfpt_env_type
     TYPE(tddfpt_control_type), POINTER                :: control
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER    :: x      ! right eigenvector     
     TYPE(fm_p_spins_p_type), DIMENSION(:), POINTER    :: y      ! left eigenvector
     TYPE(real_2d_array_p_type), DIMENSION(:), POINTER :: eval_x ! the eigenvalues       
     TYPE(real_2d_array_p_type), DIMENSION(:), POINTER :: eval_y ! the eigenvalues       
  END TYPE tddfpt_env_type
    
  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_types"
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .FALSE.
  
  PRIVATE

  ! TYPES
  PUBLIC :: tddfpt_control_type, tddfpt_env_type
  PUBLIC :: real_2d_array_p_type, fm_p_spins_p_type ! must be public so other modules
                                                    ! can access their fields

  ! METHODS
  PUBLIC :: tddfpt_env_allocate, tddfpt_env_deallocate, &
            tddfpt_env_init, check_tddfpt_env
  
CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_allocate *
!!
!!   FUNCTION
!!     Allocate enough memory for the tddfpt_environment_type structure
!!
!!***
  SUBROUTINE tddfpt_env_allocate(t_env, t_control, p_env, qs_env, glob_env)
  
    IMPLICIT NONE

  	! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(tddfpt_control_type), POINTER        :: t_control
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! LOCALS
    INTEGER                     :: stat          ! 
    INTEGER                     :: n_ev, n_spins ! (convenience values)
    INTEGER                     :: ev, spin      ! loop indices

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_allocate', &
                                   routineP = moduleN//'/'//routineN
  
    ! -------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(t_control)) THEN
       CALL stop_program(routineP, "t_control is not associated")
    END IF

    n_ev    = t_control%n_ev
    n_spins = t_control%n_spins

    ALLOCATE(t_env%x(n_ev), t_env%y(n_ev), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, "t_env%x, t_env%y", n_ev)

    ALLOCATE(t_env%eval_x(n_spins), t_env%eval_y(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "t_env%eval_x, t_env%eval_y", n_spins)

    DO ev=1, n_ev
       ALLOCATE(t_env%x(ev)%spins(n_spins), t_env%y(ev)%spins(n_spins), stat=stat)
       IF (stat/=0) CALL stop_memory(routineP, "t_env%x%spins, t_env%y%spins", n_spins)
    END DO
    
    DO spin=1, n_spins

       ! arpack needs one field more, hence the +1
       ! the 2 is because the eigenvalues can, in principle,  be complex
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

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_init *
!!
!!   FUNCTION
!!     Initialize t_env with meaningfull values.
!!
!!***
  SUBROUTINE tddfpt_env_init(t_env, t_control, p_env, qs_env, glob_env)

    IMPLICIT NONE
  
  	! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(tddfpt_control_type), POINTER        :: t_control
    TYPE(qs_p_env_type),POINTER               :: p_env
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! LOCALS
    INTEGER                            :: stat
    INTEGER                            :: ev 

    CHARACTER(len=*), PARAMETER        :: routineN = 'tddfpt_env_init', &
                                          routineP = moduleN//'/'//routineN

    ! -------------------------------------------------------------------

    IF (ASSOCIATED(t_control)) THEN
       t_env%control => t_control
    ELSE
       CALL stop_program(routineP, "t_control not associated")
    END IF

    !----------------------------!
    ! fill x and y of the t_env  !
    ! with something meaningfull !
    !----------------------------!
    IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A)') &
         routineP//": Initialization of x and y not yet implemented"
    
    ! finally make them orthogonal to the p_env%psi0d
    DO ev=1, t_env%control%n_ev
       CALL p_preortho(p_env, qs_env, t_env%x(ev)%spins)
       CALL p_preortho(p_env, qs_env, t_env%y(ev)%spins)
    END DO

  END SUBROUTINE tddfpt_env_init

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

    CPAAssert(check_tddfpt_env(t_env))

    n_ev    = t_env%control%n_ev
    n_spins = t_env%control%n_spins
    
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

    DO ev=1, n_ev
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
    NULLIFY(t_env%x); NULLIFY(t_env%y)

    IF (ASSOCIATED(t_env%control)) THEN
       DEALLOCATE(t_env%control, stat=stat)
       IF (stat/=0) CALL stop_memory(routineP, "t_env%control")
    END IF
 
  END SUBROUTINE tddfpt_env_deallocate

! *****************************************************************************

  FUNCTION check_tddfpt_env(t_env) RESULT (valid)
    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(in) :: t_env
    LOGICAL                           :: valid

    ! LOCALS
    INTEGER                           :: n_ev, n_spins
    INTEGER                           :: ev, spin

    CHARACTER(LEN=*), PARAMETER       :: routineN = 'check_tddfpt_env', &
                                         routineP = moduleN//'/'//routineN

    ! -------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(t_env%control)) &
         CALL stop_program(routineP, "t_env%control is not associated")
    n_ev    = t_env%control%n_ev
    n_spins = t_env%control%n_spins

    valid = .TRUE.
    IF (.NOT.(n_ev>=1)) valid = .FALSE.
    IF (.NOT.(ASSOCIATED(t_env%x))) valid = .FALSE.
    IF (.NOT.(n_ev==SIZE(t_env%x))) valid = .FALSE.
    IF (.NOT.(ASSOCIATED(t_env%y))) valid = .FALSE.
    IF (.NOT.(n_ev==SIZE(t_env%y))) valid = .FALSE.

    DO ev=1, n_ev
       IF (.NOT.(ASSOCIATED(t_env%x(ev)%spins))) valid = .FALSE.
       IF (.NOT.(SIZE(t_env%x(ev)%spins)==n_spins)) valid = .FALSE.
       IF (.NOT.(ASSOCIATED(t_env%y(ev)%spins))) valid = .FALSE.
       IF (.NOT.(SIZE(t_env%y(ev)%spins)==n_spins)) valid = .FALSE.
    END DO

    IF (.NOT.(ASSOCIATED(t_env%eval_x))) valid = .FALSE.
    IF (.NOT.(SIZE(t_env%eval_x)==n_spins)) valid = .FALSE.
    IF (.NOT.(ASSOCIATED(t_env%eval_y))) valid = .FALSE.
    IF (.NOT.(SIZE(t_env%eval_y)==n_spins)) valid = .FALSE.

    DO spin=1, n_spins
       IF (.NOT.(ASSOCIATED(t_env%eval_x(spin)%values))) valid = .FALSE.
       IF (.NOT.(ASSOCIATED(t_env%eval_y(spin)%values))) valid = .FALSE.
    END DO

  END FUNCTION check_tddfpt_env

! *****************************************************************************

END MODULE qs_tddfpt_types

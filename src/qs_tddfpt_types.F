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

  USE global_types,   ONLY : global_environment_type
  USE kinds,          ONLY : dbl
  USE qs_blacs,       ONLY : cp_full_matrix_p_type, &
                             cp_fm_vect_dealloc, &
                             replicate_blacs_matrix
  USE qs_p_utils,     ONLY : p_preortho
  USE qs_p_types,     ONLY : qs_p_env_type
  USE termination,    ONLY : stop_memory

  IMPLICIT NONE

  TYPE real_2d_array_p_type
     REAL(dbl), DIMENSION(:,:), POINTER :: values
  END TYPE real_2d_array_p_type

  TYPE tddfpt_env_type
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: x    ! right eigenvector
     TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: y    ! left eigenvector
     TYPE(real_2d_array_p_type), DIMENSION(:), POINTER  :: eval_x ! the eigenvalues
     TYPE(real_2d_array_p_type), DIMENSION(:), POINTER  :: eval_y ! the eigenvalues
  END TYPE tddfpt_env_type
    
  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_types"

  PRIVATE
  PUBLIC :: tddfpt_env_type, tddfpt_env_init, tddfpt_env_deallocate, real_2d_array_p_type
  
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
    CALL tddfpt_env_allocate(t_env, p_env, glob_env)

    ! fill the c_plus and c_minus (x and y) with something meaningfull
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, '(A)') routineP//": Initialization of x and y not yet implemented"
    END IF

    ! finally make them orthogonal to the p_env%psi0d
    CALL p_preortho(p_env, t_env%x, glob_env)
    CALL p_preortho(p_env, t_env%y, glob_env)

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
  SUBROUTINE tddfpt_env_allocate(t_env, p_env, glob_env)
  
    IMPLICIT NONE

  	! arguments
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(qs_p_env_type), INTENT(in)           :: p_env
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! locals
    INTEGER                            :: stat
    INTEGER                            :: n_spins, spin
    CHARACTER(len=*), PARAMETER        :: routineN = 'tddfpt_env_allocate', &
                                          routineP = moduleN//'/'//routineN
  
  
  	! code
    n_spins = p_env%qs_env%dft_control%nspins

    ALLOCATE(t_env%x(n_spins), t_env%y(n_spins), &
             t_env%eval_x(n_spins), t_env%eval_y(n_spins), &
             stat=stat)
    IF (stat /= 0) THEN
       CALL stop_memory(routineP, "t_env%x / t_env%y / t_env%eval_x / t_env%eal_y", n_spins)
    END IF
    
    DO spin=1, n_spins
       
       ALLOCATE(t_env%eval_x(spin)%values(p_env%n_mo(spin), 2), &
                t_env%eval_y(spin)%values(p_env%n_mo(spin), 2), stat=stat)
       IF (stat /= 0) THEN
          CALL stop_memory(routineP, "t_env%eval(spin)%values", 2*p_env%n_mo(spin))
       END IF
       t_env%eval_x(spin)%values = 0.0_dbl
       t_env%eval_y(spin)%values = 0.0_dbl

       CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, &
                                   t_env%x(spin)%matrix, &
                                   name = "X_matrix")
       CALL replicate_blacs_matrix(p_env%psi0d(spin)%matrix, &
                                   t_env%y(spin)%matrix, &
                                   name = "Y_matrix")
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
    INTEGER                                   :: n_spins, spin, stat
    CHARACTER(len=*), PARAMETER               :: routineN = 'tddfpt_env_deallocate', &
                                                 routineP = moduleN//'/'//routineN
  
  	! code
    CALL cp_fm_vect_dealloc(t_env%x)
    CALL cp_fm_vect_dealloc(t_env%y)
    
    DO spin=1, SIZE(t_env%eval_x)
       DEALLOCATE(t_env%eval_x(spin)%values, t_env%eval_y(spin)%values, stat=stat)
       IF (stat /= 0) THEN
          IF (glob_env%ionode) WRITE (glob_env%scr, '(A)') &
               routineP//": deallocation of t_env%eval_x%values / t_env%eval_y%values failed"
       END IF
       NULLIFY(t_env%eval_x(spin)%values)
       NULLIFY(t_env%eval_y(spin)%values)       
    END DO
    DEALLOCATE(t_env%eval_x, t_env%eval_y, stat=stat)
    IF (stat /= 0) THEN
       IF (glob_env%ionode) WRITE (glob_env%scr, '(A)') &
            routineP//": deallocation of t_env%eval_x / t_env%eval_y failed"
    END IF
    NULLIFY(t_env%eval_x)
    NULLIFY(t_env%eval_y)

  END SUBROUTINE tddfpt_env_deallocate

! *****************************************************************************

END MODULE qs_tddfpt_types

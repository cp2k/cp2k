!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_tddfpt_utils *
!!
!!   NAME
!!     qs_tddfpt_utils
!!
!!   FUNCTION
!!
!!   SYNOPSIS
!!     
!!   MODIFICATION HISTORY
!!     09-JUL-2002, TCH, development started
!!
!!   SEE ALSO
!!     
!!****

MODULE qs_tddfpt_utils

  USE cp_b_matrix_structure, ONLY : cp_b_matrix_struct_type
  USE cp_block_matrix,       ONLY : cp_block_matrix_type, &
                                    cp_sp_create
  USE cp_error_handling,     ONLY : cp_error_type, &
                                    cp_error_init
  USE cp_log_handling,       ONLY : cp_failure_level
  USE global_types,          ONLY : global_environment_type
  USE kinds,                 ONLY : dbl
  USE qs_blacs,              ONLY : blacs_matrix_type, blacs_matrix_p_type
  USE qs_build_ks_matrix,    ONLY : qs_ks_env_type
  USE qs_environment_types,  ONLY : qs_environment_type, get_qs_env
  USE qs_mo_types,           ONLY : get_mo_set
  USE qs_p_build_kernel,     ONLY : kpp1_init
  USE qs_p_types,            ONLY : p_env_init,&
                                    qs_p_env_type, &
                                    qs_K_P_P1_env_type
  USE qs_tddfpt_types,       ONLY : tddfpt_env_type, tddfpt_env_init

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: tddfpt_testing, &
       tddfpt_init, &
       tddfpt_write_banner

CONTAINS

! *****************************************************************************

  SUBROUTINE tddfpt_init(p_env, t_env, qs_env, ks_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(out)          :: p_env
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env 
    TYPE(qs_environment_type), INTENT(in)     :: qs_env
    TYPE(qs_ks_env_type), INTENT(in)          :: ks_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    TYPE(qs_K_P_P1_env_type)                  :: kpp1_env
    TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: psi0d
    TYPE(blacs_matrix_type), POINTER          :: eigenvectors
    INTEGER                                   :: n_spins, spin
    
    ! do intializations
    CALL kpp1_init(kpp1_env, qs_env, ks_env, glob_env)

    IF (glob_env%ionode) WRITE(glob_env%scr, *) &
         "returned from kpp1_init"

    n_spins = SIZE(qs_env%k)
    ALLOCATE(psi0d(n_spins))

    DO spin=1, n_spins
       CALL get_mo_set(qs_env%C(spin)%mo_set, eigenvectors=eigenvectors)
       psi0d(spin)%blacs_matrix => eigenvectors
    END DO

    IF (glob_env%ionode) WRITE(glob_env%scr, *) &
         "finished getting eigenvectors"

    CALL p_env_init(p_env, qs_env, ks_env, kpp1_env, glob_env, &
                    psi0d, .TRUE.)
    IF (glob_env%ionode) WRITE(glob_env%scr, *) &
         "returned from p_env_init"

    CALL tddfpt_env_init(t_env, p_env, glob_env)
    IF (glob_env%ionode) WRITE(glob_env%scr, *) &
         "returned from tddfpt_env_init"
    
  END SUBROUTINE tddfpt_init

! *****************************************************************************

  SUBROUTINE tddfpt_write_banner(globenv)

    IMPLICIT NONE
    
    ! arguments
    TYPE(global_environment_type), INTENT(in) :: globenv
    
    IF (.NOT.globenv%ionode) RETURN
    IF (globenv%print%program_banner) THEN
       
       WRITE (unit=globenv%scr, fmt="(T2,A79)")&
            "",&
            "-------------------------------------------------------------------------------",&
            "-                     Time-dependent calculation started                      -",&
            "-------------------------------------------------------------------------------",&
            ""    
    END IF
  END SUBROUTINE tddfpt_write_banner

! *****************************************************************************

  SUBROUTINE tddfpt_testing(p_env, t_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(in)      :: p_env
    TYPE(tddfpt_env_type), INTENT(inout) :: t_env

    ! locals
    
  END SUBROUTINE tddfpt_testing

! *****************************************************************************

END MODULE qs_tddfpt_utils

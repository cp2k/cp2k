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
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_block_matrix,                 ONLY: cp_block_matrix_type,&
                                             cp_sp_create
  USE cp_error_handling,               ONLY: cp_error_init,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_blacs,                        ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             replicate_blacs_matrix
  USE qs_build_ks_matrix,              ONLY: qs_ks_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_build_kernel,               ONLY: kpp1_create
  USE qs_p_types,                      ONLY: p_env_init,&
                                             qs_K_P_P1_env_type,&
                                             qs_p_env_type
  USE qs_p_utils,                      ONLY: p_op_l1
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_allocate,&
                                             tddfpt_env_init,&
                                             tddfpt_env_type
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: moduleN="qs_tddfpt_utils"
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  PRIVATE
  PUBLIC :: tddfpt_init, &
            tddfpt_write_banner, &
            update_p_env

CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_utils/tddfpt_init *
!!
!!   NAME
!!     tddfpt_init
!!
!!   FUNCTION
!!     Initialize some necessary structures for a tddfpt calculation.
!!
!!   ARGUMENTS
!!     - p_env    : perturbation environment to be initialized
!!     - t_env    : tddfpt environment to be initialized
!!     - qs_env   : Quickstep environment with the results of a
!!                  ground state calcualtion
!!     - ks_env   : Kohn-Sham environment already initialized
!!     - n_ev     : number of eigenvalues to calculate
!!     - glob_env : global environment
!!
!!   NOTES
!!     -
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, created
!!
!!****
  SUBROUTINE tddfpt_init(p_env, t_env, qs_env, ks_env, n_ev, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), INTENT(out)          :: p_env
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env 
    TYPE(qs_environment_type), INTENT(in)     :: qs_env
    TYPE(qs_ks_env_type), INTENT(in)          :: ks_env
    INTEGER                                   :: n_ev
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    TYPE(qs_K_P_P1_env_type), pointer         :: kpp1_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: psi0d
    TYPE(cp_full_matrix_type), POINTER        :: eigenvectors
    INTEGER                                   :: n_spins, spin
    CHARACTER(len=*), PARAMETER               :: routineN = "tddfpt_init", &
                                                 routineP = moduleN//"/"//routineN
    
    ! code

    CALL kpp1_create(kpp1_env, qs_env=qs_env)

    n_spins = SIZE(qs_env%k)
    ALLOCATE(psi0d(n_spins))

    DO spin=1, n_spins
       CALL get_mo_set(qs_env%C(spin)%mo_set, eigenvectors=eigenvectors)
       psi0d(spin)%matrix => eigenvectors
    END DO

    CALL p_env_init(p_env, qs_env, ks_env, kpp1_env, glob_env, &
                    psi0d, .TRUE.)

    CALL tddfpt_env_allocate(t_env, p_env, n_ev, glob_env)
    CALL tddfpt_env_init(t_env, p_env, n_ev, glob_env)

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

  SUBROUTINE update_p_env(p_env, t_env, glob_env)
    IMPLICIT NONE
    ! arguments
    TYPE(qs_p_env_type), INTENT(inout)        :: p_env
    TYPE(tddfpt_env_type), INTENT(in)         :: t_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    CHARACTER(len=*), PARAMETER  :: routineN = "update_p_env", &
                                    routineP = moduleN//"/"//routineN

    ! code
    IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A)') routineP//": Not yet implemented"

  END SUBROUTINE update_p_env

! *****************************************************************************

END MODULE qs_tddfpt_utils

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
  USE parser,                          ONLY: p_error, parser_end, parser_init, get_next, &
                                             test_next, read_line, stop_parser, search_label
  USE qs_blacs,                        ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             replicate_blacs_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_build_kernel,               ONLY: kpp1_create
  USE qs_p_types,                      ONLY: qs_kpp1_env_type,&
                                             qs_p_env_type
  USE qs_p_utils,                      ONLY: p_env_create,&
                                             p_op_l1
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_allocate,&
                                             tddfpt_env_init,&
                                             tddfpt_env_type, &
                                             tddfpt_control_type
  USE termination,                     ONLY: stop_memory
  USE string_utilities,                ONLY: uppercase
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: moduleN="qs_tddfpt_utils"
  INTEGER, PARAMETER          :: TOKEN_LENGTH=8
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  PRIVATE
  
  ! METHODS
  PUBLIC :: tddfpt_init, &
            tddfpt_write_banner, &
            update_p_env, &
            report_tddfpt_results, &
            read_tddfpt_control

CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_utils/tddfpt_init *
!!
!!   NAME
!!     tddfpt_init
!!
!!   SYNOPSIS
!!     Subroutine tddfpt_init(p_env, t_env, qs_env, n_ev, glob_env)
!!       Implicit None
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(tddfpt_env_type), Intent (OUT):: t_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Integer:: n_ev
!!       Type(global_environment_type), Intent (IN):: glob_env
!!     End Subroutine tddfpt_init
!!
!!   FUNCTION
!!     Initialize some necessary structures for a tddfpt calculation.
!!
!!   ARGUMENTS
!!     - p_env    : perturbation environment to be initialized
!!     - t_env    : tddfpt environment to be initialized
!!     - qs_env   : Quickstep environment with the results of a
!!                  ground state calcualtion
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
  SUBROUTINE tddfpt_init(p_env, t_env, qs_env, n_ev, glob_env)

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env 
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    INTEGER, INTENT(IN)                       :: n_ev
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! LOCALS
    TYPE(tddfpt_control_type), POINTER        :: t_control
    TYPE(qs_kpp1_env_type), POINTER           :: kpp1_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: psi0d
    TYPE(cp_full_matrix_type), POINTER        :: eigenvectors
    INTEGER                                   :: spin, stat
    CHARACTER(LEN=*), PARAMETER               :: routineN = "tddfpt_init", &
                                                 routineP = moduleN//"/"//routineN
    
    ! --------------------------------------------------------------------

    ALLOCATE(t_control, stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, "t_control", 1)
    t_control%n_ev    = n_ev ! set default value
    t_control%n_spins = qs_env%dft_control%nspins
    CALL read_tddfpt_control(glob_env, t_control)

    CALL kpp1_create(kpp1_env, qs_env=qs_env)

    ALLOCATE(psi0d(t_control%n_spins), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "psi0d", t_control%n_spins)

    DO spin=1, t_control%n_spins
       CALL get_mo_set(qs_env%C(spin)%mo_set, eigenvectors=psi0d(spin)%matrix)
    END DO

    CALL p_env_create(p_env, qs_env, kpp1_env, psi0d, .TRUE.)

    CALL tddfpt_env_allocate(t_env, t_control, p_env, qs_env, glob_env)
    CALL tddfpt_env_init(t_env, t_control, p_env, qs_env, glob_env)

    DEALLOCATE(psi0d, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "psi0d")

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
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(tddfpt_env_type), INTENT(in)         :: t_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    CHARACTER(len=*), PARAMETER  :: routineN = "update_p_env", &
                                    routineP = moduleN//"/"//routineN

    ! code
    IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A)') routineP//": Not yet implemented"

  END SUBROUTINE update_p_env

! *****************************************************************************

  SUBROUTINE report_tddfpt_results(t_env, glob_env)
    
    IMPLICIT NONE

    ! arguments
    TYPE(tddfpt_env_type), INTENT(in)         :: t_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    INTEGER :: ev, spin

    IF (.NOT.glob_env%ionode) RETURN ! only the io node may report

    WRITE (glob_env%scr, '(T2,A79)') &
         "", &
         "-------------------------------------------------------------------------------", &
         "-                               TDDFT Results                                 -", &
         "-------------------------------------------------------------------------------", &
         ""

    DO spin=1, t_env%control%n_spins
       WRITE (glob_env%scr, '(3X,A5,I3)') "SPIN ", spin
       WRITE (glob_env%scr, '(3X,A,T60,A)') "Nr", " numerical value"
       DO ev=1, t_env%control%n_ev
          WRITE (glob_env%scr, '(3X,I5,T60,F20.10)') ev, t_env%eval_x(spin)%values(ev,1)
       END DO
       WRITE (glob_env%scr, *)
    END DO
  
  END SUBROUTINE report_tddfpt_results

! *****************************************************************************

  SUBROUTINE read_tddfpt_control(glob_env, t_control)

    IMPLICIT NONE

    TYPE(global_environment_type), INTENT(in) :: glob_env
    TYPE(tddfpt_control_type), POINTER        :: t_control

    ! LOCALS
    CHARACTER(LEN=TOKEN_LENGTH)               :: label, token
    INTEGER                                   :: error, length

    CHARACTER(len=*), PARAMETER               :: routineN = "read_tddfpt_control", &
                                                 routineP = moduleN//"/"//routineN

    ! -------------------------------------------------------------------

    CALL parser_init(glob_env%input_file_name, glob_env)

    label = "&TDDFPT"
    CALL search_label(label, error)
    IF (error/=0) RETURN ! no &TDDFPT section found
    
    CALL read_line()

    length = TOKEN_LENGTH

    DO WHILE (test_next() /= 'X')
       CALL get_next(token, length)
       CALL uppercase(token)

       SELECT CASE (token)

       CASE DEFAULT

       CASE ('NEV', 'N_EV', 'EV') ! how many eigenvalues to compute
          IF (test_next() == 'N') THEN
             CALL get_next(t_control%n_ev)
             IF (t_control%n_ev < 1) &
                  CALL stop_parser(routineP, "number of eignvalues must be >= 1")
          ELSE
             CALL stop_parser(routineP, "NEV requires an integer argument")
          END IF

       END SELECT

       CALL read_line()
    END DO
    
    CALL parser_end()
    
  END SUBROUTINE read_tddfpt_control

! *****************************************************************************

END MODULE qs_tddfpt_utils

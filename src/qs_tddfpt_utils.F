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
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE qs_blacs,                        ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             replicate_blacs_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kpp1_env_methods,             ONLY: kpp1_create
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create, &
                                             p_op_l1, &
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: qs_p_env_type, p_env_release
  USE qs_tddfpt_types
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: moduleN="qs_tddfpt_utils"
  INTEGER, PARAMETER          :: TOKEN_LENGTH=8
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  PRIVATE
  
  ! METHODS
  PUBLIC :: tddfpt_cleanup, &
            tddfpt_init, &
            tddfpt_write_banner, &
            tddfpt_report_results, &
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
!!       Integer, Intent (IN):: n_ev
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

    CALL p_env_create(p_env, qs_env, kpp1_env=kpp1_env, psi0d=psi0d, orthogonal_orbitals=.TRUE.)
    ! update the m_epsilon matrix
    CALL p_env_psi0_changed(p_env, qs_env)

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

  SUBROUTINE tddfpt_report_results(t_env, glob_env)
    
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

    WRITE (glob_env%scr, '(3X,A,T35,A)') "Nr", " numerical value"
    DO ev=1, t_env%control%n_ev
       WRITE (glob_env%scr, '(3X,I5,T35,F20.10,A3,F20.10,A1)') &
            ev, t_env%ev(ev,1), " + ", t_env%ev(ev,2), "i"
    END DO
  
  END SUBROUTINE tddfpt_report_results

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

  SUBROUTINE tddfpt_cleanup(t_env, p_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(tddfpt_env_type), INTENT(inout)      :: t_env
    TYPE(global_environment_type), INTENT(IN) :: glob_env

    ! locals
    CHARACTER(LEN=*), PARAMETER :: routineN = "cleanup", &
                                   routineP = moduleN//'/'//routineN

    ! -------------------------------------------------------------------

    CALL tddfpt_env_deallocate(t_env, glob_env)
    CALL p_env_release(p_env)
    
  END SUBROUTINE tddfpt_cleanup

! *****************************************************************************

  FUNCTION check_eigenvalues(t_env, spin, glob_env, print_values_opt) RESULT (are_ok)

    IMPLICIT NONE

    ! FUNCTION DECLARATION
    LOGICAL                           :: are_ok

    ! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(IN)         :: t_env
    INTEGER, INTENT(IN)                       :: spin
    TYPE(global_environment_type), INTENT(IN) :: glob_env
    LOGICAL, INTENT(IN), OPTIONAL             :: print_values_opt

    ! LOCALS
    INTEGER                            :: ev
    LOGICAL                            :: print_values
    REAL(KIND=dbl)                     :: relative_difference
    REAL(DBL), DIMENSION(:,:), POINTER :: val_x, val_y
                                       
    CHARACTER(LEN=*), PARAMETER        :: routineN = 'check_eigenvalues', &
                                          routineP = moduleN//'/'//routineN

    ! -------------------------------------------------------------------

    val_x => t_env%eval_x(spin)%values
    val_y => t_env%eval_y(spin)%values

    IF (PRESENT(print_values_opt)) THEN
       print_values = print_values_opt
    ELSE
       print_values = .FALSE.
    END IF
    are_ok = .TRUE.
    
    DO ev=1, t_env%control%n_ev
       IF (ABS(val_x(ev,2)) < 1.0E-2_dbl) THEN

          !-------------!
          ! REAL NUMBER !
          !-------------!
          IF (ABS(val_x(ev,1)) > 1.0E-12_dbl) THEN
             relative_difference = ABS((val_x(ev,1)-val_y(ev,1))/val_x(ev,1))
          ELSEIF (ABS(val_y(ev,1)) > 1.0E-12_dbl) THEN
             relative_difference = ABS((val_x(ev,1)-val_y(ev,1))/val_y(ev,1))
          ELSE
             relative_difference = 0._dbl ! both are smaller than 1.0E-12 !!!
          END IF

          IF (relative_difference > 0.05_dbl) THEN
             are_ok = .FALSE.
             IF (print_values) THEN
                IF (glob_env%ionode) &
                     WRITE (glob_env%scr, '(1X,A,I4,A,ES11.4,2X,ES11.4)') &
                     "eigenvalues differ ", ev, ": ", val_x(ev,1), val_y(ev,1)
             END IF
          END IF

       ELSE

          !----------------!
          ! COMPLEX NUMBER !
          !----------------!
          are_ok = .FALSE.
          IF (print_values) THEN
             IF (glob_env%ionode) &
                  WRITE (glob_env%scr, '(1X,A,I4,A,ES11.4,A,ES11.4,A)') &
                  "complex eigenvalue ", ev, " is complex: ", &
                  val_x(ev,1), " + ", val_x(ev,2),"i"
          END IF

       END IF

    END DO

  END FUNCTION check_eigenvalues

! *****************************************************************************

END MODULE qs_tddfpt_utils

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_tddfpt_module *
!!
!!   NAME
!!     qs_tddfpt_module
!!
!!   FUNCTION
!!     Performs density functional perturbation theory (tddfpt) calculations.
!!     Uses the self consistent approach. The tddfpt calculation uses the ground
!!     state of the unperturbed system as the initial state.
!!
!!****

MODULE qs_tddfpt_module
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type,&
                                             optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_fm_utils,                     ONLY: copy_fm_col_to_fm
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_arpack,                       ONLY: arpack_diag_n
  USE qs_arpack_tddfpt,                ONLY: op_tddfpt_type,&
                                             setup_op_tddfpt
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_tddfpt_types,                 ONLY: real_2d_array_p_type,&
                                             tddfpt_control_type,&
                                             tddfpt_env_allocate,&
                                             tddfpt_env_deallocate,&
                                             tddfpt_env_type
  USE qs_tddfpt_utils
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  INTEGER, PARAMETER :: DEFAULT_N_EV = 1 ! default number of eigenvalues to calculate
  LOGICAL, PARAMETER :: DEBUG_THIS_MODULE = .TRUE.

  PRIVATE
  PUBLIC :: tddfpt_calculation

  CHARACTER(len=*), PARAMETER :: moduleN='qs_tddfpt_module'

CONTAINS

! *****************************************************************************

!!****f* cp2k/tddfpt_module/qs_tddfpt_calculation
!!
!!   NAME
!!     tddfpt_calculation
!!
!!   SYNOPSIS
!!     Subroutine tddfpt_calculation(qs_env, glob_env)
!!       Implicit None
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (IN):: glob_env
!!     End Subroutine tddfpt_calculation
!!
!!   FUNCTION
!!     Performs the perturbation calculation
!!
!!   ARGUMENTS
!!     - qs_env      - in    - qs environment with ground state
!!     - glob_env    - in    - global environment
!!
!!****
  SUBROUTINE tddfpt_calculation(qs_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in)        :: glob_env

    ! locals
    TYPE(qs_p_env_type), POINTER           :: p_env         ! perturbation environment
    TYPE(tddfpt_env_type)                  :: t_env         ! tddfpt environment
    TYPE(op_tddfpt_type)                   :: op            ! operator for tddfpt calculation
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: v ! temporary storage

    INTEGER                                :: stat
    INTEGER                                :: n_ev, n_spins ! cache values
    INTEGER                                :: ev, spin      ! loop indices
    LOGICAL                                :: evstat        ! status of the eigenvalues
    REAL(DBL), DIMENSION(:,:), POINTER     :: val_x, val_y  ! 

    CHARACTER(len=*), PARAMETER            :: routineN = 'tddfpt_calculation', &
                                              routineP = moduleN//'/'//routineN

    ! abort if this calculation is not wanted
    IF (.NOT.qs_env%dft_control%do_tddfpt_calculation) RETURN

    CALL tddfpt_write_banner(glob_env)

    CALL tddfpt_init(p_env, t_env, qs_env, DEFAULT_N_EV, glob_env)

    ! cache values
    n_ev    = t_env%control%n_ev
    n_spins = t_env%control%n_spins
    
    CALL create_v_matrix(v, t_env%control, p_env, qs_env, glob_env)

    CALL setup_op_tddfpt(op, p_env, qs_env, 'x', 1, glob_env)

    self_consistency: DO
       
       CALL update_p_env(p_env, t_env, glob_env)

       DO spin=1, n_spins

          op%spin = spin

          val_x => t_env%eval_x(spin)%values
          val_y => t_env%eval_y(spin)%values

          !--------------------!
          ! right eigenvectors !
          !--------------------!
          op%mode = 'x'; op%spin = spin
          CALL arpack_diag_n(n_ev, val_x, v(spin)%matrix, glob_env, op)
          DO ev=1, n_ev
             CALL copy_fm_col_to_fm(v(spin)%matrix, ev, t_env%x(ev)%spins(spin)%matrix)
          END DO

          !-------------------!
          ! left eigenvectors !
          !-------------------!
          op%mode = 'y'; op%spin = spin
          CALL arpack_diag_n(n_ev, val_y, v(spin)%matrix, glob_env, op)
          DO ev=1, n_ev
             CALL copy_fm_col_to_fm(v(spin)%matrix, ev, t_env%y(ev)%spins(spin)%matrix)
          END DO

          evstat = check_eigenvalues(t_env, spin, glob_env, print_values_opt=.TRUE.)
          
       END DO

       !---------------------------!
       ! test for self-consistency !
       !---------------------------!
       IF (.TRUE.) THEN
          IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A)') &
               routineP//": Convergence criteria not yet implemented"
          EXIT self_consistency
       END IF

    END DO self_consistency

    !-------------------!
    ! report the values !
    !-------------------!
    CALL report_tddfpt_results(t_env, glob_env)

    ! print the kohn-sham and overlap matrix

    !-------------------------!
    ! release some structures !
    !-------------------------!
    DO spin=1, n_spins
       CALL cp_fm_release(v(spin)%matrix)
    END DO
    DEALLOCATE(v, stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, 'v')
    CALL cleanup(t_env, p_env, glob_env)

    

  END SUBROUTINE tddfpt_calculation

! *****************************************************************************

  SUBROUTINE create_v_matrix(v, t_control, p_env, qs_env, glob_env)
    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: v
    TYPE(tddfpt_control_type), POINTER                 :: t_control
    TYPE(qs_p_env_type), INTENT(IN)                    :: p_env
    TYPE(qs_environment_type), INTENT(in), TARGET      :: qs_env
    TYPE(global_environment_type), INTENT(IN)          :: glob_env

    ! LOCALS
    TYPE(cp_fm_struct_type), POINTER :: v_struct
    INTEGER                          :: v_rows, v_cols
    INTEGER                          :: n_ev, n_spins
    INTEGER                          :: spin, stat
    CHARACTER(LEN=*), PARAMETER      :: routineN='create_v_matrix', &
                                        routineP=moduleN//'/'//routineN
    
    ! -------------------------------------------------------------------

    n_ev = t_control%n_ev
    n_spins = t_control%n_spins

    ALLOCATE(v(n_spins), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, 'v', n_spins)

    DO spin=1, n_spins
       v_rows = p_env%n_mo(spin) * p_env%n_ao(spin)
       v_cols = n_ev
       CALL cp_fm_struct_create(v_struct, glob_env%para_env, glob_env%blacs_env, &
            v_rows, v_cols, &
            optimal_blacs_row_block_size, optimal_blacs_col_block_size)
       CALL cp_fm_create(v(spin)%matrix, v_struct, 'tddfpt temporary v matrix')
    
    END DO

    ! Do not release v_struct !!!! It'll be released automatically
    ! by the call to cp_fm_release of the v matrix

  END SUBROUTINE create_v_matrix

! *****************************************************************************

  SUBROUTINE cleanup(t_env, p_env, glob_env)

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
    
  END SUBROUTINE cleanup

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

END MODULE qs_tddfpt_module

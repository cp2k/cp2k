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
!!     state of the unperturbed system as the initial state. So the ground state
!!     density has either to be provided to the subroutine or the be calculated
!!     before the tddfpt calculation can be started. This is why there are two 
!!     subroutines.
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     09-JUL-2002, TCH, development started
!!
!!   SEE ALSO
!!     
!!
!!****

MODULE qs_tddfpt_module

  USE cp_f_matrix_struct,   ONLY : cp_f_matrix_struct_type, &
                                   cp_fmstruct_create, cp_fmstruct_release, &
                                   optimal_blacs_col_block_size, optimal_blacs_row_block_size
  USE cp_full_matrix,       ONLY : cp_full_matrix_type, cp_full_matrix_p_type, &
                                   cp_fm_create, cp_fm_release
  USE global_types,         ONLY : global_environment_type
  USE kinds,                ONLY : dbl
  USE qs_arpack,            ONLY : arpack_diag_n
  USE qs_arpack_tddfpt,     ONLY : op_tddfpt_type, setup_op_tddfpt
  USE qs_build_ks_matrix,   ONLY : qs_ks_env_type
  USE qs_environment_types, ONLY : qs_environment_type
  USE qs_p_types,           ONLY : qs_p_env_type
  USE qs_tddfpt_utils
  USE qs_tddfpt_types,      ONLY : tddfpt_env_deallocate, tddfpt_env_allocate, &
                                   tddfpt_env_type,  &
                                   real_2d_array_p_type
  USE termination,          ONLY : stop_memory

  IMPLICIT NONE

  integer, parameter :: n_ev = 1 ! the number of eigenvalues to calculate
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
!!   FUNCTION
!!     Performs the perturbation calculation
!!
!!   ARGUMENTS
!!     - qs_env      - in    - qs environment with ground state
!!     - ks_env      - in    - KS matrix environment
!!     - glob_env    - in    - global environment
!!
!!   MODIFICATION HISTORY
!!     , TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE tddfpt_calculation(qs_env, ks_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_ks_env_type), INTENT(in)          :: ks_env
    TYPE(qs_environment_type), INTENT(in)     :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    TYPE(qs_p_env_type)                    :: p_env
    TYPE(tddfpt_env_type)                  :: t_env
    INTEGER                                :: i, n_spins, spin, stat
    TYPE(op_tddfpt_type)                   :: op
    REAL(DBL), DIMENSION(:,:), POINTER     :: val_x, val_y
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: v
    CHARACTER(len=*), PARAMETER            :: routineN = 'tddfpt_calculation', &
                                              routineP = moduleN//'/'//routineN

    ! abort if this calculation is not wanted
    IF (.NOT.qs_env%dft_control%do_tddfpt_calculation) RETURN

    CALL tddfpt_write_banner(glob_env)

    CALL tddfpt_init(p_env, t_env, qs_env, ks_env, n_ev, glob_env)

    CALL create_v_matrix(v, p_env, n_ev, glob_env)

    n_spins = p_env%qs_env%dft_control%nspins
    
    CALL setup_op_tddfpt(op, p_env, 'x', n_ev, glob_env)

    self_consistency: DO
       
       CALL update_p_env(p_env, t_env, glob_env)

       ! solve eigenequation
       DO spin=1, n_spins

          op%spin = spin

          val_x => t_env%eval_x(spin)%values
          val_y => t_env%eval_y(spin)%values

          ! right eigenvectors
          op%mode = 'x'; op%spin = spin
          CALL arpack_diag_n(n_ev, val_x, v(spin)%matrix, glob_env, op)
          ! copy v back to x

          ! left eigenvectors
          op%mode = 'y'; op%spin = spin
          CALL arpack_diag_n(n_ev, val_y, v(spin)%matrix, glob_env, op)
          ! copy v back to y
          
          IF (DEBUG_THIS_MODULE) THEN
          ! compare the eigenvalues
          DO i=1, SIZE(val_x, 1)
             IF (ABS(val_x(i,2)) < 1.0E-2_dbl) THEN
                ! real number
                IF (ABS(val_x(i,1)) > 1.0E-12_dbl) THEN
                   IF (ABS((val_x(i,1)-val_y(i,1))/val_x(i,1)) > 0.05_dbl) THEN
                      IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A,I4,ES11.4,2X,ES11.4)') &
                           routineP//": Eigenvalue ", i, ": ", val_x(i,1), val_y(i,1)
                   END IF
                END IF
             ELSE
                ! complex number!
                IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A,I4,A,ES11.4,A,ES11.4,A)') &
                     routineP//": Eigenvalue ", i, " is complex: ", &
                     val_x(i,1), " + ", val_x(i,2),"i"
             END IF
          END DO
          END IF ! DEBUG_THIS_MODULE
       END DO

       ! test for convergen
       IF (.TRUE.) THEN
          IF (glob_env%ionode) WRITE (glob_env%scr, '(1X,A)') &
               routineP//": Convergence criteria not yet implemented"
          EXIT self_consistency
       END IF

    END DO self_consistency

    DO spin=1, n_spins
       CALL cp_fm_release(v(spin)%matrix)
    END DO
    DEALLOCATE(v, stat=stat)
    if (stat/=0) CALL stop_memory(routineP, 'v')

  END SUBROUTINE tddfpt_calculation

! *****************************************************************************

  SUBROUTINE create_v_matrix(v, p_env, n_ev, glob_env)
    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: v
    TYPE(qs_p_env_type), INTENT(IN)                    :: p_env
    INTEGER, INTENT(IN)                                :: n_ev
    TYPE(global_environment_type), INTENT(IN)          :: glob_env

    ! locals
    TYPE(cp_f_matrix_struct_type), POINTER :: v_struct
    INTEGER                                :: v_rows, v_cols
    INTEGER                                :: spin, n_spins, stat
    CHARACTER(LEN=*), PARAMETER            :: routineN='create_v_matrix', &
                                              routineP=moduleN//'/'//routineN
    
    ! code
    IF (DEBUG_THIS_MODULE) WRITE(*,*) routineP//": entered"

    n_spins = p_env%qs_env%dft_control%nspins

    ALLOCATE(v(n_spins), stat=stat)
    IF (stat/=0) CALL stop_memory(routineP, 'v and/or v_struc', n_spins)

    DO spin=1, n_spins
       v_rows = p_env%n_mo(spin) * p_env%n_ao(spin)
       v_cols = n_ev
       call cp_fmstruct_create(v_struct, glob_env%para_env, glob_env%blacs_env, &
            v_rows, v_cols, &
            optimal_blacs_row_block_size, optimal_blacs_col_block_size)
       call cp_fm_create(v(spin)%matrix, v_struct, 'tddfpt temporary v matrix')
    
    END DO

    ! Do not release v_struct !!!! It'll be released automatically
    ! by the call to cp_fm_release of the v matrix

    IF (DEBUG_THIS_MODULE) WRITE(*,*) routineP//": finished"

  END SUBROUTINE create_v_matrix

! *****************************************************************************

  SUBROUTINE cleanup(t_env, p_env, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(tddfpt_env_type), POINTER            :: t_env
    TYPE(global_environment_type), INTENT(IN) :: glob_env

    ! locals
    CHARACTER(LEN=*), PARAMETER :: routineN = "cleanup", &
                                   routineP = moduleN//'/'//routineN

    ! code
    IF (DEBUG_THIS_MODULE) WRITE(*,*) routineP//": entered"

    CALL tddfpt_env_deallocate(t_env, glob_env)

    ! clean p_env also?
    
    IF (DEBUG_THIS_MODULE) WRITE(*,*) routineP//": finished"

  END SUBROUTINE cleanup

! *****************************************************************************

END MODULE qs_tddfpt_module

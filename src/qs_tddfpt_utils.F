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
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose, &
                                             cp_fm_cholesky_invert
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
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type, &
                                             cp_full_matrix_type, &
                                             cp_fm_get_element, &
                                             cp_fm_init_random
  USE qs_blacs,                        ONLY: replicate_blacs_matrix, &
                                             copy_full_matrix_to_real_matrix, &
                                             copy_real_matrix_to_full_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kpp1_env_methods,             ONLY: kpp1_create
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create, &
                                             p_env_psi0_changed, p_preortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type, p_env_release
  USE qs_tddfpt_types
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory, stop_program
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: moduleN="qs_tddfpt_utils"
  INTEGER, PARAMETER          :: TOKEN_LENGTH=8
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  PRIVATE
  
  ! METHODS
  PUBLIC :: tddfpt_cleanup, &
            tddfpt_init, &
            tddfpt_write_banner, &
            print_fm_matrix, &
            tddfpt_env_init

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
!!     - glob_env : global environment
!!
!!   NOTES
!!     -
!!
!!****
  SUBROUTINE tddfpt_init(p_env, t_env, qs_env, glob_env)

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(qs_p_env_type), POINTER              :: p_env
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env 
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! LOCALS
    TYPE(qs_kpp1_env_type), POINTER           :: kpp1_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: psi0d
    TYPE(cp_full_matrix_type), POINTER        :: eigenvectors
    INTEGER                                   :: spin, n_spins, stat
    CHARACTER(LEN=*), PARAMETER               :: routineN = "tddfpt_init", &
                                                 routineP = moduleN//"/"//routineN
    
    ! --------------------------------------------------------------------

    n_spins   = qs_env%dft_control%nspins

    !------------------!
    ! create the p_env !
    !------------------!
    CALL kpp1_create(kpp1_env, qs_env=qs_env)
    ALLOCATE(psi0d(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "psi0d", n_spins)
    DO spin=1, n_spins
       CALL get_mo_set(qs_env%C(spin)%mo_set, eigenvectors=psi0d(spin)%matrix)
    END DO
    CALL p_env_create(p_env, qs_env, kpp1_env=kpp1_env, psi0d=psi0d, &
                      orthogonal_orbitals=.TRUE.)
    CALL p_env_psi0_changed(p_env, qs_env) ! update the m_epsilon matrix


    !------------------!
    ! create the t_env !
    !------------------!
    CALL tddfpt_env_allocate(t_env, p_env, qs_env, glob_env)
    CALL tddfpt_env_init(t_env, p_env, qs_env, glob_env)


    !----------!
    ! clean up !
    !----------!
    DEALLOCATE(psi0d, stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, "psi0d")

  END SUBROUTINE tddfpt_init

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_init *
!!
!!   FUNCTION
!!     Initialize t_env with meaningfull values.
!!
!!***
  SUBROUTINE tddfpt_env_init(t_env, p_env, qs_env, glob_env)

    IMPLICIT NONE
  
  	! ARGUMENTS
    TYPE(tddfpt_env_type), INTENT(out)        :: t_env
    TYPE(qs_p_env_type),POINTER               :: p_env
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! LOCALS
    INTEGER                            :: stat
    INTEGER                            :: ev, spin, n_spins
    TYPE(cp_full_matrix_type), POINTER :: fm
    

    CHARACTER(len=*), PARAMETER        :: routineN = 'tddfpt_env_init', &
                                          routineP = moduleN//'/'//routineN

    ! -------------------------------------------------------------------

    n_spins = qs_env%dft_control%nspins

    !------------------!
    ! initialize the x !
    !------------------!
    ev = 1 ! we only need the first vector as starting point
    DO spin=1, n_spins
       CALL cp_fm_init_random(t_env%X(ev)%spins(spin)%matrix, p_env%n_mo(spin))
    END DO

    CALL p_preortho(p_env, qs_env, t_env%X(ev)%spins)

    !-------------------------------------------------------!
    ! Compute the inverse of S using cholesky factorization !
    !-------------------------------------------------------!
    DO spin=1, n_spins
       CALL copy_real_matrix_to_full_matrix(qs_env%s(spin)%matrix, t_env%invS(spin)%matrix)
       CALL cp_fm_cholesky_decompose(t_env%invS(spin)%matrix)
       CALL cp_fm_cholesky_invert(t_env%invS(spin)%matrix)
    END DO

  END SUBROUTINE tddfpt_env_init

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

  SUBROUTINE print_fm_matrix(matrix, imax, jmax, name, sym)

    IMPLICIT NONE

    TYPE(cp_full_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)                :: imax, jmax
    CHARACTER(*), INTENT(IN)           :: name
    CHARACTER, OPTIONAL, INTENT(IN)    :: sym

    INTEGER   :: i, j, ioerr
    REAL(DBL) :: element
    CHARACTER :: symmetric

    symmetric = 'F'
    IF (PRESENT(sym)) symmetric=sym

    OPEN(UNIT=8, FILE=name, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)

    DO i=1, imax
       DO j=1, jmax
          IF ((symmetric == 'T') .AND. (i>j)) THEN
             CALL cp_fm_get_element(matrix, j, i, element)
          ELSE
             CALL cp_fm_get_element(matrix, i, j, element)
          END IF
          WRITE (8, '(1X,F20.10)') element
       END DO
    END DO

    CLOSE(8)

  END SUBROUTINE print_fm_matrix

! *****************************************************************************

END MODULE qs_tddfpt_utils

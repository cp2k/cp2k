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

  USE cp_fm_basic_linalg,              ONLY: cp_sm_fm_multiply, &
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose, &
                                             cp_fm_cholesky_invert
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type, &
                                             cp_full_matrix_type, &
                                             cp_fm_get_element, &
                                             cp_fm_get_submatrix, &
                                             cp_fm_init_random, &
                                             cp_fm_set_all, &
                                             cp_fm_set_submatrix
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_blacs,                        ONLY: replicate_blacs_matrix, &
                                             copy_full_matrix_to_real_matrix, &
                                             copy_real_matrix_to_full_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kpp1_env_methods,             ONLY: kpp1_create
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create, &
                                             p_preortho, &
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: qs_p_env_type, &
                                             p_env_release
  USE qs_tddfpt_types
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_memory, stop_program
  USE timings,                         ONLY: timeset, timestop

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
            tddfpt_env_init, &
            calc_initial_guess, &
            calc_norm, &
            reorthogonalize

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
       
       WRITE (unit=globenv%scr, fmt="(2X,A)")&
       "",&
       "-----------------------------------------------------------------------------",&
       "-                    Time-dependent calculation started                     -",&
       "-----------------------------------------------------------------------------",&
       ""    

    END IF

  END SUBROUTINE tddfpt_write_banner

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
    LOGICAL   :: local

    symmetric = 'F'
    IF (PRESENT(sym)) symmetric=sym

    OPEN(UNIT=8, FILE=name, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)

    DO i=1, imax
       DO j=1, jmax
          IF ((symmetric == 'T') .AND. (i>j)) THEN
             CALL cp_fm_get_element(matrix, j, i, element, local)
             ! this is to get the old meaning of the call, unclear if this is really meant to be the case ! Joost
             IF (.NOT. local) element = 0.0_dbl
          ELSE
             CALL cp_fm_get_element(matrix, i, j, element, local)
             ! this is to get the old meaning of the call, unclear if this is really meant to be the case ! Joost
             IF (.NOT. local) element = 0.0_dbl
          END IF
          WRITE (8, '(1X,F20.10)') element
       END DO
    END DO

    CLOSE(8)

  END SUBROUTINE print_fm_matrix

! *****************************************************************************

  SUBROUTINE calc_initial_guess(matrices, qs_env, p_env, n_kv)

    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER  :: matrices
    TYPE(qs_environment_type), INTENT(IN),TARGET        :: qs_env
    TYPE(qs_p_env_type), POINTER                        :: p_env
    INTEGER, INTENT(IN)                                 :: n_kv

    ! locals
    TYPE(tddfpt_control_type), POINTER          :: tddfpt_control
    REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: lumos, one_lumo
    REAL(KIND=DBL)                              :: scale_factor
    INTEGER                                     :: spin, n_spins
    INTEGER                                     :: n_rows, n_cols ! of the matrices
    INTEGER                                     :: n_ev, n_lumos
    INTEGER                                     :: stat           ! error flags and such
    INTEGER                                     :: i, j           ! loop indices

    CHARACTER(len=*), PARAMETER :: routineN="calc_initial_guess", &
                                   routineP=moduleN//"/"//routineN


    tddfpt_control => qs_env%dft_control%tddfpt_control
    n_spins = qs_env%dft_control%nspins
    n_ev    = tddfpt_control%n_ev

    DO spin=1, n_spins

       n_cols = matrices(spin)%matrix%matrix_struct%ncol_global
       n_rows = matrices(spin)%matrix%matrix_struct%nrow_global

       CALL cp_fm_set_all(matrices(spin)%matrix, 0.0_dbl)

       IF (ASSOCIATED(tddfpt_control%lumos)) THEN       

          n_lumos = tddfpt_control%lumos(spin)%matrix%matrix_struct%ncol_global

          ALLOCATE(lumos(n_rows, n_lumos), one_lumo(n_rows, 1), stat=stat)

          CALL cp_fm_get_submatrix(tddfpt_control%lumos(spin)%matrix, lumos, &
                                   start_col=1, n_cols=n_lumos)

          DO i=1, n_ev
             CALL dcopy(n_rows, 0.0_dbl, 0, one_lumo(:,1), 1)
             DO j=1, n_lumos
                scale_factor = 0.8_dbl
                IF (j==i) scale_factor = REAL(n_ev,dbl)
                scale_factor = scale_factor/(REAL(i,dbl)**2)
                CALL daxpy(n_rows, scale_factor, lumos(:,j), 1, one_lumo(:,1), 1)
             END DO
             CALL cp_fm_set_submatrix(matrices(spin)%matrix, &
                                      one_lumo(:,1:1), &
                                      start_col=(n_cols-i+1))
          END DO

          DEALLOCATE(lumos, one_lumo)

       ELSE

          CALL cp_fm_init_random(matrices(spin)%matrix,  n_kv, start_col=n_cols-n_kv+1)

       END IF

    END DO
    CALL p_preortho(p_env, qs_env, matrices)
    ! normalization not necessary since it is done in the krylov routine

  END SUBROUTINE calc_initial_guess

! *****************************************************************************

  !---------------------------------------!
  ! x must not be changed in this routine !
  ! tmp_vec may be changed                !
  !---------------------------------------!
  SUBROUTINE calc_norm(norm, x, tmp_vec, qs_env, glob_env)

    IMPLICIT NONE

    REAL(KIND=dbl), INTENT(OUT)                        :: norm
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: x, tmp_vec
    TYPE(qs_environment_type)                          :: qs_env
    TYPE(global_environment_type)                      :: glob_env

    INTEGER        :: spin, n_spins
    REAL(KIND=dbl) :: spin_norm

    n_spins = SIZE(x)

    !-------------------------!
    ! calc beta=transp(X)*S*X !
    !-------------------------!
    norm = 0.0_dbl
    DO spin=1, n_spins
       CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, x(spin)%matrix, &
            tmp_vec(spin)%matrix, &
            x(spin)%matrix%matrix_struct%ncol_global, &
            glob_env%para_env, 1.0_dbl, 0.0_dbl)
       CALL cp_fm_trace(tmp_vec(spin)%matrix, x(spin)%matrix, spin_norm)
       norm = norm + spin_norm
    END DO

  END SUBROUTINE calc_norm

! *****************************************************************************

  ! X        : the vector to reorthogonalize
  ! kv       : set of vectors to which to reorthogonalize
  ! S        : metric of the vectorspaceA
  ! tmp_mat    : temporary matrix
  ! j        : up to which vector of kv to reorthogonalize
  ! glob_env : for communication

  SUBROUTINE reorthogonalize(X, kv, S, tmp_mat, j, glob_env)

    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER   :: X, tmp_mat
    TYPE(cp_full_matrix_p_type), DIMENSION(:,:), POINTER :: kv
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER      :: S
    INTEGER, INTENT(IN)                                  :: j
    TYPE(global_environment_type), INTENT(in)            :: glob_env

    ! locals
    INTEGER :: i, index, spin, time_handle, n_spins
    REAL(DBL) :: tmp, dot_product

    ! parameters
    CHARACTER(LEN=*), PARAMETER :: routineN = "reorthogonalize", &
                                   routineP = moduleN//"/"//routineN

    CALL timeset(routineN,"I","",time_handle)

    n_spins = SIZE(X)

    IF (j > 1) THEN
       DO i=1, j
          dot_product = 0.0_dbl
          DO spin=1, n_spins
             CALL cp_sm_fm_multiply(S(spin)%matrix, X(spin)%matrix, &
                                    tmp_mat(spin)%matrix, &
                                    X(spin)%matrix%matrix_struct%ncol_global, &
                                    glob_env%para_env, 1.0_dbl, 0.0_dbl)
             CALL cp_fm_trace(tmp_mat(spin)%matrix, kv(i,spin)%matrix, tmp)
             dot_product = dot_product + tmp
          END DO
          DO spin=1, n_spins
             CALL cp_fm_scale_and_add(1.0_dbl, X(spin)%matrix, &
                                      -1.0_dbl*dot_product, kv(i,spin)%matrix)
          END DO
       END DO
    END IF

    CALL timestop(0.0_dbl, time_handle)

  END SUBROUTINE reorthogonalize

! *****************************************************************************

END MODULE qs_tddfpt_utils

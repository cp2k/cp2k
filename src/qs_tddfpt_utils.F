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

  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_kpp1_env_methods,             ONLY: kpp1_create
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_tddfpt_types
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: moduleN="qs_tddfpt_utils"
  INTEGER, PARAMETER          :: TOKEN_LENGTH=8
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  TYPE simple_solution_sorter
     INTEGER                               :: orbit
     INTEGER                               :: lumo
     REAL(KIND=DBL)                        :: value
     TYPE(simple_solution_sorter), POINTER :: next
  END TYPE simple_solution_sorter

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


    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(tddfpt_env_type), INTENT(out)       :: t_env
    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: glob_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "tddfpt_init", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: n_spins, spin, stat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi0d
    TYPE(cp_fm_type), POINTER                :: eigenvectors
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env

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

  
    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: glob_env

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_init', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: n_spins, spin

! -------------------------------------------------------------------

    n_spins = qs_env%dft_control%nspins

    !-------------------------------------------------------!
    ! Compute the inverse of S using cholesky factorization !
    !-------------------------------------------------------!
    DO spin=1, n_spins
       CALL copy_sm_to_fm(qs_env%s(spin)%matrix, t_env%invS(spin)%matrix)
       CALL cp_fm_cholesky_decompose(t_env%invS(spin)%matrix)
       CALL cp_fm_cholesky_invert(t_env%invS(spin)%matrix)
    END DO

  END SUBROUTINE tddfpt_env_init

! *****************************************************************************

  SUBROUTINE tddfpt_write_banner(globenv)

    
    TYPE(global_environment_type), &
      INTENT(in)                             :: globenv

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


    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: glob_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "cleanup", &
      routineP = moduleN//'/'//routineN

! -------------------------------------------------------------------

    CALL tddfpt_env_deallocate(t_env, glob_env)
    CALL p_env_release(p_env)
    
  END SUBROUTINE tddfpt_cleanup

! *****************************************************************************

  SUBROUTINE print_fm_matrix(matrix, imax, jmax, name, sym)


    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN)                      :: imax, jmax
    CHARACTER(*), INTENT(IN)                 :: name
    CHARACTER, INTENT(IN), OPTIONAL          :: sym

    CHARACTER                                :: symmetric
    INTEGER                                  :: i, ioerr, j
    LOGICAL                                  :: local
    REAL(DBL)                                :: element

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

  SUBROUTINE calc_initial_guess(matrices, qs_env, p_env, energies)


    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrices
    TYPE(qs_environment_type), INTENT(IN), &
      TARGET                                 :: qs_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    REAL(kind=DBL), DIMENSION(:, 0:), &
      INTENT(OUT)                            :: energies

    CHARACTER(len=*), PARAMETER :: routineN = "calc_initial_guess", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: i, j, n_cols, n_ev, n_lumos, &
                                                n_orbits, n_rows, n_spins, &
                                                spin, stat
    REAL(KIND=DBL)                           :: scale_factor, evd
    REAL(KIND=DBL), ALLOCATABLE, &
      DIMENSION(:, :)                        :: lumos, tmp_guess
    REAL(KIND=DBL), DIMENSION(:), POINTER    :: orbital_eigenvalues
    TYPE(simple_solution_sorter), POINTER    :: sorter_iterator, &
                                                sorter_pointer, sorter_start
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

    tddfpt_control => qs_env%dft_control%tddfpt_control
    n_spins = qs_env%dft_control%nspins
    n_ev    = tddfpt_control%n_ev

    DO spin=1, n_spins

       n_cols = matrices(spin)%matrix%matrix_struct%ncol_global
       n_rows = matrices(spin)%matrix%matrix_struct%nrow_global

       CALL cp_fm_set_all(matrices(spin)%matrix, 0.0_dbl)

       IF (ASSOCIATED(tddfpt_control%lumos)) THEN       

          CALL get_mo_set(qs_env%C(spin)%mo_set, eigenvalues=orbital_eigenvalues)

          n_lumos = tddfpt_control%lumos(spin)%matrix%matrix_struct%ncol_global

          n_orbits = SIZE(orbital_eigenvalues)

          ! first element
          evd = tddfpt_control%lumos_eigenvalues(1,spin)-orbital_eigenvalues(n_orbits)
          ALLOCATE(sorter_start)
          sorter_start%orbit=n_orbits
          sorter_start%lumo=1
          sorter_start%value=evd
          NULLIFY(sorter_start%next)
         
          DO i=n_orbits, 1, -1
             DO j=1, MIN(n_ev, n_lumos)

                IF (i==n_orbits .AND. j==1) CYCLE ! already in list

                evd = tddfpt_control%lumos_eigenvalues(j,spin)-orbital_eigenvalues(i)

                sorter_iterator => sorter_start
                NULLIFY(sorter_pointer)
                DO WHILE (ASSOCIATED(sorter_iterator%next))
                   IF (sorter_iterator%next%value > evd) THEN
                      sorter_pointer => sorter_iterator%next
                      EXIT
                   END IF
                   sorter_iterator => sorter_iterator%next
                END DO

                ALLOCATE(sorter_iterator%next)
                sorter_iterator%next%orbit = i
                sorter_iterator%next%lumo  = j
                sorter_iterator%next%value = evd
                sorter_iterator%next%next  => sorter_pointer

             END DO
          END DO

          ALLOCATE(lumos(n_rows, n_lumos), tmp_guess(n_rows, n_orbits), stat=stat)

          CALL cp_fm_get_submatrix(tddfpt_control%lumos(spin)%matrix, lumos, &
                                   start_col=1, n_cols=n_lumos)

          
          CALL dcopy(n_rows*n_orbits, 0.0_dbl, 0, tmp_guess, 1)          
          sorter_iterator=>sorter_start

          scale_factor = 1/SQRT(REAL(n_ev,DBL))
          DO i=1, MIN(n_ev, n_lumos)
             CALL daxpy(n_rows, scale_factor/REAL(i,dbl), &
                        lumos(:,sorter_iterator%lumo), 1, &
                        tmp_guess(:,sorter_iterator%orbit),1)
             energies(i,0) = sorter_iterator%value
             sorter_iterator=>sorter_iterator%next
          END DO
          CALL cp_fm_set_submatrix(matrices(spin)%matrix, &
                                   tmp_guess(:,1:n_orbits))

          DEALLOCATE(lumos, tmp_guess)

       ELSE

          CALL cp_fm_init_random(matrices(spin)%matrix, n_ev, start_col=n_cols-n_Ev+1)
          DO i=1, n_ev
             energies(i,0)=9.999E99_dbl
          END DO

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


    REAL(KIND=dbl), INTENT(OUT)              :: norm
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: x, tmp_vec
    TYPE(qs_environment_type)                :: qs_env
    TYPE(global_environment_type)            :: glob_env

    INTEGER                                  :: n_spins, spin
    REAL(KIND=dbl)                           :: spin_norm

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


    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: X
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: kv
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: S
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: tmp_mat
    INTEGER, INTENT(IN)                      :: j
    TYPE(global_environment_type), &
      INTENT(in)                             :: glob_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "reorthogonalize", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: i, n_spins, spin, time_handle
    REAL(DBL)                                :: dot_product, tmp

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

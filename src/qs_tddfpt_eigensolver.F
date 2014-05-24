!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE qs_tddfpt_eigensolver
  USE cp_control_types,                ONLY: tddfpt_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_trace
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_p_type,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm
  USE f77_blas
  USE input_constants,                 ONLY: tddfpt_davidson,&
                                             tddfpt_lanczos
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_p_env_methods,                ONLY: p_op_l1,&
                                             p_op_l2,&
                                             p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_type
  USE qs_tddfpt_utils,                 ONLY: co_initial_guess,&
                                             normalize,&
                                             reorthogonalize
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_eigensolver'

  PRIVATE

  PUBLIC :: apply_op, eigensolver

CONTAINS

! *****************************************************************************
  SUBROUTINE eigensolver(p_env, qs_env, t_env, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(tddfpt_env_type), INTENT(INOUT)     :: t_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eigensolver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_ev, nspins, &
                                                output_unit, restarts, stat
    LOGICAL                                  :: do_kernel_save, failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ievals
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    NULLIFY(logger)
    failure   = .FALSE.

    logger    => cp_error_get_logger(error)
    output_unit  = cp_logger_get_default_io_unit(logger)

    n_ev      = qs_env%dft_control%tddfpt_control%n_ev
    nspins    = qs_env%dft_control%nspins

    ALLOCATE(ievals(n_ev), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !---------------!
    ! initial guess !
    !---------------!
    do_kernel_save = qs_env%dft_control%tddfpt_control%do_kernel
    qs_env%dft_control%tddfpt_control%do_kernel = .FALSE.
    IF (output_unit>0) THEN
       WRITE (output_unit,*) " Generating initial guess"
       WRITE (output_unit,*)
    END IF
    IF (ASSOCIATED(qs_env%dft_control%tddfpt_control%lumos)) THEN
       CALL co_initial_guess(t_env%evecs, ievals, n_ev, qs_env, error=error)
    ELSE
       IF (output_unit>0) WRITE (output_unit,*) "LUMOS are needed in TDDFPT!"
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    DO restarts=1, qs_env%dft_control%tddfpt_control%n_restarts
       IF (iterative_solver(ievals, t_env, p_env, qs_env, ievals, error=error)) EXIT
       IF (output_unit>0) THEN
          WRITE (output_unit,*) " Restarting"
          WRITE (output_unit,*)
       END IF
    END DO
    qs_env%dft_control%tddfpt_control%do_kernel = do_kernel_save

    !-----------------!
    ! call the solver !
    !-----------------!
    IF (output_unit>0) THEN
       WRITE (output_unit,*)
       WRITE (output_unit,*) " Doing TDDFPT calculation"
       WRITE (output_unit,*)
    END IF
    DO restarts=1, qs_env%dft_control%tddfpt_control%n_restarts
          IF (iterative_solver(ievals, t_env, p_env, qs_env, t_env%evals, error=error)) EXIT
       IF (output_unit>0) THEN
          WRITE (output_unit,*) " Restarting"
          WRITE (output_unit,*)
       END IF
    END DO

    !---------!
    ! cleanup !
    !---------!
    DEALLOCATE(ievals, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE eigensolver

  ! in_evals  : approximations to the eigenvalues for the preconditioner
  ! t_env     : TD-DFT environment values
  ! p_env     : perturbation environment values
  ! qs_env    : general Quickstep environment values
  ! out_evals : the resulting eigenvalues
  ! error     : used for error handling
  !
  ! res       : the function will return wheter the eigenvalues are converged or not

! *****************************************************************************
  FUNCTION iterative_solver(in_evals, &
                            t_env, p_env, qs_env, &
                            out_evals, error) RESULT(res)

    REAL(KIND=dp), DIMENSION(:)              :: in_evals
    TYPE(tddfpt_env_type), INTENT(INOUT)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(kind=dp), DIMENSION(:), OPTIONAL    :: out_evals
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'iterative_solver', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mode
    INTEGER :: col, handle, i, iev, iter, j, k, max_krylovspace_dim, max_kv, &
      n_ev, n_kv, nspins, output_unit, row, spin, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: must_improve
    LOGICAL                                  :: failure
    REAL(dp)                                 :: Atilde_ij, convergence, tmp, &
                                                tmp2
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_difference, evals_tmp
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evals
    TYPE(cp_fm_p_type)                       :: Atilde, Us
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: R, X
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: Ab, b, Sb
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_p_type), &
      DIMENSION(:), POINTER                  :: kv_fm_struct
    TYPE(cp_fm_struct_type), POINTER         :: tilde_fm_struct
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

    failure   = .FALSE.
    res       = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(ao_mo_fm_pools, X, R, b, Ab, Sb, tddfpt_control, logger, &
            tilde_fm_struct, Atilde%matrix, Us%matrix)

    CALL get_qs_env(qs_env, tddfpt_control=tddfpt_control,error=error)
    logger    => cp_error_get_logger(error)
    output_unit  = cp_logger_get_default_io_unit(logger)
    n_ev      = tddfpt_control%n_ev
    nspins    = qs_env%dft_control%nspins

    IF (qs_env%dft_control%tddfpt_control%diag_method == tddfpt_lanczos) THEN
       mode = 'L'
    ELSE IF (qs_env%dft_control%tddfpt_control%diag_method == tddfpt_davidson) THEN
       mode = 'D'
    END IF

    !-----------------------------------------!
    ! determine the size of the problem       !
    ! and how many krylov space vetors to use !
    !-----------------------------------------!
    max_krylovspace_dim = SUM(p_env%n_ao(1:nspins)*p_env%n_mo(1:nspins))
    max_kv = tddfpt_control%max_kv
    IF (max_krylovspace_dim <= max_kv) THEN
       max_kv = max_krylovspace_dim
       IF (output_unit>0) THEN
          WRITE (output_unit,*) "  Setting the maximum number of krylov vectors to ", max_kv, "!!"
       END IF
    END IF
!    CPPostcondition(max_krylovspace_dim>=max_kv,cp_failure_level,routineP,error,failure)

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, X, name=routineP//":X",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, R, name=routineP//":R",error=error)

    ALLOCATE(evals_difference(n_ev), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(must_improve(n_ev), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(evals(max_kv, 0:max_kv), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(evals_tmp(max_kv), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(b(max_kv,nspins), Ab(max_kv,nspins), &
             Sb(max_kv,nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(kv_fm_struct(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO spin=1, nspins
       NULLIFY (kv_fm_struct(spin)%struct)
       CALL cp_fm_struct_create(kv_fm_struct(spin)%struct, qs_env%para_env, &
                                qs_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_mo(spin),error=error)
    END DO
    DO spin=1, nspins
       DO i=1, max_kv
          NULLIFY(b(i,spin)%matrix, Ab(i,spin)%matrix, Sb(i,spin)%matrix)
       END DO
    END DO

    IF (output_unit>0) THEN
       WRITE (output_unit, '(2X,A,T69,A)') &
              "nvec", "Convergence"
       WRITE (output_unit, '(2X,A)') &
              "-----------------------------------------------------------------------------"
    END IF

    iter = 1
    k    = 0
    n_kv = n_ev
    iteration: DO

       CALL allocate_krylov_vectors(b,  "b-",  k+1, n_kv, nspins, kv_fm_struct,error=error)
       CALL allocate_krylov_vectors(Ab, "Ab-", k+1, n_kv, nspins, kv_fm_struct,error=error)
       CALL allocate_krylov_vectors(Sb, "Sb-", k+1, n_kv, nspins, kv_fm_struct,error=error)

       DO i=1, n_kv
          k=k+1

          IF (k <= SIZE(t_env%evecs,1)) THEN ! the first iteration

             ! take the initial guess
             DO spin=1, nspins
                CALL cp_fm_to_fm(t_env%evecs(k,spin)%matrix, b(k,spin)%matrix,error=error)
             END DO

          ELSE

             ! create a new vector
             IF (mode == 'L') THEN

                DO spin=1, nspins
                   IF (tddfpt_control%invert_S) THEN
                      CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                           1.0_dp, t_env%invS(spin)%matrix, Ab(k-1,spin)%matrix, &
                           0.0_dp, b(k,spin)%matrix,error=error)
                   ELSE
                      CALL cp_fm_to_fm(Ab(k-1,spin)%matrix, b(k,spin)%matrix,error=error)
                   END IF
                END DO

             ELSE IF (mode == 'D') THEN

                iev = must_improve(i)
                ! create the new davidson vector
                DO spin=1, nspins

                   CALL cp_fm_set_all(R(spin)%matrix, 0.0_dp, error=error)
                   DO j=1, k-i
                      CALL cp_fm_to_fm(Ab(j,spin)%matrix, X(spin)%matrix,error=error)
                      CALL cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                               -evals(iev,iter-1), Sb(j,spin)%matrix,error=error)
                      CALL cp_fm_get_element(Us%matrix, j, iev, tmp)
                      CALL cp_fm_scale_and_add(1.0_dp, R(spin)%matrix, &
                                               tmp, X(spin)%matrix,error=error)
                   END DO

                   IF (tddfpt_control%invert_S) THEN
                      CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                           1.0_dp, t_env%invS(spin)%matrix, R(spin)%matrix, &
                           0.0_dp, X(spin)%matrix,error=error)
                   ELSE
                      CALL cp_fm_to_fm(R(spin)%matrix, X(spin)%matrix,error=error)
                   END IF

                   !----------------!
                   ! preconditioner !
                   !----------------!
                   IF (qs_env%dft_control%tddfpt_control%precond) THEN
                      DO col=1, p_env%n_mo(spin)
                         IF (col<=n_ev) THEN
                            tmp2 = ABS(evals(iev,iter-1) - in_evals(col))
                         ELSE
                            tmp2 = ABS(evals(iev,iter-1) - (in_evals(n_ev)+10.0_dp))
                         END IF
                         ! protect against division by 0 by a introducing a cutoff.
                         tmp2=MAX(tmp2,100*EPSILON(1.0_dp))
                         DO row=1, p_env%n_ao(spin)
                            CALL cp_fm_get_element(X(spin)%matrix, row, col, tmp)
                            CALL cp_fm_set_element(b(k,spin)%matrix, row, col, tmp/tmp2,error=error)
                         END DO
                      END DO
                   ELSE
                      CALL cp_fm_to_fm(X(spin)%matrix, b(k,spin)%matrix,error=error)
                   END IF

                END DO

             ELSE
                IF (output_unit>0) WRITE (output_unit,*) "unknown mode"
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF

          END IF

          CALL p_preortho(p_env, qs_env, b(k,:),error=error)
          DO j=1, tddfpt_control%n_reortho
             CALL reorthogonalize(b(k,:), b, Sb, R, k-1,error=error) ! R is temp
          END DO
          CALL normalize(b(k,:), R, qs_env%matrix_s,error=error) ! R is temp
          DO spin=1, nspins
             CALL cp_fm_to_fm(b(k,spin)%matrix, X(spin)%matrix,error=error)
          END DO
          CALL apply_op(X, Ab(k,:),  p_env, qs_env, &
               qs_env%dft_control%tddfpt_control%do_kernel,error=error)
          CALL p_postortho(p_env, qs_env, Ab(k,:),error=error)
          DO spin=1, nspins
             CALL cp_dbcsr_sm_fm_multiply(qs_env%matrix_s(1)%matrix, &
                  b(k,spin)%matrix, &
                  Sb(k,spin)%matrix, &
                  p_env%n_mo(spin),error=error)
          END DO
       END DO

       !--------------------------------------------!
       ! deallocate memory for the reduced matrices !
       !--------------------------------------------!
       IF (ASSOCIATED(Atilde%matrix)) CALL cp_fm_release(Atilde%matrix,error=error)
       IF (ASSOCIATED(Us%matrix)) CALL cp_fm_release(Us%matrix,error=error)
       IF (ASSOCIATED(tilde_fm_struct)) CALL cp_fm_struct_release(tilde_fm_struct,error=error)

       !------------------------------------------!
       ! allocate memory for the reduced matrices !
       !------------------------------------------!
       CALL cp_fm_struct_create(tilde_fm_struct, &
            qs_env%para_env, &
            qs_env%blacs_env, &
            k, k,error=error)
       CALL cp_fm_create(Atilde%matrix, &
            tilde_fm_struct, &
            routineP//"Atilde",error=error)
       CALL cp_fm_create(Us%matrix, &
            tilde_fm_struct, &
            routineP//"Us",error=error)

       !---------------------------------------!
       ! calc the matrix Atilde = transp(b)*Ab !
       !---------------------------------------!
       DO i=1, k
          DO j=1, k
             Atilde_ij= 0.0_dp
             DO spin=1, nspins
                CALL cp_fm_trace(b(i,spin)%matrix, Ab(j,spin)%matrix, tmp,error=error)
                Atilde_ij = Atilde_ij + tmp
             END DO
             CALL cp_fm_set_element(Atilde%matrix, i, j, Atilde_ij,error=error)
          END DO
       END DO

       !--------------------!
       ! diagonalize Atilde !
       !--------------------!
       evals_tmp(:) = evals(:,iter)
       CALL cp_fm_syevd(Atilde%matrix, Us%matrix, evals_tmp(:),error=error)
       evals(:,iter) =  evals_tmp(:)

       !-------------------!
       ! check convergence !
       !-------------------!
       evals_difference =  1.0_dp
       IF (iter/=1) THEN

          evals_difference(:) = ABS((evals(1:n_ev,iter-1)-evals(1:n_ev,iter)))
          ! For debugging
          IF (output_unit>0) THEN
             WRITE (output_unit,*)
             DO i=1, n_ev
                WRITE (output_unit,'(2X,F10.7,T69,ES11.4)') evals(i,iter)*evolt, evals_difference(i)
             END DO
             WRITE (output_unit,*)
          END IF

          convergence = MAXVAL(evals_difference)
          IF (output_unit>0) WRITE (output_unit, '(2X,I4,T69,ES11.4)') k, convergence

          IF (convergence < tddfpt_control%tolerance) THEN
             res = .TRUE.
             EXIT iteration
          END IF
       END IF

       IF (mode == 'L') THEN
          n_kv = 1
       ELSE
          must_improve = 0
          DO i=1, n_ev
             IF (evals_difference(i) > tddfpt_control%tolerance) must_improve(i) = 1
          END DO
!! Set must_improve to 1 if all the vectors should
!! be updated in one iteration.
!!          must_improve = 1
          n_kv = SUM(must_improve)
          j = 1
          DO i=1, n_ev
             IF (must_improve(i) == 1) THEN
                must_improve(j) = i
                j = j + 1
             END IF
          END DO
       END IF

       IF (k+n_kv>max_kv) EXIT iteration

       iter = iter + 1

    END DO iteration

    IF (PRESENT(out_evals)) THEN
       out_evals(1:n_ev) = evals(1:n_ev,iter)
    END IF

    DO spin=1, nspins
       DO j=1, n_ev
          CALL cp_fm_set_all(t_env%evecs(j,spin)%matrix, 0.0_dp,error=error)
          DO i=1, k
             CALL cp_fm_get_element(Us%matrix, i, j, tmp)
             CALL cp_fm_scale_and_add(1.0_dp, t_env%evecs(j,spin)%matrix, &
                                      tmp, b(i,spin)%matrix,error=error)
          END DO
       END DO
    END DO

    DO spin=1, nspins
       DO i=1, max_kv
          IF (ASSOCIATED(b(i,spin)%matrix)) &
               CALL cp_fm_release(b(i,spin)%matrix,error=error)
          IF (ASSOCIATED(Ab(i,spin)%matrix)) &
               CALL cp_fm_release(Ab(i,spin)%matrix,error=error)
          IF (ASSOCIATED(Sb(i,spin)%matrix)) &
               CALL cp_fm_release(Sb(i,spin)%matrix,error=error)
       END DO
    END DO

    !----------!
    ! clean up !
    !----------!
    IF (ASSOCIATED(Atilde%matrix)) CALL cp_fm_release(Atilde%matrix,error=error)
    IF (ASSOCIATED(Us%matrix)) CALL cp_fm_release(Us%matrix,error=error)
    IF (ASSOCIATED(tilde_fm_struct)) CALL cp_fm_struct_release(tilde_fm_struct,error)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, X,error=error)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, R,error=error)
    DO spin=1, nspins
       CALL cp_fm_struct_release(kv_fm_struct(spin)%struct,error=error)
       DO i=1, max_kv
          IF (ASSOCIATED(b(i,spin)%matrix)) &
               CALL cp_fm_release(b(i,spin)%matrix,error=error)
          IF (ASSOCIATED(Ab(i,spin)%matrix)) &
               CALL cp_fm_release(Ab(i,spin)%matrix,error=error)
          IF (ASSOCIATED(Sb(i,spin)%matrix)) &
               CALL cp_fm_release(Sb(i,spin)%matrix,error=error)
       END DO
    END DO
    DEALLOCATE(b, Ab, Sb, evals, evals_tmp, evals_difference, must_improve, kv_fm_struct, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(X, R, b, Ab, Sb, kv_fm_struct)

    CALL timestop(handle)

  END FUNCTION iterative_solver

  ! X        : the vector on which to apply the op
  ! R        : the result
  ! t_env    : td-dft environment (mainly control information)
  ! p_env    : perturbation environment (variables)
  !            both of these carry info for the tddfpt calculation
  ! qs_env   : info about a quickstep ground state calculation

! *****************************************************************************
  SUBROUTINE apply_op(X, R,  p_env, qs_env, do_kernel, error)

    TYPE(cp_fm_p_type), DIMENSION(:)         :: X, R
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: do_kernel
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_op', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspins, spin
    INTEGER, SAVE                            :: counter = 0

    CALL timeset(routineN,handle)

    counter = counter + 1

    nspins = qs_env%dft_control%nspins

    !------------!
    ! R = HX-SXL !
    !------------!
    CALL p_op_l1(p_env, qs_env, X, R,error=error)  ! acts on both spins, result in R

    !-----------------!
    ! calc P1 and     !
    ! R = R + K(P1)*C !
    !-----------------!
    IF (do_kernel) THEN
       DO spin=1, nspins
          CALL cp_dbcsr_set(p_env%p1(spin)%matrix, 0.0_dp, error=error) ! optimize?
          CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=p_env%psi0d(spin)%matrix,&
                                  matrix_g=X(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=X(spin)%matrix,&
                                  matrix_g=p_env%psi0d(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dp,error=error)
       END DO
       DO spin=1, nspins
          CALL cp_fm_set_all(X(spin)%matrix, 0.0_dp,error=error)
       END DO
       CALL p_op_l2(p_env, qs_env, p_env%p1, X, &
                    alpha=1.0_dp, beta=0.0_dp,error=error) ! X = beta*X + alpha*K(P1)*C
       DO spin=1, nspins
          CALL cp_fm_scale_and_add(1.0_dp, R(spin)%matrix, &
                                   1.0_dp, X(spin)%matrix,error=error) ! add X to R
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE apply_op

! *****************************************************************************
  SUBROUTINE allocate_krylov_vectors(vectors, vectors_name, &
                                     startv, n_v, nspins, fm_struct, error)

    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: vectors
    CHARACTER(LEN=*), INTENT(IN)             :: vectors_name
    INTEGER, INTENT(IN)                      :: startv, n_v, nspins
    TYPE(cp_fm_struct_p_type), &
      DIMENSION(:), POINTER                  :: fm_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_krylov_vectors', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: mat_name
    INTEGER                                  :: index, spin

    DO spin=1, nspins
       DO index=startv, startv + n_v - 1
          NULLIFY(vectors(index,spin)%matrix)
          mat_name = routineP//vectors_name//TRIM(cp_to_string(index))&
                 //","//TRIM(cp_to_string(spin))
          CALL cp_fm_create(vectors(index,spin)%matrix,  &
               fm_struct(spin)%struct, mat_name,error=error)
          IF (.NOT.ASSOCIATED(vectors(index,spin)%matrix)) &
               CALL stop_program(routineN,moduleN,__LINE__,&
                                 "Could not allocate "//TRIM(mat_name)//".")
       END DO
    END DO

  END SUBROUTINE allocate_krylov_vectors

END MODULE qs_tddfpt_eigensolver

!-----------------------------------------------------------------------------!
!   CP2K: A GENERAL PROGRAM TO PERFORM MOLECULAR DYNAMICS SIMULATIONS         !
!   COPYRIGHT (C) 2000  CP2K DEVELOPERS GROUP                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****M* CP2K/qs_tddfpt_eigensolver *
!!
!!   NAME
!!     qs_tddfpt_eigensolver
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     
!!     
!!****
MODULE qs_tddfpt_eigensolver
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_symm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release, &
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_get_element,&
                                             cp_fm_set_element,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE cp_log_handling,                 ONLY: cp_to_string
  USE qs_mo_types,                     ONLY: get_mo_set
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_p_env_methods,                ONLY: p_op_l1,&
                                             p_op_l2,&
                                             p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_tddfpt_types
  USE qs_tddfpt_utils,                 ONLY: calc_initial_guess,&
                                             calc_norm,&
                                             davidson_initial_guess, &
                                             normalize, &
                                             reorthogonalize, &
                                             dot_prod
  USE sparse_matrix_types,             ONLY: set_matrix, &
                                             get_matrix_info, &
                                             get_matrix_diagonal
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  TYPE cp_fm_struct_p_type
     TYPE(cp_fm_struct_type), POINTER :: struct
  END TYPE cp_fm_struct_p_type

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_eigensolver"
  INTEGER, PARAMETER          :: MAX_ITER_GLOBAL = 20
  REAL(kind=dbl), PARAMETER   :: TOLERANCE = 1.0E-6

  PRIVATE

  PUBLIC :: eigensolver

CONTAINS

! *****************************************************************************

  SUBROUTINE eigensolver(t_env, p_env, qs_env, glob_env)

    TYPE(tddfpt_env_type), INTENT(INOUT)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), INTENT(IN)    :: qs_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: glob_env

    CHARACTER(len=*), PARAMETER :: routineN = "davidson", &
      routineP = moduleN//"/"//routineN

    INTEGER :: col, handle, i, index, iter, j, max_iter, &
      max_krylovspace_dim, max_krylovvec, n_ev, n_kv, n_reortho, n_spins, &
      row, spin, stat
    INTEGER :: K_rows
    REAL(dbl)                                :: Atilde_ij, tmp, R_norm, convergence
    REAL(kind=dbl), DIMENSION(:), ALLOCATABLE :: K_diagonal, S_diagonal
    REAL(kind=dbl), DIMENSION(:), POINTER :: occ_evals
    REAL(kind=dbl), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evals
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: R, X
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: b, Ab, Sb
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_p_type), &
      DIMENSION(:), POINTER                  :: kv_fm_struct
    TYPE(cp_fm_p_type)                       :: Atilde, Us
    TYPE(cp_fm_struct_type), POINTER         :: tilde_fm_struct, Us_fm_struct

    CALL timeset(routineN,"I","",handle)
    NULLIFY(ao_mo_fm_pools, X, R, b, Ab, Sb )
    n_ev      = qs_env%dft_control%tddfpt_control%n_ev
    n_spins   = qs_env%dft_control%nspins
    n_reortho = qs_env%dft_control%tddfpt_control%n_reortho
    iter = 1
    DO spin=1, n_spins
       CALL get_matrix_info(qs_env%matrix_ks(spin)%matrix, nrow=K_rows)
       ALLOCATE(K_diagonal(K_rows), S_diagonal(K_rows))
       CALL get_matrix_diagonal(qs_env%matrix_ks(spin)%matrix, K_diagonal)
       CALL get_matrix_diagonal(qs_env%matrix_s(spin)%matrix, S_diagonal)
       CALL get_mo_set(qs_env%mos(spin)%mo_set, eigenvalues=occ_evals)
    END DO


    !-----------------------------------------!
    ! determine the size of the problem       !
    ! and how many krylov space vetors to use ! 
    !-----------------------------------------!
    max_krylovspace_dim = SUM(p_env%n_ao(1:n_spins)*p_env%n_mo(1:n_spins))
    n_kv = MAX(1, n_ev)
    max_iter = MIN(( max_krylovspace_dim / n_kv ) * n_kv, MAX_ITER_GLOBAL)
    max_krylovvec = max_iter * n_kv

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, X, name=routineP//":X")
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, R, name=routineP//":R")

    ALLOCATE(evals(max_krylovvec, 0:max_iter))
    ALLOCATE(b(max_krylovvec,n_spins), Ab(max_krylovvec,n_spins), &
             Sb(max_krylovvec,n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of kv/opkv failed", glob_env)
    ALLOCATE(kv_fm_struct(n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of kv_fm_struct failed", &
                                     glob_env)
    DO spin=1, n_spins
       NULLIFY (kv_fm_struct(spin)%struct)
       CALL cp_fm_struct_create(kv_fm_struct(spin)%struct, qs_env%para_env, &
                                qs_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_mo(spin))
    END DO
    DO spin=1, n_spins
       DO i=1, max_krylovvec
          NULLIFY(b(i,spin)%matrix, Ab(i,spin)%matrix, Sb(i,spin)%matrix)
       END DO
    END DO
    CALL allocate_krylov_vectors(b, "b-", iter, n_kv, n_spins, kv_fm_struct)
    CALL allocate_krylov_vectors(Ab, "Ab-", iter, n_kv, n_spins, kv_fm_struct)
    CALL allocate_krylov_vectors(Sb, "Sb-", iter, n_kv, n_spins, kv_fm_struct)


    !-----------------------!
    ! provide initial guess !
    !-----------------------!
    CALL davidson_initial_guess(b, n_kv, qs_env, p_env, evals)
    DO i=1, n_kv
       DO j=1, n_reortho
          CALL p_preortho(p_env, qs_env, b(i,:))
          CALL reorthogonalize(b(i,:), b, Sb, i-1)
       END DO
       CALL normalize(b(i,:), R, qs_env%matrix_s, glob_env)
       DO spin=1, n_spins
          CALL cp_fm_to_fm(b(i,spin)%matrix, X(spin)%matrix)
       END DO
       CALL apply_op(X, Ab(i,:), t_env, p_env, qs_env)
       DO spin=1, n_spins
          CALL cp_sm_fm_multiply(qs_env%matrix_s(spin)%matrix, &
                                 b(i,spin)%matrix, &
                                 Sb(i,spin)%matrix, &
                                 p_env%n_mo(spin), &
                                 qs_env%para_env, 1.0_dbl, 0.0_dbl)
       END DO
    END DO

    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, &
              '(2X,A,2X,A,T69,A)') &
              "Step", "Vectors", "Convergence"
       WRITE (glob_env%scr, '(2X,A)') &
              "-----------------------------------------------------------------------------"
    END IF

    iteration: DO 
       
       !------------------------------------------!
       ! allocate memory for the reduced matrices !
       !------------------------------------------!
       CALL cp_fm_struct_create(tilde_fm_struct, &
                                qs_env%para_env, &
                                qs_env%blacs_env, &
                                iter*n_kv, iter*n_kv)
       CALL cp_fm_create(Atilde%matrix, &
                         tilde_fm_struct, &
                         routineP//"Atilde")
       CALL cp_fm_create(Us%matrix, &
                         tilde_fm_struct, &
                         routineP//"Us")

       !---------------------------------------!
       ! calc the matrix Atilde = transp(b)*Ab !
       !---------------------------------------!
       DO i=1, iter*n_kv
          DO j=1, iter*n_kv
             Atilde_ij= 0.0_dbl
             DO spin=1, n_spins
                CALL cp_fm_trace(b(i,spin)%matrix, Ab(j,spin)%matrix, tmp)
                Atilde_ij = Atilde_ij + tmp
             END DO
             CALL cp_fm_set_element(Atilde%matrix, i, j, Atilde_ij)
          END DO
       END DO
      
       !--------------------!
       ! diagonalize Atilde !
       !--------------------!
       CALL cp_fm_syevd(Atilde%matrix, Us%matrix, evals(:,iter))
       DO i=1, n_ev
          IF (glob_env%ionode) WRITE (*,*) evals(i,iter)
       END DO

       !-------------------!
       ! check convergence !
       !-------------------!
       convergence = ABS(evals(n_ev,iter)-evals(n_ev,iter-1))
       IF (glob_env%ionode) &
            WRITE (glob_env%scr, &
            '(2X,I4,2X,I7,T69,ES11.4)') &
            iter, iter*n_kv, convergence
       IF (convergence<TOLERANCE) EXIT iteration
       IF (iter == max_iter) EXIT iteration
       iter = iter + 1

       !-----------------!
       ! add new vectors !
       !-----------------!
       CALL allocate_krylov_vectors(b, "b-", iter, n_kv, n_spins, kv_fm_struct)
       CALL allocate_krylov_vectors(Ab, "Ab-", iter, n_kv, n_spins, kv_fm_struct)
       CALL allocate_krylov_vectors(Sb, "Sb-", iter, n_kv, n_spins, kv_fm_struct)
       DO i=(iter-1)*n_kv+1, iter*n_kv
          IF (qs_env%dft_control%tddfpt_control%diag_method .EQ. "DAVIDSON") THEN
             DO spin=1, n_spins
                CALL cp_fm_set_all(b(i,spin)%matrix, 0.0_dbl)
                DO j=1, n_kv
                   CALL cp_fm_get_element(Us%matrix, i-(iter-1)*n_kv, j, tmp)
                   CALL cp_fm_scale_and_add(1.0_dbl, b(i,spin)%matrix, &
                                            tmp, Ab(j,spin)%matrix)
                   CALL cp_fm_scale_and_add(1.0_dbl, b(i,spin)%matrix, &
                                            -1.0_dbl*evals(j,iter)*tmp, Sb(j,spin)%matrix)
                END DO
                DO row=1, p_env%n_ao(spin)
                   DO col=1, p_env%n_mo(spin)
                      CALL cp_fm_get_element(b(i,spin)%matrix, row, col, tmp)
                      tmp = tmp / (K_diagonal(row) - (evals(i-(iter-1)*n_kv,iter-1)+occ_evals(col)))
                      CALL cp_fm_set_element(b(i,spin)%matrix, row, col, tmp)
                   END DO
                END DO
                IF (qs_env%dft_control%tddfpt_control%invert_S) THEN
                   CALL cp_fm_to_fm(b(i,spin)%matrix, R(spin)%matrix)
                   CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                        1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
                        0.0_dbl, b(i,spin)%matrix)
                END IF
             END DO
          ELSE IF (qs_env%dft_control%tddfpt_control%diag_method .EQ. "LANCZOS") THEN
             DO spin=1, n_spins
                IF (qs_env%dft_control%tddfpt_control%invert_S) THEN
                   CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                        1.0_dbl, t_env%invS(spin)%matrix, Ab(i-n_kv,spin)%matrix, &
                        0.0_dbl, b(i,spin)%matrix)
                ELSE
                   CALL cp_fm_to_fm(Ab(i-n_kv,spin)%matrix, b(i,spin)%matrix)
                END IF
             END DO
          ELSE
             CALL stop_program(routineP, "Unknown method.");
          END IF
          DO j=1, n_reortho
             CALL p_preortho(p_env, qs_env, b(i,:))
             CALL reorthogonalize(b(i,:), b, Sb, i-1)
          END DO
          CALL normalize(b(i,:), R, qs_env%matrix_s, glob_env)
          DO spin=1, n_spins
             CALL cp_fm_to_fm(b(i,spin)%matrix, X(spin)%matrix)
          END DO
          CALL apply_op(X, Ab(i,:), t_env, p_env, qs_env)
          DO spin=1, n_spins
             CALL cp_sm_fm_multiply(qs_env%matrix_s(spin)%matrix, &
                  b(i,spin)%matrix, &
                  Sb(i,spin)%matrix, &
                  p_env%n_mo(spin), &
                  qs_env%para_env, 1.0_dbl, 0.0_dbl)
          END DO
       END DO

       !--------------------------------------------!
       ! deallocate memory for the reduced matrices !
       !--------------------------------------------!
       IF (ASSOCIATED(Atilde%matrix)) CALL cp_fm_release(Atilde%matrix)
       IF (ASSOCIATED(Us%matrix)) CALL cp_fm_release(Us%matrix)
       IF (ASSOCIATED(tilde_fm_struct)) CALL cp_fm_struct_release(tilde_fm_struct)


    END DO iteration

    !-------------------!
    ! Print the results !
    !-------------------!
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, *)
       WRITE (glob_env%scr, *)
       WRITE (glob_env%scr, '(1X,A,T19,A,T45,A)') &
            "Nr.", "Eigenvalue (hartree)", "Eigenvalue (eV)"
       DO i=1, MIN(n_ev,max_krylovvec)
          WRITE (glob_env%scr, '(1X,I3,T20,F19.10,2X,F19.10)') &
                                i, evals(i, iter), &
                                evals(i, iter)*evolt
       END DO
       WRITE (glob_env%scr, *)
    END IF

    !----------!
    ! clean up !
    !----------!
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, X)
    CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, R)
    DO spin=1, n_spins
       DO i=1, max_krylovvec
          IF (ASSOCIATED(b(i,spin)%matrix)) &
             CALL cp_fm_release(b(i,spin)%matrix)
          IF (ASSOCIATED(Ab(i,spin)%matrix)) &
             CALL cp_fm_release(Ab(i,spin)%matrix)
          IF (ASSOCIATED(Sb(i,spin)%matrix)) &
             CALL cp_fm_release(Sb(i,spin)%matrix)
       END DO
    END DO
    DEALLOCATE(b, Ab, Sb, evals, K_diagonal, S_diagonal)
    NULLIFY(X, R, b, Ab, Sb)
    
    CALL timestop(0.0_dbl, handle)

  END SUBROUTINE eigensolver

! *****************************************************************************

  ! X        : the vector on which to apply the op
  ! R        : the result
  ! t_env    : td-dft environment (mainly control information)
  ! p_env    : perturbation environment (variables)
  !            both of these carry info for the tddfpt calculation
  ! qs_env   : info about a quickstep ground state calculation

  SUBROUTINE apply_op(X, R, t_env, p_env, qs_env)


    TYPE(cp_fm_p_type), DIMENSION(:)         :: X, R
    TYPE(tddfpt_env_type), INTENT(INOUT)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type)                :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "apply_op", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: handle, n_spins, spin
    INTEGER, SAVE                            :: counter = 0

    CALL timeset(routineN,"I","",handle)

    counter = counter + 1
    
    n_spins = qs_env%dft_control%nspins

    !------------!
    ! R = HX-SXL !
    !------------!
    CALL p_op_l1(p_env, qs_env, X, R)  ! acts on both spins, result in R

    !-----------------!
    ! calc P1 and     !
    ! R = R + K(P1)*C !
    !-----------------!
    IF (qs_env%dft_control%tddfpt_control%do_kernel) THEN
       DO spin=1, n_spins
          CALL set_matrix(p_env%p1(spin)%matrix, value=0.0_dbl) ! optimize?
          CALL cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=p_env%psi0d(spin)%matrix,&
                                  matrix_g=X(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dbl) 
          CALL cp_sm_plus_fm_fm_t(p_env%p1(spin)%matrix,&
                                  matrix_v=X(spin)%matrix,&
                                  matrix_g=p_env%psi0d(spin)%matrix,&
                                  ncol=p_env%n_mo(spin), &
                                  alpha=0.5_dbl) 
       END DO
       CALL p_op_l2(p_env, qs_env, p_env%p1, R, &
                    alpha=1.0_dbl, beta=1.0_dbl) ! R = beta*R + alpha*K(P1)*C
    END IF

    CALL timestop(0.0_dbl, handle)

  END SUBROUTINE apply_op

! *****************************************************************************

  SUBROUTINE allocate_krylov_vectors(vectors, vectors_name, &
                                     iteration, n_v, n_spins, fm_struct)

    TYPE(cp_fm_p_type), DIMENSION(:,:), &
      POINTER                                :: vectors
    CHARACTER(LEN=*), INTENT(IN)             :: vectors_name
    INTEGER, INTENT(IN)                      :: iteration, n_v, n_spins
    TYPE(cp_fm_struct_p_type), &
      DIMENSION(:), POINTER                  :: fm_struct

    CHARACTER(LEN=*), PARAMETER :: routineN = "allocate_krylov_vectors", &
                                   routineP = moduleN//"/"//routineN

    CHARACTER(LEN=80)                        :: mat_name
    INTEGER                                  :: i, index, spin

    DO spin=1, n_spins
       DO i=1, n_v
          index = (iteration-1)*n_v+i
          NULLIFY(vectors(index,spin)%matrix)
          mat_name = routineP//vectors_name//TRIM(cp_to_string(index))&
                 //","//TRIM(cp_to_string(spin))
          CALL cp_fm_create(vectors(index,spin)%matrix,  &
               fm_struct(spin)%struct, mat_name)
          IF (.NOT.ASSOCIATED(vectors(index,spin)%matrix)) &
               CALL stop_program(routineP, &
               "Could not allocate "//TRIM(mat_name)//".")
       END DO
    END DO

  END SUBROUTINE allocate_krylov_vectors

! *****************************************************************************

END MODULE qs_tddfpt_eigensolver

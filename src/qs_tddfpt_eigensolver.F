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
                                             cp_fm_trace
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm
  USE cp_log_handling,                 ONLY: cp_to_string
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_p_env_methods,                ONLY: p_op_l1,&
                                             p_op_l2,&
                                             p_postortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_tddfpt_types
  USE qs_tddfpt_utils,                 ONLY: calc_initial_guess,&
                                             calc_norm,&
                                             reorthogonalize
  USE sparse_matrix_types,             ONLY: set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  TYPE cp_fm_struct_p_type
     TYPE(cp_fm_struct_type), POINTER :: struct
  END TYPE cp_fm_struct_p_type

  CHARACTER(len=*), PARAMETER :: moduleN = "qs_tddfpt_eigensolver"
  INTEGER, PARAMETER          :: MAX_ITER_GLOBAL = 20
  REAL(kind=dbl), PARAMETER   :: TOLERANCE = 1.0E-5

  PRIVATE

  PUBLIC :: tddfpt_es

CONTAINS

! *****************************************************************************

  SUBROUTINE tddfpt_es(t_env, p_env, qs_env, glob_env)
    

    TYPE(tddfpt_env_type), INTENT(INOUT)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), INTENT(IN)    :: qs_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: glob_env

    CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_es_lanczos", &
      routineP = moduleN//"/"//routineN

    INTEGER :: col, handle, i, index, iteration, j, max_iter, &
      max_krylovspace_dim, max_krylovvec, n_ev, n_kv, n_reortho, n_spins, &
      row, spin, stat
    REAL(dbl)                                :: beta, difference, &
                                                dot_product, tmp, tol
    REAL(kind=dbl), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tddfpt_eigenvalues
    TYPE(cp_fm_p_type)                       :: T, Tev
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: R, X
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: kv, opkv, residues
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_p_type), &
      DIMENSION(:), POINTER                  :: kv_fm_struct
    TYPE(cp_fm_struct_type), POINTER         :: T_fm_struct, Tev_fm_struct

    CALL timeset(routineN,"I","",handle)

    NULLIFY(ao_mo_fm_pools, X, R, kv )
    n_ev      = qs_env%dft_control%tddfpt_control%n_ev
    n_spins   = qs_env%dft_control%nspins
    tol       = 1.0E-10_dbl
    n_reortho = qs_env%dft_control%tddfpt_control%n_reortho
    iteration = 1

    !-----------------------------------------!
    ! determine the size of the problem       !
    ! and how many krylov space vetors to use ! 
    !-----------------------------------------!
    max_krylovspace_dim = 0
    DO spin=1, n_spins
       max_krylovspace_dim = max_krylovspace_dim + p_env%n_ao(spin)*p_env%n_mo(spin)
    END DO
    n_kv = MAX(2, n_ev) ! how many vectors per iteration
    max_iter = MIN(( max_krylovspace_dim / n_kv ) * n_kv, MAX_ITER_GLOBAL)
    max_krylovvec = max_iter * n_kv

    !----------------------!
    ! allocate the vectors !
    !----------------------!
    ALLOCATE(tddfpt_eigenvalues(max_krylovvec, 0:max_iter))
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, X, &
                                 name=routineP//":X")
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, R, &
                                 name=routineP//":R")
    ALLOCATE(kv(max_krylovvec,n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of kv failed", glob_env)
    ALLOCATE(opkv(max_krylovvec,n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of opkv failed", glob_env)
    DO spin=1, n_spins
       DO i=1, max_krylovvec
          NULLIFY(kv(i,spin)%matrix, opkv(i,spin)%matrix)
       END DO
    END DO
    ALLOCATE(residues(n_ev,n_spins), stat=stat)
    IF (stat /= 0) CALL stop_program(routineP, &
                                     "allocation of residues failed", glob_env)
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
    CALL allocate_krylov_vectors(kv, "kv-", iteration, n_kv, n_spins, kv_fm_struct)
    CALL allocate_krylov_vectors(opkv, "opkv-", iteration, n_kv, n_spins, kv_fm_struct)
    CALL allocate_krylov_vectors(residues, "residues-", iteration, n_kv, n_spins, kv_fm_struct)


    !-----------------------!
    ! provide initial guess !
    !-----------------------!
    CALL calc_initial_guess(X, qs_env, p_env, tddfpt_eigenvalues)
    

    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, &
              '(2X,A,2X,A,T69,A)') &
              "Step", "Vectors", "Convergence"
       WRITE (glob_env%scr, '(2X,A)') &
              "-----------------------------------------------------------------------------"
    END IF

    restart: DO 

       krylov: DO j=1, n_kv
       
          index = (iteration-1)*n_kv+j
          
          !---------------------------------------!
          ! apply the operator to working vec and !
          ! store the result for later use (T)    !
          !---------------------------------------!
          IF (j /= 1) THEN
             CALL apply_op(X, R, t_env, p_env, qs_env) ! result in X
             DO spin=1, n_spins
                CALL cp_fm_to_fm(X(spin)%matrix, &
                                                opkv(index-1,spin)%matrix)
             END DO
          END IF

          DO i=1, n_reortho
             CALL reorthogonalize(X, kv, qs_env%s, R, index, glob_env)
          END DO

          CALL calc_norm(norm=beta, x=X, tmp_vec=R, &
                         qs_env=qs_env, glob_env=glob_env)

          !-------------!
          ! normalize X ! 
          !-------------!
          IF (ABS(beta) > tol) THEN
             DO spin=1, n_spins
                CALL cp_fm_scale(SQRT(1._dbl/beta), X(spin)%matrix)
             END DO
          ELSE 
             IF (index==1) THEN
                CALL stop_program(routineP, "Initial guess to close to zero.")
             ELSE
                ! already converged
                index = index-1
                max_krylovvec = index
                IF (glob_env%ionode) THEN
                   WRITE (*,*)
                   WRITE (*,*) "Converged with ", index," krylov vectors!"
                   WRITE (*,*)
                END IF
                EXIT krylov
             END IF
          END IF

          !---------------------------!
          ! write current vector back !
          !---------------------------!
          DO spin=1, n_spins
             CALL cp_fm_to_fm(X(spin)%matrix, &
                                             kv(index,spin)%matrix)
          END DO

       END DO krylov

       !---------------------------------------!
       ! compute the last of the Op*kv vectors !
       !---------------------------------------!
       DO spin=1, n_spins
          CALL cp_fm_to_fm(kv(index,spin)%matrix, &
                                          X(spin)%matrix)
       END DO
       CALL apply_op(X, R, t_env, p_env, qs_env) ! result in X
       DO spin=1, n_spins
          CALL cp_fm_to_fm(X(spin)%matrix, &
                                         opkv(index,spin)%matrix)
       END DO
       

       !--------------------------------------------!
       ! allocate memory for the T and Tev matrices !
       !--------------------------------------------!
       IF (ASSOCIATED(T%matrix)) THEN
          CALL cp_fm_release(T%matrix)
       END IF
       NULLIFY(T_fm_struct)
       CALL cp_fm_struct_create(T_fm_struct, &
                                qs_env%para_env, &
                                qs_env%blacs_env, &
                                index, index)
       CALL cp_fm_create(T%matrix, &
                         T_fm_struct, &
                         routineP//"T")
       IF (ASSOCIATED(Tev%matrix)) THEN
          CALL cp_fm_release(Tev%matrix)
       END IF
       NULLIFY(Tev_fm_struct)
       CALL cp_fm_struct_create(Tev_fm_struct, &
                                qs_env%para_env, &
                                qs_env%blacs_env, &
                                index, index)
       CALL cp_fm_create(Tev%matrix, &
                         Tev_fm_struct, &
                         routineP//"Tev")


       !-----------------------------------------!
       ! calc the matrix T = transp(kv)*S*OP(kv) !
       !-----------------------------------------!
       DO col=1, index
          DO spin=1, n_spins
             CALL cp_sm_fm_multiply(qs_env%S(spin)%matrix, &
                                    opkv(col,spin)%matrix, &
                                    R(spin)%matrix, &
                                    p_env%n_mo(spin), &
                                    qs_env%para_env, 1.0_dbl, 0.0_dbl)
          END DO
          DO row=1, index
             dot_product = 0.0_dbl
             DO spin=1, n_spins
                CALL cp_fm_trace(kv(row,spin)%matrix, R(spin)%matrix, tmp)
                dot_product = dot_product + tmp
             END DO
             CALL cp_fm_set_element(T%matrix, row, col, dot_product)
          END DO
       END DO


       !---------------------------!
       ! calc the eigenvalues of T !
       !---------------------------!
       CALL cp_fm_syevd(T%matrix, Tev%matrix, tddfpt_eigenvalues(:,iteration))

       IF (glob_env%ionode) THEN
          DO i=1, n_ev
             WRITE (*,*) tddfpt_eigenvalues(i,iteration)
          END DO
       END IF


       !---------------------!
       ! check, if converged !
       !---------------------!
       difference = ABS(tddfpt_eigenvalues(n_ev,iteration) - &
                        tddfpt_eigenvalues(n_ev,iteration-1))
       IF (glob_env%ionode) &
       WRITE (glob_env%scr, &
              '(2X,I4,2X,I7,T69,ES11.4)') &
              iteration, index, difference
       IF (difference < TOLERANCE) THEN
          EXIT restart
       END IF

       IF (index == max_krylovvec) THEN
          WRITE (*,*) "MAXIMUM NUMBER OF KRYLOV SPACE VECTORS REACHED"
          EXIT restart
       END IF
       
       iteration = iteration + 1

       !--------------------------------------!
       ! allocate matrizes for next iteration !
       !--------------------------------------!
       CALL allocate_krylov_vectors(kv, "kv-", iteration, n_kv, n_spins, kv_fm_struct)
       CALL allocate_krylov_vectors(opkv, "opkv-", iteration, n_kv, n_spins, kv_fm_struct)

    END DO restart


    !-------------------!
    ! Print the results !
    !-------------------!
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, *)
       WRITE (glob_env%scr, *)
       WRITE (glob_env%scr, '(1X,A,T19,A,T45,A)') &
            "Nr.", "Eigenvalue (hartree)", "Eigenvalue (eV)"
       DO i=1, MIN(n_ev,max_krylovvec)
          ! this should use the physcons library !
          WRITE (glob_env%scr, '(1X,I3,T20,F19.10,2X,F19.10)') &
                                i, tddfpt_eigenvalues(i, iteration), &
                                tddfpt_eigenvalues(i, iteration)*27.211_dbl  
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
          IF (ASSOCIATED(kv(i,spin)%matrix)) &
             CALL cp_fm_release(kv(i,spin)%matrix)
          IF (ASSOCIATED(opkv(i,spin)%matrix)) &
             CALL cp_fm_release(opkv(i,spin)%matrix)
       END DO
       DO i=1, n_ev
          IF (ASSOCIATED(residues(i,spin)%matrix)) &
               CALL cp_fm_release(residues(i,spin)%matrix)
       END DO
    END DO
    IF (ASSOCIATED(T%matrix)) CALL cp_fm_release(T%matrix)
    IF (ASSOCIATED(Tev%matrix)) CALL cp_fm_release(Tev%matrix)
    DEALLOCATE(kv, opkv, residues, tddfpt_eigenvalues)
    NULLIFY(X, R, kv, opkv, residues)
    
    CALL timestop(0.0_dbl, handle)

  END SUBROUTINE tddfpt_es

! *****************************************************************************

  ! X        : on input the vector on which to apply the op
  !            on output the result
  ! R        : a temporary vector
  ! t_env    : td-dft environment (mainly control information)
  ! p_env    : perturbation environment (variables)
  !            both of these carry info for the tddfpt calculation
  ! qs_env   : info about a quickstep ground state calculation

  SUBROUTINE apply_op(X, R, t_env, p_env, qs_env)


    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: X, R
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

    !----------------------!
    ! orthogonalize R to C !
    !----------------------!
    CALL p_postortho(p_env, qs_env, R)

    !----------!
    ! X=invS*R !
    !----------!
    DO spin=1, n_spins
       CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                       1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
                       0.0_dbl, X(spin)%matrix)
    END DO

#if 0
    CALL p_op_l1(p_env, qs_env, X, R) ! acts on both spins

    !----------!
    ! X=invS*R !
    !----------!
    DO spin=1, n_spins
       CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                       1.0_dbl, t_env%invS(spin)%matrix, R(spin)%matrix, &
                       0.0_dbl, X(spin)%matrix)
    END DO
#endif

    CALL timestop(0.0_dbl, handle)

  END SUBROUTINE apply_op

! *****************************************************************************

  SUBROUTINE allocate_krylov_vectors(vectors, vectors_name, &
                                     iteration, n_ev, n_spins, fm_struct)


    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: vectors
    CHARACTER(LEN=*), INTENT(IN)             :: vectors_name
    INTEGER, INTENT(IN)                      :: iteration, n_ev, n_spins
    TYPE(cp_fm_struct_p_type), &
      DIMENSION(:), POINTER                  :: fm_struct

    CHARACTER(LEN=*), PARAMETER :: routineN = "allocate_krylov_vectors", &
      routineP = moduleN//"/"//routineN

    CHARACTER(LEN=80)                        :: mat_name
    INTEGER                                  :: i, index, spin

    DO spin=1, n_spins
       DO i=1, n_ev
          index = (iteration-1)*n_ev+i
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

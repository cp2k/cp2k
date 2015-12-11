!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Common framework for using eigenvectors of a Fock matrix as PAO basis.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_param_fock
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_create, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_reserve_diag_blocks, cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE pao_types,                       ONLY: pao_env_type
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_fock'

  PUBLIC :: pao_param_init_fock, pao_param_finalize_fock, pao_calc_U_fock

CONTAINS

! *****************************************************************************
!> \brief Initialize fock parametrization
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_param_init_fock(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_fock', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set)

    IF(dft_control%nspins /= 1) CPABORT("open shell not yet implemented")

    ! fill U0 with random noise
    CALL cp_dbcsr_init(pao%matrix_U0)
    CALL cp_dbcsr_create(pao%matrix_U0, template=matrix_s(1)%matrix, matrix_type="N")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_U0)
    CALL fill_randomly(pao%matrix_U0, pao%random_seed)

    CALL timestop(handle)
  END SUBROUTINE pao_param_init_fock


! *****************************************************************************
!> \brief Finalize fock parametrization
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_param_finalize_fock(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize_fock', &
      routineP = moduleN//':'//routineN

    CALL cp_dbcsr_release(pao%matrix_U0)

  END SUBROUTINE pao_param_finalize_fock


! *****************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param iatom ...
!> \param V ...
!> \param U ...
!> \param M1 ...
!> \param G ...
! *****************************************************************************
  SUBROUTINE pao_calc_U_fock(pao, iatom, V, U, M1, G)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER                                  :: iatom
    REAL(dp), DIMENSION(:, :), POINTER       :: V, U
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: M1, G

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_fock', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, m, n
    LOGICAL                                  :: found
    REAL(dp)                                 :: denom, invsqrt, nom
    REAL(dp), DIMENSION(:), POINTER          :: V_evals, W_evals
    REAL(dp), DIMENSION(:, :), POINTER       :: C, D4, D8, M2, M3, M4, M5, &
                                                M6, M7, P, U0, V_evecs, W, &
                                                W_evecs, Y, Z

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U0, row=iatom, col=iatom, block=U0, found=found)
    CPASSERT(ASSOCIATED(U0))
    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=Y, found=found)
    CPASSERT(ASSOCIATED(Y))

    IF(MAXVAL(ABS(V-TRANSPOSE(V)))>1e-14_dp) CPABORT("Expect symmetric matrix")
    n = SIZE(V, 1)
    m = SIZE(Y, 2)
    IF(SIZE(Y, 1) /= n) CPABORT("strange size of block Y")

    ! calculate P = sign(V)
    ALLOCATE(V_evals(n), V_evecs(n,n), P(n,n))
    V_evecs = V
    CALL diamat_all(V_evecs, V_evals)
    P(:,:) = 0.0_dp
    DO k=1, m ! First m states have occupation 1, remaining have 0 .
       DO i=1, n
          DO j=1, n
             P(i,j) = P(i,j) + V_evecs(i,k)*V_evecs(j,k)
          ENDDO
       ENDDO
    ENDDO
    IF(MAXVAL(ABS(P-TRANSPOSE(P)))>1e-14_dp) CPABORT("block P not symmetric")

    ! Apply projector P onto random matrix U0
    ALLOCATE(C(n,m))
    C = MATMUL(MATMUL(P, U0), Y)

    ! calc W = C^T * C
    ALLOCATE(W(m,m))
    W = MATMUL(TRANSPOSE(C), C)

    ! calc Z = 1 / sqrt(CC)
    ALLOCATE(Z(m,m), W_evals(m), W_evecs(m,m))
    W_evecs = W
    CALL diamat_all(W_evecs, W_evals)
    Z=0.0_dp
    DO k=1, m
       IF(W_evals(k) <= 0.0_dp) CYCLE !TODO: is this a clean solution?
       invsqrt = 1.0_dp / SQRT(W_evals(k))
       DO i=1, m
          DO j=1, m
             Z(i,j) = Z(i,j) + invsqrt * W_evecs(i,k)*W_evecs(j,k)
          ENDDO
       ENDDO
    ENDDO

    ! pice together block U
    U = MATMUL(MATMUL(C, Z), TRANSPOSE(Y))

    ! TURNING POINT (if calc grad) ------------------------------------------
    IF(PRESENT(G))THEN
       CPASSERT(PRESENT(M1))

       ! dU / dZ
       ALLOCATE(M2(m,m))
       M2 = MATMUL(MATMUL(TRANSPOSE(C), M1), Y)

       ! dU / dC
       ALLOCATE(M3(n,m))
       M3 = MATMUL(MATMUL(M1, Y), TRANSPOSE(Z))

       ! dZ / dW
       ALLOCATE(D4(m,m), M4(m,m))
       DO i=1, m
          IF(W_evals(i)<=0.0_dp) CYCLE
          DO j=1, m
             IF(W_evals(j)<=0.0_dp) CYCLE
             IF(i==j) THEN
                 D4(i,i) = -0.5_dp / (W_evals(i)*SQRT(W_evals(i)))
             ELSE
                 nom = 1.0_dp/SQRT(W_evals(i)) - 1.0_dp/SQRT(W_evals(j))
                 denom = W_evals(i) - W_evals(j)
                 D4(i,j) = nom / denom
             ENDIF
          ENDDO
       ENDDO
       M4 = fold_derivatives(M2, D4, W_evecs, +1)

       ! dW / dC
       ALLOCATE(M5(n,m))
       M5 = MATMUL(C, TRANSPOSE(M4) + M4)

       ALLOCATE(M6(n,m))
       M6 = M3 + M5

       ! dC / dP
       ALLOCATE(M7(n,n))
       M7 = MATMUL(MATMUL(M6, TRANSPOSE(Y)), TRANSPOSE(U0))

       ! dP / dV
       ALLOCATE(D8(n,n))
       DO i=1, n
          DO j=1, n
             IF(i==j) THEN
                  D8(i,i) = 0.0 ! derivative of step-function
             ELSE
                 nom = 0.0_dp
                 IF(i<=m) nom = nom + 1.0_dp
                 IF(j<=m) nom = nom - 1.0_dp
                 denom = V_evals(i) - V_evals(j)
                 IF(nom==0.0_dp) THEN ! avoid division by zero if denom==0 too.
                    D8(i,j) = 0.0_dp
                 ELSE
                    D8(i,j) = nom / denom
                 ENDIF
             ENDIF
          ENDDO
       ENDDO

       ! the final gradient
       G = fold_derivatives(M7, D8, V_evecs, +1)

       DEALLOCATE(M2, M3, M4, M5, M6, M7)
       DEALLOCATE(D4, D8)
    ENDIF

    DEALLOCATE(V_evals, V_evecs, W_evals, W_evecs)
    DEALLOCATE(Z, P, C, W)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_U_fock


! *****************************************************************************
!> \brief Helper routine for calculating derivatives
!> \param M ...
!> \param D ...
!> \param R ...
!> \param sym ...
!> \retval G ...
! *****************************************************************************
  FUNCTION fold_derivatives(M, D, R, sym) RESULT(G)
    REAL(dp), DIMENSION(:, :)                :: M, D, R
    INTEGER                                  :: sym
    REAL(dp), &
      DIMENSION(SIZE(M, 1), SIZE(M, 1))      :: G

    CHARACTER(len=*), PARAMETER :: routineN = 'fold_derivatives', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n
    REAL(dp), DIMENSION(:, :), POINTER       :: F, RF, RFR, RM, RMR

     IF(ABS(sym) /= 1) CPABORT("expected +1 or -1 as sym")

     n = SIZE(R, 1)

     ALLOCATE(RM(n,n), RMR(n,n), F(n,n), RF(n,n), RFR(n,n))

     RM  = MATMUL(TRANSPOSE(R), TRANSPOSE(M))
     RMR = MATMUL(RM, R)
     F   = RMR * D !Hadamard product
     RF  = MATMUL(R, F)
     RFR = REAL(MATMUL(RF, TRANSPOSE(R)))

     ! gradient dE/dX has to be (anti)symmetric
     G = 0.5_dp * (TRANSPOSE(RFR) + REAL(sym,dp) * RFR)

     DEALLOCATE(RM, RMR, F, RF, RFR)
  END FUNCTION fold_derivatives


! *****************************************************************************
!> \brief Fills a DBCSR matrix with random noise.
!> \param matrix ...
!> \param seed ...
! *****************************************************************************
  SUBROUTINE fill_randomly(matrix, seed)
    TYPE(cp_dbcsr_type)                      :: matrix
    INTEGER                                  :: seed

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_randomly', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, j
    REAL(dp), DIMENSION(:, :), POINTER       :: block
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(rng_stream_type), POINTER           :: rng_stream

    NULLIFY(rng_stream)
    CALL timeset(routineN,handle)

    CALL cp_dbcsr_iterator_start(iter, matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block)
       ! a "unique" reproducible seed for each block
       initial_seed(1,:) = REAL(seed*arow,KIND=dp)
       initial_seed(2,:) = REAL(seed*acol,KIND=dp)
       initial_seed(3,:) = REAL(seed*SIZE(block),KIND=dp)
       CALL create_rng_stream(rng_stream, seed=initial_seed, name="pao")
       DO i=1, SIZE(block, 1)
          DO j=1, SIZE(block, 2)
             block(i,j) = next_random_number(rng_stream)
          ENDDO
       ENDDO
       CALL delete_rng_stream(rng_stream)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE fill_randomly

END MODULE pao_param_fock

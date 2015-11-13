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
       cp_dbcsr_reserve_diag_blocks
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE pao_types,                       ONLY: pao_env_type
  USE pao_utils,                       ONLY: fill_randomly
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_initial_guess,                ONLY: calculate_atomic_fock_matrix
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

    INTEGER                                  :: acol, arow, handle
    REAL(dp), DIMENSION(:, :), POINTER       :: block_H0
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
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
    CALL fill_randomly(pao%matrix_U0)

    ! make a copy of the atomic guess
    CALL cp_dbcsr_init(pao%matrix_H0)
    CALL cp_dbcsr_create(pao%matrix_H0, template=matrix_s(1)%matrix)
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_H0)

    ! Can't use matrix_ks from ls_scf_qs_atomic_guess(), it's not rotationally invariant
    ! getting H0 directly from the atomic code
    CALL calculate_atomic_fock_matrix(pao%matrix_H0,&
                                      particle_set,&
                                      atomic_kind_set,&
                                      qs_kind_set,&
                                      output_unit=pao%iw)

    ! symmetrize matrix_H0 exactly
    CALL cp_dbcsr_iterator_start(iter, pao%matrix_H0)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_H0)
       block_H0 = 0.5_dp*( block_H0 + TRANSPOSE(block_H0) )
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

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
    CALL cp_dbcsr_release(pao%matrix_H0)

  END SUBROUTINE pao_param_finalize_fock


! *****************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param iatom ...
!> \param X ...
!> \param U ...
!> \param M1 ...
!> \param G ...
! *****************************************************************************
  SUBROUTINE pao_calc_U_fock(pao, iatom, X, U, M1, G)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER                                  :: iatom
    REAL(dp), DIMENSION(:, :), POINTER       :: X, U
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: M1, G

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_fock', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, m, n
    LOGICAL                                  :: found
    REAL(dp)                                 :: denom, nom, v
    REAL(dp), DIMENSION(:), POINTER          :: H_evals, W_evals
    REAL(dp), DIMENSION(:, :), POINTER       :: C, D4, D8, H, H0, H_evecs, &
                                                M2, M3, M4, M5, M6, M7, P, &
                                                U0, W, W_evecs, Y, Z

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_H0, row=iatom, col=iatom, block=H0, found=found)
    CPASSERT(ASSOCIATED(H0))
    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U0, row=iatom, col=iatom, block=U0, found=found)
    CPASSERT(ASSOCIATED(U0))
    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=Y, found=found)
    CPASSERT(ASSOCIATED(Y))

    IF(MAXVAL(ABS(X-TRANSPOSE(X)))>1e-14_dp) CPABORT("Expect symmetric matrix")
    n = SIZE(X, 1)
    m = SIZE(Y, 2)
    IF(SIZE(Y, 1) /= n) CPABORT("strange size of block Y")

    ! calculate P-block from H-block
    ALLOCATE(H_evals(n), H_evecs(n,n), P(n,n), H(n,n))
    H = H0 + X

    IF(MAXVAL(ABS(H-TRANSPOSE(H)))>1e-14_dp) CPABORT("Expect symmetric matrix")

    H_evecs = H
    CALL diamat_all(H_evecs, H_evals)

    P(:,:) = 0.0_dp
    DO k=1, m ! First m states have occupation 1, remaining have 0 .
       DO i=1, n
          DO j=1, n
             P(i,j) = P(i,j) + H_evecs(i,k)*H_evecs(j,k)
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
       v = 1.0_dp / SQRT(W_evals(k))
       DO i=1, m
          DO j=1, m
             Z(i,j) = Z(i,j) + v * W_evecs(i,k)*W_evecs(j,k)
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

       ! dP / dX
       ALLOCATE(D8(n,n))
       DO i=1, n
          DO j=1, n
             IF(i==j) THEN
                  D8(i,i) = 0.0 ! derivative of step-function
             ELSE
                 nom = 0.0_dp
                 IF(i<=m) nom = nom + 1.0_dp
                 IF(j<=m) nom = nom - 1.0_dp
                 denom = H_evals(i) - H_evals(j)
                 IF(nom==0.0_dp) THEN ! avoid division by zero if denom==0 too.
                    D8(i,j) = 0.0_dp
                 ELSE
                    D8(i,j) = nom / denom
                 ENDIF
             ENDIF
          ENDDO
       ENDDO

       ! the final gradient
       G = fold_derivatives(M7, D8, H_evecs, +1)

       DEALLOCATE(M2, M3, M4, M5, M6, M7)
       DEALLOCATE(D4, D8)
    ENDIF

    DEALLOCATE(H, H_evals, H_evecs, W_evals, W_evecs)
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

END MODULE pao_param_fock

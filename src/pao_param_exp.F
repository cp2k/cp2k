!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Original matrix exponential parametrization
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_param_exp
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_create, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_reserve_diag_blocks, cp_dbcsr_set, cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE pao_types,                       ONLY: pao_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_exp'

  PUBLIC :: pao_param_init_exp, pao_param_finalize_exp, pao_calc_U_exp
  PUBLIC :: pao_param_xblock_dims_exp, pao_param_initguess_exp

CONTAINS

! *****************************************************************************
!> \brief Initialize matrix exponential parametrization
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_param_init_exp(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_exp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, natoms
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:), POINTER          :: evals
    REAL(dp), DIMENSION(:, :), POINTER       :: block_H, block_N, block_U0
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    natom=natoms,&
                    dft_control=dft_control)

    IF(dft_control%nspins /= 1) CPABORT("open shell not yet implemented")

    ! allocate matrix_U0
    CALL cp_dbcsr_init(pao%matrix_U0)
    CALL cp_dbcsr_create(pao%matrix_U0, template=matrix_s(1)%matrix, matrix_type="N")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_U0)

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_U0)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_U0)
       CALL cp_dbcsr_get_block_p(matrix=matrix_ks(1)%matrix, row=arow, col=acol, block=block_H, found=found)
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_N, row=arow, col=acol, block=block_N, found=found)
       CPASSERT(ASSOCIATED(block_H) .AND. ASSOCIATED(block_N))
       block_U0 = MATMUL(MATMUL(block_N,block_H),block_N)  ! transform H into orthonormal basis
       ALLOCATE(evals(SIZE(block_H,1)))
       CALL diamat_all(block_U0, evals)
       DEALLOCATE(evals)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE pao_param_init_exp


! *****************************************************************************
!> \brief Finalize exponential parametrization
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_param_finalize_exp(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize_exp', &
      routineP = moduleN//':'//routineN

    CALL cp_dbcsr_release(pao%matrix_U0)

  END SUBROUTINE pao_param_finalize_exp


! *****************************************************************************
!> \brief Returns the dimension of Xblock for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param rows ...
!> \param cols ...
! *****************************************************************************
  SUBROUTINE pao_param_xblock_dims_exp(qs_env, ikind, rows, cols)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: ikind
    INTEGER, INTENT(OUT)                     :: rows, cols

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_xblock_dims_exp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: pao_basis_size, pri_basis_size
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
    CALL get_qs_kind(qs_kind_set(ikind),&
                     basis_set=basis_set,&
                     pao_basis_size=pao_basis_size)
    pri_basis_size = basis_set%nsgf

    ! we only consider rotations between occupied and virtuals
    rows = pao_basis_size
    cols = pri_basis_size - pao_basis_size

  END SUBROUTINE pao_param_xblock_dims_exp


! *****************************************************************************
!> \brief Fills matrix_X with an initial guess
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_param_initguess_exp(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initguess_exp', &
      routineP = moduleN//':'//routineN

    ! matrix_U0 contains actual intial guess
    CALL cp_dbcsr_set(pao%matrix_X, 0.0_dp)

  END SUBROUTINE pao_param_initguess_exp


! *****************************************************************************
!> \brief Calculate new matrix U and optionally its gradient G
!> \param pao ...
!> \param matrix_M ...
!> \param matrix_G ...
! *****************************************************************************
  SUBROUTINE pao_calc_U_exp(pao, matrix_M, matrix_G)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(cp_dbcsr_type), OPTIONAL            :: matrix_M, matrix_G

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_exp', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp)                              :: denom
    COMPLEX(dp), DIMENSION(:), POINTER       :: evals
    COMPLEX(dp), DIMENSION(:, :), POINTER    :: block_D, evecs
    INTEGER                                  :: acol, arow, handle, i, iatom, &
                                                j, k, M, N
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER :: block_G, block_G_full, block_M, &
      block_tmp, block_U, block_U0, block_X, block_X_full, block_Y
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
       IF(arow /= acol) CPABORT("pao_calc_U_exp: encountered off-diagonal block")
       iatom = arow
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=iatom, col=iatom, block=block_U, found=found)
       CPASSERT(ASSOCIATED(block_U))
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U0, row=iatom, col=iatom, block=block_U0, found=found)
       CPASSERT(ASSOCIATED(block_U0))
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=block_Y, found=found)
       CPASSERT(ASSOCIATED(block_Y))

       N=SIZE(block_Y, 1)
       M=SIZE(block_Y, 2)

       ! block_X stores only rotations between occupied and virtuals
       ! hence, we first have to build the full anti-symmetric exponent block
       ALLOCATE(block_X_full(N,N))
       block_X_full(:,:) = 0.0_dp
       block_X_full(1:M, M+1:N) = block_X
       block_X_full(M+1:N, 1:M) = -TRANSPOSE(block_X)

       ! diagonalize block_X_full
       ALLOCATE(evals(N), evecs(N,N))
       CALL diag_antisym(block_X_full, evecs, evals)

       ! construct rotation matrix
       block_U(:,:) = 0.0_dp
       DO k=1, N
          DO i=1, N
             DO j=1, N
                block_U(i,j) = block_U(i,j) + REAL(EXP(evals(k)) * evecs(i,k) * CONJG(evecs(j,k)), dp)
             ENDDO
          ENDDO
       ENDDO

       block_U = MATMUL(block_U0, block_U) ! prepend initial guess rotation

       ! TURNING POINT (if calc grad) ------------------------------------------
       IF(PRESENT(matrix_G))THEN
          CPASSERT(PRESENT(matrix_M))
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
          CPASSERT(ASSOCIATED(block_G))
          CALL cp_dbcsr_get_block_p(matrix=matrix_M, row=iatom, col=iatom, block=block_M, found=found)
          ! don't check ASSOCIATED(block_M), it might have been filtered out.

          ALLOCATE(block_D(N,N), block_tmp(N,N), block_G_full(N,N))
          DO i=1, N
             DO j=1, N
                denom = evals(i) - evals(j)
                IF(i==j) THEN
                   block_D(i,i) = EXP(evals(i)) ! diagonal elements
                ELSE IF(ABS(denom) > 1e-10_dp) THEN
                   block_D(i,j) = (EXP(evals(i)) - EXP(evals(j))) / denom
                ELSE
                   block_D(i,j) = 1.0_dp ! limit according to L'Hospital's rule
                ENDIF
             ENDDO
          ENDDO

          IF (ASSOCIATED(block_M)) THEN
             block_tmp = MATMUL(TRANSPOSE(block_U0), block_M)
          ELSE
             block_tmp = 0.0_dp
          ENDIF
          block_G_full = fold_derivatives(block_tmp, block_D, evecs)

          ! return only gradient for rotations between occupied and virtuals
          IF(M < N)& ! catch special case N==M
             block_G = block_G_full(1:M, M+1:N)

          DEALLOCATE(block_D, block_tmp, block_G_full)
       ENDIF

       DEALLOCATE(block_X_full, evals, evecs)

    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_U_exp


! *****************************************************************************
!> \brief Helper routine, for calculating derivatives
!> \param M ...
!> \param D ...
!> \param R ...
!> \retval G ...
! *****************************************************************************
  FUNCTION fold_derivatives(M, D, R) RESULT(G)
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: M
    COMPLEX(dp), DIMENSION(:, :), INTENT(IN) :: D, R
    REAL(dp), &
      DIMENSION(SIZE(M, 1), SIZE(M, 1))      :: G

    CHARACTER(len=*), PARAMETER :: routineN = 'fold_derivatives', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:, :), POINTER    :: F, RF, RM, RMR
    INTEGER                                  :: n
    REAL(dp), DIMENSION(:, :), POINTER       :: RFR

     n = SIZE(M, 1)

     ALLOCATE(RM(n,n), RMR(n,n), F(n,n), RF(n,n), RFR(n,n))

     RM  = MATMUL(TRANSPOSE(CONJG(R)), TRANSPOSE(M))
     RMR = MATMUL(RM, R)
     F   = RMR * D !Hadamard product
     RF  = MATMUL(R, F)
     RFR = REAL(MATMUL(RF, TRANSPOSE(CONJG(R))))

     ! gradient dE/dX has to be anti-symmetric
     G = TRANSPOSE(RFR) - RFR

     DEALLOCATE(RM, RMR, F, RF, RFR)
  END FUNCTION fold_derivatives


! *****************************************************************************
!> \brief Helper routine for diagonalizing anti symmetric matrices
!> \param matrix ...
!> \param evecs ...
!> \param evals ...
! *****************************************************************************
   SUBROUTINE diag_antisym(matrix, evecs, evals)
    REAL(dp), DIMENSION(:, :)                :: matrix
    COMPLEX(dp), DIMENSION(:, :)             :: evecs
    COMPLEX(dp), DIMENSION(:)                :: evals

    CHARACTER(len=*), PARAMETER :: routineN = 'diag_antisym', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:, :), POINTER    :: matrix_c
    INTEGER                                  :: N
    REAL(dp), DIMENSION(:), POINTER          :: evals_r

    IF(MAXVAL(ABS(matrix+TRANSPOSE(matrix)))>1e-14_dp) CPABORT("Expected anti-symmetric matrix")
    N = SIZE(matrix, 1)
    ALLOCATE(matrix_c(N,N), evals_r(N))

    matrix_c = CMPLX(0.0_dp, -matrix, kind=dp)
    CALL zheevd_wrapper(matrix_c, evecs, evals_r)
    evals = CMPLX(0.0_dp, evals_r, kind=dp)

    DEALLOCATE(matrix_c, evals_r)
   END SUBROUTINE diag_antisym


! *****************************************************************************
!> \brief Helper routine for calling LAPACK zheevd
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! *****************************************************************************
   SUBROUTINE zheevd_wrapper(matrix,eigenvectors,eigenvalues)
    COMPLEX(dp), DIMENSION(:, :)             :: matrix, eigenvectors
    REAL(dp), DIMENSION(:)                   :: eigenvalues

    CHARACTER(len=*), PARAMETER :: routineN = 'zheevd_wrapper', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: work
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: A
    INTEGER                                  :: handle, info, liwork, lrwork, &
                                                lwork, n
    INTEGER, DIMENSION(:), POINTER           :: iwork
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rwork

    CALL timeset(routineN,handle)

    IF(SIZE(matrix, 1) /= SIZE(matrix, 2)) CPABORT("expected square matrix")
    IF(MAXVAL(ABS(matrix-CONJG(TRANSPOSE(matrix))))>1e-14_dp) CPABORT("Expect hermitian matrix")

    n = SIZE(matrix, 1)
    ALLOCATE(iwork(1), rwork(1), work(1), A(n,n))

    A(:,:) = matrix ! ZHEEVD will overwrite A
    ! work space query
    lwork  = -1
    lrwork = -1
    liwork = -1

    CALL ZHEEVD('V','U',n,A(1,1),n,eigenvalues(1),&
                work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)
    lwork  = INT(REAL(work(1)))
    lrwork = INT(REAL(rwork(1)))
    liwork = iwork(1)

    DEALLOCATE(iwork, rwork, work)
    ALLOCATE(iwork(liwork))
    iwork(:) = 0
    ALLOCATE(rwork(lrwork))
    rwork(:) = 0.0_dp
    ALLOCATE(work(lwork))
    work(:) = CMPLX(0.0_dp,0.0_dp,KIND=dp)

    CALL ZHEEVD('V','U',n,A(1,1),n,eigenvalues(1), &
                work(1),lwork,rwork(1),lrwork,iwork(1),liwork,info)

    eigenvectors = A

    IF (info /= 0) CPABORT("diagonalization failed")

    DEALLOCATE(iwork, rwork, work, A)

    CALL timestop(handle)

  END SUBROUTINE zheevd_wrapper

END MODULE pao_param_exp

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Common framework for using eigenvectors of a Fock matrix as PAO basis.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_param_fock
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_create,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_reserve_diag_blocks
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE pao_types,                       ONLY: pao_env_type
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

    INTEGER                                  :: handle
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set)

    ! allocate matrix_H0
    CALL cp_dbcsr_init(pao%matrix_H0)
    CALL cp_dbcsr_create(pao%matrix_H0, template=matrix_s(1)%matrix)
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_H0)

    ! calculate inital atomic fock matrix H0
    ! Can't use matrix_ks from ls_scf_qs_atomic_guess(), it's not rotationally invariant
    ! getting H0 directly from the atomic code
    CALL calculate_atomic_fock_matrix(pao%matrix_H0,&
                                      particle_set,&
                                      atomic_kind_set,&
                                      qs_kind_set,&
                                      output_unit=pao%iw)

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

    CALL cp_dbcsr_release(pao%matrix_H0)

  END SUBROUTINE pao_param_finalize_fock


! *****************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param iatom ...
!> \param penalty ...
!> \param min_gap ...
!> \param V ...
!> \param U ...
!> \param M1 ...
!> \param G ...
! *****************************************************************************
  SUBROUTINE pao_calc_U_fock(pao, iatom, penalty, min_gap, V, U, M1, G)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER                                  :: iatom
    REAL(dp), INTENT(INOUT)                  :: penalty, min_gap
    REAL(dp), DIMENSION(:, :), POINTER       :: V, U
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: M1, G

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_fock', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, m, n
    LOGICAL                                  :: found
    REAL(dp)                                 :: alpha, beta, denom, diff, gap
    REAL(dp), DIMENSION(:), POINTER          :: H_evals
    REAL(dp), DIMENSION(:, :), POINTER       :: block_N, D1, D2, H, H0, &
                                                H_evecs, M2, M3, M4, M5, M6, &
                                                M7, Y

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_H0, row=iatom, col=iatom, block=H0, found=found)
    CPASSERT(ASSOCIATED(H0))
    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_N, row=iatom, col=iatom, block=block_N, found=found)
    CPASSERT(ASSOCIATED(block_N))
    CALL cp_dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=Y, found=found)
    CPASSERT(ASSOCIATED(Y))

    IF(MAXVAL(ABS(V-TRANSPOSE(V)))>1e-14_dp) CPABORT("Expect symmetric matrix")
    n = SIZE(Y, 1) ! size of primary basis
    m = SIZE(Y, 2) ! size of pao basis

    ! calculate H in the orthonormal basis
    ALLOCATE(H(n,n))
    H = MATMUL(MATMUL(block_N, H0 + V), block_N)

    ! diagonalize H
    ALLOCATE(H_evals(n), H_evecs(n,n))
    H_evecs = H
    CALL diamat_all(H_evecs, H_evals)

    ! normalize sign of eigenvectors
    ! The total enegies are unaffacted by this,
    ! however pao_check_grad_param() gets confused by spontaneous sign flips.
    DO i=1, n
       H_evecs(:,i) =  H_evecs(:,i) * SIGN(1.0_dp, H_evecs(i,i))
    ENDDO

    ! the eigenvectors of H become the rotation matrix U
    U = H_evecs

    ! calculate homo-lumo gap (it's useful for detecting numerical issues)
    IF(m < n) THEN ! catch special case n==m
      gap = H_evals(m+1) - H_evals(m)
      min_gap = MIN(min_gap, gap)
    ENDIF

    ! penalty terms: occupied and virtual eigenvalues repel each other
    alpha = pao%penalty_strength
    beta  = pao%penalty_dist
    DO i=1, m
    DO j=m+1, n
       diff = H_evals(i) - H_evals(j)
       penalty = penalty + alpha*EXP(-(diff/beta)**2)
    ENDDO
    ENDDO

    IF(PRESENT(G))THEN ! TURNING POINT (if calc grad) -------------------------

       CPASSERT(PRESENT(M1))

       ! calculate derivatives between eigenvectors of H
       ALLOCATE(D1(n,n), M2(n,n), M3(n,n), M4(n,n))
       DO i=1, n
       DO j=1, n
          ! ignore changes among occupied or virtual eigenvectors
          ! They will get filtered out by M2*D1 anyways, however this early
          ! intervention might stabilize numerics in the case of level-crossings.
          IF(i<=m .EQV. j<=m) THEN
             D1(i,j) = 0.0_dp
          ELSE
             denom = H_evals(i) - H_evals(j)
             IF(ABS(denom) > 1e-9_dp) THEN ! avoid division by zero
                D1(i,j) = 1.0_dp / denom
             ELSE
                D1(i,j) = SIGN(1e+9_dp, denom)
             ENDIF
          ENDIF
       ENDDO
       ENDDO
       IF (ASSOCIATED(M1)) THEN
          M2 = MATMUL(TRANSPOSE(M1), H_evecs)
       ELSE
          M2 = 0.0_dp
       ENDIF
       M3 = M2 * D1 ! Hadamard product
       M4 = MATMUL(MATMUL(H_evecs, M3), TRANSPOSE(H_evecs))

       ! gradient contribution from penalty terms
       ALLOCATE(D2(n,n), M5(n,n))
       D2 = 0.0_dp
       DO i=1, n
       DO j=1, n
          IF (i<=m .EQV. j<=m) CYCLE
          diff = H_evals(i) - H_evals(j)
          D2(i,i) = D2(i,i) - 2.0_dp * alpha * diff / beta**2 * EXP(-(diff/beta)**2)
       ENDDO
       ENDDO
       M5 = MATMUL(MATMUL(H_evecs, D2), TRANSPOSE(H_evecs))

       ! combine M1*dH/dH_evecs with dE_penalty/dH_evecs
       ALLOCATE(M6(n,n))
       M6 = M4 + M5

       ! dH / dV, return to non-orthonormal basis
       ALLOCATE(M7(n,n))
       M7 = MATMUL(MATMUL(block_N, M6), block_N)

       ! symmetrize
       G = 0.5_dp * (M7 + TRANSPOSE(M7)) ! the final gradient

       DEALLOCATE(D1, D2, M2, M3, M4, M5, M6, M7)
    ENDIF

    DEALLOCATE(H, H_evals, H_evecs)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_U_fock

END MODULE pao_param_fock

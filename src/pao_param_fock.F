!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common framework for using eigenvectors of a Fock matrix as PAO basis.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_fock
   USE dbcsr_api,                       ONLY: dbcsr_col_block_sizes,&
                                              dbcsr_get_block_p,&
                                              dbcsr_row_block_sizes
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE pao_types,                       ONLY: pao_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_fock'

   PUBLIC :: pao_calc_U_block_fock

CONTAINS

! **************************************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param iatom ...
!> \param V ...
!> \param U ...
!> \param penalty ...
!> \param gap ...
!> \param evals ...
!> \param M1 ...
!> \param G ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U_block_fock(pao, iatom, V, U, penalty, gap, evals, M1, G)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :), POINTER                 :: V, U
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty
      REAL(dp), INTENT(OUT)                              :: gap
      REAL(dp), DIMENSION(:), INTENT(OUT), OPTIONAL      :: evals
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: M1, G

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_block_fock', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, m, n
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      LOGICAL                                            :: found
      REAL(dp)                                           :: alpha, beta, denom, diff
      REAL(dp), DIMENSION(:), POINTER                    :: H_evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_N, D1, D2, H, H0, H_evecs, M2, M3, &
                                                            M4, M5

      CALL timeset(routineN, handle)

      CALL dbcsr_get_block_p(matrix=pao%matrix_H0, row=iatom, col=iatom, block=H0, found=found)
      CPASSERT(ASSOCIATED(H0))
      CALL dbcsr_get_block_p(matrix=pao%matrix_N_diag, row=iatom, col=iatom, block=block_N, found=found)
      CPASSERT(ASSOCIATED(block_N))
      IF (MAXVAL(ABS(V-TRANSPOSE(V))) > 1e-14_dp) CPABORT("Expect symmetric matrix")

      ! figure out basis sizes
      blk_sizes_pri => dbcsr_row_block_sizes(pao%matrix_Y)
      blk_sizes_pao => dbcsr_col_block_sizes(pao%matrix_Y)
      n = blk_sizes_pri(iatom) ! size of primary basis
      m = blk_sizes_pao(iatom) ! size of pao basis

      ! calculate H in the orthonormal basis
      ALLOCATE (H(n, n))
      H = MATMUL(MATMUL(block_N, H0+V), block_N)

      ! diagonalize H
      ALLOCATE (H_evals(n), H_evecs(n, n))
      H_evecs = H
      CALL diamat_all(H_evecs, H_evals)

      ! the eigenvectors of H become the rotation matrix U
      U = H_evecs

      ! copy eigenvectors around the gap from H_evals into evals array
      IF (PRESENT(evals)) THEN
         CPASSERT(MOD(SIZE(evals), 2) == 0) ! gap will be exactely in the middle
         i = SIZE(evals)/2
         j = MIN(m, i)
         evals(1+i-j:i) = H_evals(1+m-j:m) ! eigenvalues below gap
         j = MIN(n-m, i)
         evals(i:i+j) = H_evals(m:m+j) ! eigenvalues above gap
      ENDIF

      ! calculate homo-lumo gap (it's useful for detecting numerical issues)
      gap = HUGE(dp)
      IF (m < n) & ! catch special case n==m
         gap = H_evals(m+1)-H_evals(m)

      IF (PRESENT(penalty)) THEN
         ! penalty terms: occupied and virtual eigenvalues repel each other
         alpha = pao%penalty_strength
         beta = pao%penalty_dist
         DO i = 1, m
         DO j = m+1, n
            diff = H_evals(i)-H_evals(j)
            penalty = penalty+alpha*EXP(-(diff/beta)**2)
         ENDDO
         ENDDO

         ! regularization energy
         penalty = penalty+pao%regularization*SUM(V**2)
      ENDIF

      IF (PRESENT(G)) THEN ! TURNING POINT (if calc grad) -------------------------

         CPASSERT(PRESENT(M1))

         ! calculate derivatives between eigenvectors of H
         ALLOCATE (D1(n, n), M2(n, n), M3(n, n), M4(n, n))
         DO i = 1, n
         DO j = 1, n
            ! ignore changes among occupied or virtual eigenvectors
            ! They will get filtered out by M2*D1 anyways, however this early
            ! intervention might stabilize numerics in the case of level-crossings.
            IF (i <= m .EQV. j <= m) THEN
               D1(i, j) = 0.0_dp
            ELSE
               denom = H_evals(i)-H_evals(j)
               IF (ABS(denom) > 1e-9_dp) THEN ! avoid division by zero
                  D1(i, j) = 1.0_dp/denom
               ELSE
                  D1(i, j) = SIGN(1e+9_dp, denom)
               ENDIF
            ENDIF
         ENDDO
         ENDDO
         IF (ASSOCIATED(M1)) THEN
            M2 = MATMUL(TRANSPOSE(M1), H_evecs)
         ELSE
            M2 = 0.0_dp
         ENDIF
         M3 = M2*D1 ! Hadamard product
         M4 = MATMUL(MATMUL(H_evecs, M3), TRANSPOSE(H_evecs))

         ! gradient contribution from penalty terms
         IF (PRESENT(penalty)) THEN
            ALLOCATE (D2(n, n))
            D2 = 0.0_dp
            DO i = 1, n
            DO j = 1, n
               IF (i <= m .EQV. j <= m) CYCLE
               diff = H_evals(i)-H_evals(j)
               D2(i, i) = D2(i, i)-2.0_dp*alpha*diff/beta**2*EXP(-(diff/beta)**2)
            ENDDO
            ENDDO
            M4 = M4+MATMUL(MATMUL(H_evecs, D2), TRANSPOSE(H_evecs))
            DEALLOCATE (D2)
         ENDIF

         ! dH / dV, return to non-orthonormal basis
         ALLOCATE (M5(n, n))
         M5 = MATMUL(MATMUL(block_N, M4), block_N)

         ! add regularization gradient
         IF (PRESENT(penalty)) &
            M5 = M5+2.0_dp*pao%regularization*V

         ! symmetrize
         G = 0.5_dp*(M5+TRANSPOSE(M5)) ! the final gradient

         DEALLOCATE (D1, M2, M3, M4, M5)
      ENDIF

      DEALLOCATE (H, H_evals, H_evecs)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_U_block_fock

END MODULE pao_param_fock

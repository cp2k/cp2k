!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Equivariant parametrization
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_equi
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_complete_redistribute, dbcsr_create, dbcsr_distribution_type, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_comm_type
   USE pao_param_methods,               ONLY: pao_calc_grad_lnv_wrt_AB
   USE pao_potentials,                  ONLY: pao_guess_initial_potential
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_equi'

   PUBLIC :: pao_param_init_equi, pao_param_finalize_equi, pao_calc_AB_equi
   PUBLIC :: pao_param_count_equi, pao_param_initguess_equi

CONTAINS

! **************************************************************************************************
!> \brief Initialize equivariant parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_equi(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      IF (pao%precondition) &
         CPABORT("PAO preconditioning not supported for selected parametrization.")

   END SUBROUTINE pao_param_init_equi

! **************************************************************************************************
!> \brief Finalize equivariant parametrization
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_equi()

      ! Nothing to do.

   END SUBROUTINE pao_param_finalize_equi

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_equi(qs_env, ikind, nparams)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      INTEGER                                            :: pao_basis_size, pri_basis_size
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), &
                       basis_set=basis_set, &
                       pao_basis_size=pao_basis_size)
      pri_basis_size = basis_set%nsgf

      nparams = pao_basis_size*pri_basis_size

   END SUBROUTINE pao_param_count_equi

! **************************************************************************************************
!> \brief Fills matrix_X with an initial guess
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_equi(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initguess_equi'

      INTEGER                                            :: acol, arow, handle, i, iatom, m, n
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:), POINTER                    :: H_evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: A, block_H0, block_N, block_N_inv, &
                                                            block_X, H, H_evecs, V0
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,blk_sizes_pri,blk_sizes_pao) &
!$OMP PRIVATE(iter,arow,acol,iatom,n,m,i,found) &
!$OMP PRIVATE(block_X,block_H0,block_N,block_N_inv,A,H,H_evecs,H_evals,V0)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)

         CALL dbcsr_get_block_p(matrix=pao%matrix_H0, row=iatom, col=iatom, block=block_H0, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_diag, row=iatom, col=iatom, block=block_N, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_inv_diag, row=iatom, col=iatom, block=block_N_inv, found=found)
         CPASSERT(ASSOCIATED(block_H0) .AND. ASSOCIATED(block_N) .AND. ASSOCIATED(block_N_inv))

         n = blk_sizes_pri(iatom) ! size of primary basis
         m = blk_sizes_pao(iatom) ! size of pao basis

         ALLOCATE (V0(n, n))
         CALL pao_guess_initial_potential(qs_env, iatom, V0)

         ! construct H
         ALLOCATE (H(n, n))
         H = MATMUL(MATMUL(block_N, block_H0 + V0), block_N) ! transform into orthonormal basis

         ! diagonalize H
         ALLOCATE (H_evecs(n, n), H_evals(n))
         H_evecs = H
         CALL diamat_all(H_evecs, H_evals)

         ! use first m eigenvectors as initial guess
         ALLOCATE (A(n, m))
         A = MATMUL(block_N_inv, H_evecs(:, 1:m))

         ! normalize vectors
         DO i = 1, m
            A(:, i) = A(:, i)/NORM2(A(:, i))
         END DO

         block_X = RESHAPE(A, (/n*m, 1/))
         DEALLOCATE (H, V0, A, H_evecs, H_evals)

      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_param_initguess_equi

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
!> \param penalty ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_equi(pao, qs_env, ls_scf_env, gradient, penalty)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_AB_equi'

      INTEGER                                            :: acol, arow, group_handle, handle, i, &
                                                            iatom, j, k, m, n
      LOGICAL                                            :: found
      REAL(dp)                                           :: denom, w
      REAL(dp), DIMENSION(:), POINTER                    :: ANNA_evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: ANNA, ANNA_evecs, ANNA_inv, block_A, &
                                                            block_B, block_G, block_Ma, block_Mb, &
                                                            block_N, block_X, D, G, M1, M2, M3, &
                                                            M4, M5, NN
      TYPE(dbcsr_distribution_type)                      :: main_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_G_nondiag, matrix_Ma, matrix_Mb, &
                                                            matrix_X_nondiag
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      ls_mstruct => ls_scf_env%ls_mstruct

      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)
      END IF

      ! Redistribute matrix_X from diag_distribution to distribution of matrix_s.
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix=matrix_s(1)%matrix, distribution=main_dist)
      CALL dbcsr_create(matrix_X_nondiag, &
                        name="PAO matrix_X_nondiag", &
                        dist=main_dist, &
                        template=pao%matrix_X)
      CALL dbcsr_reserve_diag_blocks(matrix_X_nondiag)
      CALL dbcsr_complete_redistribute(pao%matrix_X, matrix_X_nondiag)

      ! Compuation of matrix_G uses distr. of matrix_s, afterwards we redistribute to diag_distribution.
      IF (gradient) THEN
         CALL dbcsr_create(matrix_G_nondiag, &
                           name="PAO matrix_G_nondiag", &
                           dist=main_dist, &
                           template=pao%matrix_G)
         CALL dbcsr_reserve_diag_blocks(matrix_G_nondiag)
      END IF

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(pao,ls_mstruct,matrix_X_nondiag,matrix_G_nondiag,matrix_Ma,matrix_Mb,gradient,penalty) &
!$OMP PRIVATE(iter,arow,acol,iatom,found,n,m,w,i,j,k,denom) &
!$OMP PRIVATE(NN,ANNA,ANNA_evals,ANNA_evecs,ANNA_inv,D,G,M1,M2,M3,M4,M5) &
!$OMP PRIVATE(block_X,block_A,block_B,block_N,block_Ma, block_Mb, block_G)
      CALL dbcsr_iterator_start(iter, matrix_X_nondiag)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_A, row=iatom, col=iatom, block=block_A, found=found)
         CPASSERT(ASSOCIATED(block_A))
         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_B, row=iatom, col=iatom, block=block_B, found=found)
         CPASSERT(ASSOCIATED(block_B))
         CALL dbcsr_get_block_p(matrix=pao%matrix_N, row=iatom, col=iatom, block=block_N, found=found)
         CPASSERT(ASSOCIATED(block_N))

         n = SIZE(block_A, 1) ! size of primary basis
         m = SIZE(block_A, 2) ! size of pao basis
         block_A = RESHAPE(block_X, (/n, m/))

         ! restrain pao basis vectors to unit norm
         IF (PRESENT(penalty)) THEN
            DO i = 1, m
               w = 1.0_dp - SUM(block_A(:, i)**2)
               penalty = penalty + pao%penalty_strength*w**2
            END DO
         END IF

         ALLOCATE (NN(n, n), ANNA(m, m))
         NN = MATMUL(block_N, block_N) ! it's actually S^{-1}
         ANNA = MATMUL(MATMUL(TRANSPOSE(block_A), NN), block_A)

         ! diagonalize ANNA
         ALLOCATE (ANNA_evecs(m, m), ANNA_evals(m))
         ANNA_evecs(:, :) = ANNA
         CALL diamat_all(ANNA_evecs, ANNA_evals)
         IF (MINVAL(ABS(ANNA_evals)) < 1e-10_dp) CPABORT("PAO basis singualar.")

         ! build ANNA_inv
         ALLOCATE (ANNA_inv(m, m))
         ANNA_inv(:, :) = 0.0_dp
         DO k = 1, m
            w = 1.0_dp/ANNA_evals(k)
            DO i = 1, m
            DO j = 1, m
               ANNA_inv(i, j) = ANNA_inv(i, j) + w*ANNA_evecs(i, k)*ANNA_evecs(j, k)
            END DO
            END DO
         END DO

         !B = 1/S * A * 1/(A^T 1/S A)
         block_B = MATMUL(MATMUL(NN, block_A), ANNA_inv)

         ! TURNING POINT (if calc grad) ------------------------------------------
         IF (gradient) THEN
            CALL dbcsr_get_block_p(matrix=matrix_G_nondiag, row=iatom, col=iatom, block=block_G, found=found)
            CPASSERT(ASSOCIATED(block_G))
            CALL dbcsr_get_block_p(matrix=matrix_Ma, row=iatom, col=iatom, block=block_Ma, found=found)
            CALL dbcsr_get_block_p(matrix=matrix_Mb, row=iatom, col=iatom, block=block_Mb, found=found)
            ! don't check ASSOCIATED(block_M), it might have been filtered out.

            ALLOCATE (G(n, m))
            G(:, :) = 0.0_dp

            IF (PRESENT(penalty)) THEN
               DO i = 1, m
                  w = 1.0_dp - SUM(block_A(:, i)**2)
                  G(:, i) = -4.0_dp*pao%penalty_strength*w*block_A(:, i)
               END DO
            END IF

            IF (ASSOCIATED(block_Ma)) THEN
               G = G + block_Ma
            END IF

            IF (ASSOCIATED(block_Mb)) THEN
               G = G + MATMUL(MATMUL(NN, block_Mb), ANNA_inv)

               ! calculate derivatives dAA_inv/ dAA
               ALLOCATE (D(m, m), M1(m, m), M2(m, m), M3(m, m), M4(m, m), M5(m, m))

               DO i = 1, m
               DO j = 1, m
                  denom = ANNA_evals(i) - ANNA_evals(j)
                  IF (i == j) THEN
                     D(i, i) = -1.0_dp/ANNA_evals(i)**2 ! diagonal elements
                  ELSE IF (ABS(denom) > 1e-10_dp) THEN
                     D(i, j) = (1.0_dp/ANNA_evals(i) - 1.0_dp/ANNA_evals(j))/denom
                  ELSE
                     D(i, j) = -1.0_dp ! limit according to L'Hospital's rule
                  END IF
               END DO
               END DO

               M1 = MATMUL(MATMUL(TRANSPOSE(block_A), NN), block_Mb)
               M2 = MATMUL(MATMUL(TRANSPOSE(ANNA_evecs), M1), ANNA_evecs)
               M3 = M2*D ! Hadamard product
               M4 = MATMUL(MATMUL(ANNA_evecs, M3), TRANSPOSE(ANNA_evecs))
               M5 = 0.5_dp*(M4 + TRANSPOSE(M4))
               G = G + 2.0_dp*MATMUL(MATMUL(NN, block_A), M5)

               DEALLOCATE (D, M1, M2, M3, M4, M5)
            END IF

            block_G = RESHAPE(G, (/n*m, 1/))
            DEALLOCATE (G)
         END IF

         DEALLOCATE (NN, ANNA, ANNA_evecs, ANNA_evals, ANNA_inv)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      ! sum penalty energies across ranks
      IF (PRESENT(penalty)) THEN
         CALL dbcsr_get_info(pao%matrix_X, group=group_handle)
         CALL group%set_handle(group_handle)
         CALL group%sum(penalty)
      END IF

      CALL dbcsr_release(matrix_X_nondiag)

      IF (gradient) THEN
         CALL dbcsr_complete_redistribute(matrix_G_nondiag, pao%matrix_G)
         CALL dbcsr_release(matrix_G_nondiag)
         CALL dbcsr_release(matrix_Ma)
         CALL dbcsr_release(matrix_Mb)
      END IF

      CALL timestop(handle)

   END SUBROUTINE pao_calc_AB_equi

END MODULE pao_param_equi

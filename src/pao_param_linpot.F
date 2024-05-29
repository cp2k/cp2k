!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common framework for a linear parametrization of the potential.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_linpot
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE pao_input,                       ONLY: pao_fock_param,&
                                              pao_rotinv_param
   USE pao_linpot_full,                 ONLY: linpot_full_calc_terms,&
                                              linpot_full_count_terms
   USE pao_linpot_rotinv,               ONLY: linpot_rotinv_calc_forces,&
                                              linpot_rotinv_calc_terms,&
                                              linpot_rotinv_count_terms
   USE pao_param_fock,                  ONLY: pao_calc_U_block_fock
   USE pao_param_methods,               ONLY: pao_calc_AB_from_U,&
                                              pao_calc_grad_lnv_wrt_U
   USE pao_potentials,                  ONLY: pao_guess_initial_potential
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_param_init_linpot, pao_param_finalize_linpot, pao_calc_AB_linpot
   PUBLIC :: pao_param_count_linpot, pao_param_initguess_linpot

CONTAINS

! **************************************************************************************************
!> \brief Initialize the linear potential parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_linpot(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_linpot'

      INTEGER                                            :: acol, arow, handle, iatom, ikind, N, &
                                                            natoms, nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pri, col_blk_size, row_blk_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_V_terms
      REAL(dp), DIMENSION(:, :, :), POINTER              :: V_blocks
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      natom=natoms)

      IF (dft_control%nspins /= 1) CPABORT("open shell not yet implemented")

      ! figure out number of potential terms
      ALLOCATE (row_blk_size(natoms), col_blk_size(natoms))
      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL pao_param_count_linpot(pao, qs_env, ikind, nterms)
         col_blk_size(iatom) = nterms
      END DO

      ! allocate matrix_V_terms
      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri)
      row_blk_size = blk_sizes_pri**2
      CALL dbcsr_create(pao%matrix_V_terms, &
                        name="PAO matrix_V_terms", &
                        dist=pao%diag_distribution, &
                        matrix_type="N", &
                        row_blk_size=row_blk_size, &
                        col_blk_size=col_blk_size)
      CALL dbcsr_reserve_diag_blocks(pao%matrix_V_terms)
      DEALLOCATE (row_blk_size, col_blk_size)

      ! calculate, normalize, and store potential terms as rows of block_V_terms
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,blk_sizes_pri) &
!$OMP PRIVATE(iter,arow,acol,iatom,N,nterms,block_V_terms,V_blocks)
      CALL dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_V_terms)
         iatom = arow; CPASSERT(arow == acol)
         nterms = SIZE(block_V_terms, 2)
         IF (nterms == 0) CYCLE ! protect against corner-case of zero pao parameters
         N = blk_sizes_pri(iatom)
         CPASSERT(N*N == SIZE(block_V_terms, 1))
         ALLOCATE (V_blocks(N, N, nterms))
         CALL linpot_calc_terms(pao, qs_env, iatom, V_blocks)
         block_V_terms = RESHAPE(V_blocks, (/N*N, nterms/)) ! convert matrices into vectors
         DEALLOCATE (V_blocks)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL pao_param_linpot_regularizer(pao)

      IF (pao%precondition) &
         CALL pao_param_linpot_preconditioner(pao)

      CALL para_env%sync() ! ensure that timestop is not called too early

      CALL timestop(handle)
   END SUBROUTINE pao_param_init_linpot

! **************************************************************************************************
!> \brief Builds the regularization metric matrix_R
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_linpot_regularizer(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_linpot_regularizer'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, &
                                                            nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_nterms
      LOGICAL                                            :: found
      REAL(dp)                                           :: v, w
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: S, S_evecs
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_R, V_terms
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Building linpot regularizer"

      CALL dbcsr_get_info(pao%matrix_V_terms, col_blk_size=blk_sizes_nterms)

      ! build regularization metric
      CALL dbcsr_create(pao%matrix_R, &
                        template=pao%matrix_V_terms, &
                        matrix_type="N", &
                        row_blk_size=blk_sizes_nterms, &
                        col_blk_size=blk_sizes_nterms, &
                        name="PAO matrix_R")
      CALL dbcsr_reserve_diag_blocks(pao%matrix_R)

      ! fill matrix_R
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_R,V_terms,found,nterms,S,S_evecs,S_evals,k,i,j,v,w)
      CALL dbcsr_iterator_start(iter, pao%matrix_R)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_R)
         iatom = arow; CPASSERT(arow == acol)
         CALL dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=V_terms, found=found)
         CPASSERT(ASSOCIATED(V_terms))
         nterms = SIZE(V_terms, 2)
         IF (nterms == 0) CYCLE ! protect against corner-case of zero pao parameters

         ! build overlap matrix
         ALLOCATE (S(nterms, nterms))
         S(:, :) = MATMUL(TRANSPOSE(V_terms), V_terms)

         ! diagonalize S
         ALLOCATE (S_evals(nterms), S_evecs(nterms, nterms))
         S_evecs(:, :) = S
         CALL diamat_all(S_evecs, S_evals)

         block_R = 0.0_dp
         DO k = 1, nterms
            v = pao%linpot_regu_delta/S_evals(k)
            w = pao%linpot_regu_strength*MIN(1.0_dp, ABS(v))
            DO i = 1, nterms
            DO j = 1, nterms
               block_R(i, j) = block_R(i, j) + w*S_evecs(i, k)*S_evecs(j, k)
            END DO
            END DO
         END DO

         ! clean up
         DEALLOCATE (S, S_evals, S_evecs)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE pao_param_linpot_regularizer

! **************************************************************************************************
!> \brief Builds the preconditioner matrix_precon and matrix_precon_inv
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_linpot_preconditioner(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_linpot_preconditioner'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, &
                                                            nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_nterms
      LOGICAL                                            :: found
      REAL(dp)                                           :: eval_capped
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: S, S_evecs
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_precon, block_precon_inv, &
                                                            block_V_terms
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Building linpot preconditioner"

      CALL dbcsr_get_info(pao%matrix_V_terms, col_blk_size=blk_sizes_nterms)

      CALL dbcsr_create(pao%matrix_precon, &
                        template=pao%matrix_V_terms, &
                        matrix_type="N", &
                        row_blk_size=blk_sizes_nterms, &
                        col_blk_size=blk_sizes_nterms, &
                        name="PAO matrix_precon")
      CALL dbcsr_reserve_diag_blocks(pao%matrix_precon)

      CALL dbcsr_create(pao%matrix_precon_inv, template=pao%matrix_precon, name="PAO matrix_precon_inv")
      CALL dbcsr_reserve_diag_blocks(pao%matrix_precon_inv)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_V_terms,block_precon,block_precon_inv,found,nterms,S,S_evals,S_evecs,i,j,k,eval_capped)
      CALL dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_V_terms)
         iatom = arow; CPASSERT(arow == acol)
         nterms = SIZE(block_V_terms, 2)
         IF (nterms == 0) CYCLE ! protect against corner-case of zero pao parameters

         CALL dbcsr_get_block_p(matrix=pao%matrix_precon, row=iatom, col=iatom, block=block_precon, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_precon_inv, row=iatom, col=iatom, block=block_precon_inv, found=found)
         CPASSERT(ASSOCIATED(block_precon))
         CPASSERT(ASSOCIATED(block_precon_inv))

         ALLOCATE (S(nterms, nterms))
         S(:, :) = MATMUL(TRANSPOSE(block_V_terms), block_V_terms)

         ! diagonalize S
         ALLOCATE (S_evals(nterms), S_evecs(nterms, nterms))
         S_evecs(:, :) = S
         CALL diamat_all(S_evecs, S_evals)

         ! construct 1/Sqrt(S) and Sqrt(S)
         block_precon = 0.0_dp
         block_precon_inv = 0.0_dp
         DO k = 1, nterms
            eval_capped = MAX(pao%linpot_precon_delta, S_evals(k)) ! too small eigenvalues are hurtful
            DO i = 1, nterms
            DO j = 1, nterms
               block_precon(i, j) = block_precon(i, j) + S_evecs(i, k)*S_evecs(j, k)/SQRT(eval_capped)
               block_precon_inv(i, j) = block_precon_inv(i, j) + S_evecs(i, k)*S_evecs(j, k)*SQRT(eval_capped)
            END DO
            END DO
         END DO

         DEALLOCATE (S, S_evecs, S_evals)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE pao_param_linpot_preconditioner

! **************************************************************************************************
!> \brief Finalize the linear potential parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_linpot(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL dbcsr_release(pao%matrix_V_terms)
      CALL dbcsr_release(pao%matrix_R)

      IF (pao%precondition) THEN
         CALL dbcsr_release(pao%matrix_precon)
         CALL dbcsr_release(pao%matrix_precon_inv)
      END IF

   END SUBROUTINE pao_param_finalize_linpot

! **************************************************************************************************
!> \brief Returns the number of potential terms for given atomic kind
!> \param pao ...
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_linpot(pao, qs_env, ikind, nparams)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      INTEGER                                            :: pao_basis_size
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)

      CALL get_qs_kind(qs_kind_set(ikind), &
                       basis_set=basis_set, &
                       pao_basis_size=pao_basis_size)

      IF (pao_basis_size == basis_set%nsgf) THEN
         nparams = 0 ! pao disabled for iatom

      ELSE
         SELECT CASE (pao%parameterization)
         CASE (pao_fock_param)
            CALL linpot_full_count_terms(qs_env, ikind, nterms=nparams)
         CASE (pao_rotinv_param)
            CALL linpot_rotinv_count_terms(qs_env, ikind, nterms=nparams)
         CASE DEFAULT
            CPABORT("unknown parameterization")
         END SELECT
      END IF

   END SUBROUTINE pao_param_count_linpot

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
!> \param penalty ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_linpot(pao, qs_env, ls_scf_env, gradient, penalty, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_AB_linpot'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_M, matrix_U

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_create(matrix_U, matrix_type="N", dist=pao%diag_distribution, template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U)

      !TODO: move this condition into pao_calc_U, use matrix_N as template
      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_U(qs_env, ls_scf_env, matrix_M)
         CALL pao_calc_U_linpot(pao, qs_env, matrix_U, matrix_M, pao%matrix_G, penalty, forces)
         CALL dbcsr_release(matrix_M)
      ELSE
         CALL pao_calc_U_linpot(pao, qs_env, matrix_U, penalty=penalty)
      END IF

      CALL pao_calc_AB_from_U(pao, qs_env, ls_scf_env, matrix_U)
      CALL dbcsr_release(matrix_U)
      CALL timestop(handle)
   END SUBROUTINE pao_calc_AB_linpot

! **************************************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param qs_env ...
!> \param matrix_U ...
!> \param matrix_M ...
!> \param matrix_G ...
!> \param penalty ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U_linpot(pao, qs_env, matrix_U, matrix_M, matrix_G, penalty, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type)                                   :: matrix_U
      TYPE(dbcsr_type), OPTIONAL                         :: matrix_M, matrix_G
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U_linpot'

      INTEGER                                            :: acol, arow, group_handle, handle, iatom, &
                                                            kterm, n, natoms, nterms
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: gaps
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: evals
      REAL(dp), DIMENSION(:), POINTER                    :: vec_M2, vec_V
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_M1, block_M2, block_R, &
                                                            block_U, block_V, block_V_terms, &
                                                            block_X
      REAL(dp), DIMENSION(:, :, :), POINTER              :: M_blocks
      REAL(KIND=dp)                                      :: regu_energy
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(matrix_G) .EQV. PRESENT(matrix_M))

      CALL get_qs_env(qs_env, natom=natoms)
      ALLOCATE (gaps(natoms), evals(10, natoms)) ! printing 10 eigenvalues seems reasonable
      evals(:, :) = 0.0_dp
      gaps(:) = HUGE(1.0_dp)
      regu_energy = 0.0_dp
      CALL dbcsr_get_info(matrix_U, group=group_handle)
      CALL group%set_handle(group_handle)

      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL dbcsr_get_block_p(matrix=pao%matrix_R, row=iatom, col=iatom, block=block_R, found=found)
         CALL dbcsr_get_block_p(matrix=matrix_U, row=iatom, col=iatom, block=block_U, found=found)
         CPASSERT(ASSOCIATED(block_R) .AND. ASSOCIATED(block_U))
         n = SIZE(block_U, 1)

         ! calculate potential V
         ALLOCATE (vec_V(n*n))
         vec_V(:) = 0.0_dp
         CALL dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=block_V_terms, found=found)
         CPASSERT(ASSOCIATED(block_V_terms))
         nterms = SIZE(block_V_terms, 2)
         IF (nterms > 0) & ! protect against corner-case of zero pao parameters
            vec_V = MATMUL(block_V_terms, block_X(:, 1))
         block_V(1:n, 1:n) => vec_V(:) ! map vector into matrix

         ! symmetrize
         IF (MAXVAL(ABS(block_V - TRANSPOSE(block_V))/MAX(1.0_dp, MAXVAL(ABS(block_V)))) > 1e-12) &
            CPABORT("block_V not symmetric")
         block_V = 0.5_dp*(block_V + TRANSPOSE(block_V)) ! symmetrize exactly

         ! regularization energy
         ! protect against corner-case of zero pao parameters
         IF (PRESENT(penalty) .AND. nterms > 0) &
            regu_energy = regu_energy + DOT_PRODUCT(block_X(:, 1), MATMUL(block_R, block_X(:, 1)))

         CALL pao_calc_U_block_fock(pao, iatom=iatom, penalty=penalty, V=block_V, U=block_U, &
                                    gap=gaps(iatom), evals=evals(:, iatom))

         IF (PRESENT(matrix_G)) THEN ! TURNING POINT (if calc grad) --------------------------------
            CPASSERT(PRESENT(matrix_M))
            CALL dbcsr_get_block_p(matrix=matrix_M, row=iatom, col=iatom, block=block_M1, found=found)

            ! corner-cases: block_M1 might have been filtered out or there might be zero pao parameters
            IF (ASSOCIATED(block_M1) .AND. SIZE(block_V_terms) > 0) THEN
               ALLOCATE (vec_M2(n*n))
               block_M2(1:n, 1:n) => vec_M2(:) ! map vector into matrix
               !TODO: this 2nd call does double work. However, *sometimes* this branch is not taken.
               CALL pao_calc_U_block_fock(pao, iatom=iatom, penalty=penalty, V=block_V, U=block_U, &
                                          M1=block_M1, G=block_M2, gap=gaps(iatom), evals=evals(:, iatom))
               IF (MAXVAL(ABS(block_M2 - TRANSPOSE(block_M2))) > 1e-14_dp) &
                  CPABORT("matrix not symmetric")

               ! gradient dE/dX
               IF (PRESENT(matrix_G)) THEN
                  CALL dbcsr_get_block_p(matrix=matrix_G, row=iatom, col=iatom, block=block_G, found=found)
                  CPASSERT(ASSOCIATED(block_G))
                  block_G(:, 1) = MATMUL(vec_M2, block_V_terms)
                  IF (PRESENT(penalty)) &
                     block_G = block_G + 2.0_dp*MATMUL(block_R, block_X) ! regularization gradient
               END IF

               ! forced dE/dR
               IF (PRESENT(forces)) THEN
                  ALLOCATE (M_blocks(n, n, nterms))
                  DO kterm = 1, nterms
                     M_blocks(:, :, kterm) = block_M2*block_X(kterm, 1)
                  END DO
                  CALL linpot_calc_forces(pao, qs_env, iatom=iatom, M_blocks=M_blocks, forces=forces)
                  DEALLOCATE (M_blocks)
               END IF

               DEALLOCATE (vec_M2)
            END IF
         END IF
         DEALLOCATE (vec_V)
      END DO
      CALL dbcsr_iterator_stop(iter)

      IF (PRESENT(penalty)) THEN
         ! sum penalty energies across ranks
         CALL group%sum(penalty)
         CALL group%sum(regu_energy)
         penalty = penalty + regu_energy
      END IF

      ! print stuff, but not during second invocation for forces
      IF (.NOT. PRESENT(forces)) THEN
         ! print eigenvalues from fock-layer
         CALL group%sum(evals)
         IF (pao%iw_fockev > 0) THEN
            DO iatom = 1, natoms
               WRITE (pao%iw_fockev, *) "PAO| atom:", iatom, " fock evals around gap:", evals(:, iatom)
            END DO
            CALL m_flush(pao%iw_fockev)
         END IF
         ! print homo-lumo gap encountered by fock-layer
         CALL group%min(gaps)
         IF (pao%iw_gap > 0) THEN
            DO iatom = 1, natoms
               WRITE (pao%iw_gap, *) "PAO| atom:", iatom, " fock gap:", gaps(iatom)
            END DO
         END IF
         ! one-line summaries
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| linpot regularization energy:", regu_energy
         IF (pao%iw > 0) WRITE (pao%iw, "(A,E20.10,A,T71,I10)") " PAO| min_gap:", MINVAL(gaps), " for atom:", MINLOC(gaps)
      END IF

      DEALLOCATE (gaps, evals)
      CALL timestop(handle)

   END SUBROUTINE pao_calc_U_linpot

! **************************************************************************************************
!> \brief Internal routine, calculates terms in potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param V_blocks ...
! **************************************************************************************************
   SUBROUTINE linpot_calc_terms(pao, qs_env, iatom, V_blocks)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :, :), INTENT(OUT)          :: V_blocks

      SELECT CASE (pao%parameterization)
      CASE (pao_fock_param)
         CALL linpot_full_calc_terms(V_blocks)
      CASE (pao_rotinv_param)
         CALL linpot_rotinv_calc_terms(qs_env, iatom, V_blocks)
      CASE DEFAULT
         CPABORT("unknown parameterization")
      END SELECT

   END SUBROUTINE linpot_calc_terms

! **************************************************************************************************
!> \brief Internal routine, calculates force contributions from potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param M_blocks ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE linpot_calc_forces(pao, qs_env, iatom, M_blocks, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: M_blocks
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      SELECT CASE (pao%parameterization)
      CASE (pao_fock_param)
         ! no force contributions
      CASE (pao_rotinv_param)
         CALL linpot_rotinv_calc_forces(qs_env, iatom, M_blocks, forces)
      CASE DEFAULT
         CPABORT("unknown parameterization")
      END SELECT

   END SUBROUTINE linpot_calc_forces

! **************************************************************************************************
!> \brief Calculate initial guess for matrix_X
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_linpot(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initguess_linpot'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, n, &
                                                            nterms
      INTEGER, DIMENSION(:), POINTER                     :: pri_basis_size
      LOGICAL                                            :: found
      REAL(dp)                                           :: w
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: S, S_evecs, S_inv
      REAL(dp), DIMENSION(:), POINTER                    :: V_guess_vec
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X, V_guess, V_terms
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=pri_basis_size)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,pri_basis_size) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_X,N,nterms,V_terms,found,V_guess,V_guess_vec,S,S_evecs,S_evals,S_inv,k,i,j,w)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=V_terms, found=found)
         CPASSERT(ASSOCIATED(V_terms))
         nterms = SIZE(V_terms, 2)
         IF (nterms == 0) CYCLE ! protect against corner-case of zero pao parameters

         ! guess initial potential
         N = pri_basis_size(iatom)
         ALLOCATE (V_guess_vec(n*n))
         V_guess(1:n, 1:n) => V_guess_vec
         CALL pao_guess_initial_potential(qs_env, iatom, V_guess)

         ! build overlap matrix
         ALLOCATE (S(nterms, nterms))
         S(:, :) = MATMUL(TRANSPOSE(V_terms), V_terms)

         ! diagonalize S
         ALLOCATE (S_evals(nterms), S_evecs(nterms, nterms))
         S_evecs(:, :) = S
         CALL diamat_all(S_evecs, S_evals)

         ! calculate Tikhonov regularized inverse
         ALLOCATE (S_inv(nterms, nterms))
         S_inv(:, :) = 0.0_dp
         DO k = 1, nterms
            w = S_evals(k)/(S_evals(k)**2 + pao%linpot_init_delta)
            DO i = 1, nterms
            DO j = 1, nterms
               S_inv(i, j) = S_inv(i, j) + w*S_evecs(i, k)*S_evecs(j, k)
            END DO
            END DO
         END DO

         ! perform fit
         block_X(:, 1) = MATMUL(MATMUL(S_inv, TRANSPOSE(V_terms)), V_guess_vec)

         ! clean up
         DEALLOCATE (V_guess_vec, S, S_evecs, S_evals, S_inv)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE pao_param_initguess_linpot

END MODULE pao_param_linpot

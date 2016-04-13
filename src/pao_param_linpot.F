!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common framework for a linear parametrization of the potential.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_linpot
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_col_block_sizes, cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_frobenius_norm, &
        cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
        cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, &
        cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, cp_dbcsr_type, &
        dbcsr_distribution_mp, dbcsr_mp_group
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_min,&
                                              mp_sum
   USE pao_input,                       ONLY: pao_fock_param,&
                                              pao_rotinv_param
   USE pao_linpot_full,                 ONLY: linpot_full_calc_term,&
                                              linpot_full_count_terms
   USE pao_linpot_rotinv,               ONLY: linpot_rotinv_calc_term,&
                                              linpot_rotinv_count_terms
   USE pao_param_fock,                  ONLY: pao_calc_U_block_fock,&
                                              pao_param_finalize_fock,&
                                              pao_param_init_fock
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_param_init_linpot, pao_param_finalize_linpot, pao_calc_U_linpot
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

      INTEGER                                            :: acol, arow, handle, iatom, ikind, kterm, &
                                                            n, natoms, nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pri, col_blk_size, row_blk_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_V_term, block_V_terms
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      CALL pao_param_init_fock(pao, qs_env)

      CALL get_qs_env(qs_env, &
                      matrix_s=matrix_s, &
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
      ENDDO

      ! allocate matrix_V_terms
      blk_sizes_pri => cp_dbcsr_col_block_sizes(matrix_s(1)%matrix)
      row_blk_size = blk_sizes_pri**2
      CALL cp_dbcsr_init(pao%matrix_V_terms)
      CALL cp_dbcsr_create(pao%matrix_V_terms, &
                           template=matrix_s(1)%matrix, &
                           matrix_type="N", &
                           row_blk_size=row_blk_size, &
                           col_blk_size=col_blk_size)
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_V_terms)
      DEALLOCATE (row_blk_size, col_blk_size)

      ! calculate, normalize, and store potential terms as rows of block_V_terms
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_V_terms)
         iatom = arow
         n = blk_sizes_pri(iatom)
         CPASSERT(n*n == SIZE(block_V_terms, 1))
         nterms = SIZE(block_V_terms, 2)
         DO kterm = 1, nterms
            block_V_term(1:n, 1:n) => block_V_terms(:, kterm) ! map column into matrix
            CALL linpot_calc_term(pao, qs_env, iatom, kterm, block_V_term)
            block_V_term = block_V_term/SQRT(SUM(block_V_term*block_V_term)) ! normalize
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      IF (pao%precondition) &
         CPABORT("PAO preconditioning not supported for selected parametrization.")

      CALL timestop(handle)
   END SUBROUTINE pao_param_init_linpot

! **************************************************************************************************
!> \brief Finalize the linear potential parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_linpot(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL pao_param_finalize_fock(pao)
      CALL cp_dbcsr_release(pao%matrix_V_terms)

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

      SELECT CASE (pao%parameterization)
      CASE (pao_fock_param)
         CALL linpot_full_count_terms(qs_env, ikind, nterms=nparams)
      CASE (pao_rotinv_param)
         CALL linpot_rotinv_count_terms(qs_env, ikind, nterms=nparams)
      CASE DEFAULT
         CPABORT("unkown parameterization")
      END SELECT

   END SUBROUTINE pao_param_count_linpot

! **************************************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param penalty ...
!> \param matrix_M ...
!> \param matrix_G ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U_linpot(pao, penalty, matrix_M, matrix_G)
      TYPE(pao_env_type), POINTER                        :: pao
      REAL(dp), INTENT(INOUT)                            :: penalty
      TYPE(cp_dbcsr_type), OPTIONAL                      :: matrix_M, matrix_G

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U_linpot'

      INTEGER                                            :: acol, arow, group, handle, iatom, n
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:), POINTER                    :: vec_M2, vec_V
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_M1, block_M2, block_U, &
                                                            block_V, block_V_terms, block_X
      REAL(KIND=dp)                                      :: min_gap
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      min_gap = HUGE(dp)
      group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(pao%matrix_U)))

      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         IF (arow /= acol) CPABORT("pao_calc_U_linpot: encountered off-diagonal block")
         iatom = arow
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=iatom, col=iatom, block=block_U, found=found)
         CPASSERT(ASSOCIATED(block_U))
         n = SIZE(block_U, 1)

         ! calculate potential V
         ALLOCATE (vec_V(n*n))
         vec_V(:) = 0.0_dp
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=block_V_terms, found=found)
         IF (SIZE(block_V_terms) > 0) & ! protect against corner-case of zero pao parameters
            vec_V = MATMUL(block_V_terms, block_X(:, 1))
         block_V(1:n, 1:n) => vec_V(:) ! map vector into matrix

         ! symmetrize
         IF (MAXVAL(ABS(block_V-TRANSPOSE(block_V))/MAX(1.0_dp, MAXVAL(ABS(block_V)))) > 1e-12) &
            CPABORT("block_V not symmetric")
         block_V = 0.5_dp*(block_V+TRANSPOSE(block_V)) ! symmetrize exactly

         IF (.NOT. PRESENT(matrix_G)) THEN
            CALL pao_calc_U_block_fock(pao, iatom, penalty, min_gap, block_V, block_U)

         ELSE ! TURNING POINT (if calc grad) ------------------------------------
            CPASSERT(PRESENT(matrix_M))
            CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
            CPASSERT(ASSOCIATED(block_G))
            CALL cp_dbcsr_get_block_p(matrix=matrix_M, row=iatom, col=iatom, block=block_M1, found=found)
            ! don't check ASSOCIATED(block_M1), it might have been filtered out.

            ALLOCATE (vec_M2(n*n))
            block_M2(1:n, 1:n) => vec_M2(:) ! map vector into matrix
            CALL pao_calc_U_block_fock(pao, iatom, penalty, min_gap, block_V, block_U, block_M1, block_M2)
            IF (MAXVAL(ABS(block_M2-TRANSPOSE(block_M2))) > 1e-14_dp) CPABORT("Expect symmetric matrix")
            IF (SIZE(block_V_terms) > 0) & ! protect against corner-case of zero pao parameters
               block_G(:, 1) = MATMUL(vec_M2, block_V_terms)
            DEALLOCATE (vec_M2)
         ENDIF

         DEALLOCATE (vec_V)

      END DO
      CALL cp_dbcsr_iterator_stop(iter)

      ! sum penalty energies across ranks
      CALL mp_sum(penalty, group)

      ! print smallest homo-lumo gap encountered by fock-layer
      CALL mp_min(min_gap, group)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| min_gap:", min_gap

      CALL timestop(handle)

   END SUBROUTINE pao_calc_U_linpot

! **************************************************************************************************
!> \brief Internal routine, calculates terms in potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param kterm ...
!> \param block_V ...
! **************************************************************************************************
   SUBROUTINE linpot_calc_term(pao, qs_env, iatom, kterm, block_V)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom, kterm
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V

      SELECT CASE (pao%parameterization)
      CASE (pao_fock_param)
         CALL linpot_full_calc_term(kterm, block_V)
      CASE (pao_rotinv_param)
         CALL linpot_rotinv_calc_term(qs_env, iatom, kterm, block_V)
      CASE DEFAULT
         CPABORT("unkown parameterization")
      END SELECT

   END SUBROUTINE linpot_calc_term

! **************************************************************************************************
!> \brief Calculate intial guess for matrix_X
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_linpot(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initguess_linpot'

      INTEGER                                            :: acol, arow, handle, n
      LOGICAL                                            :: found
      REAL(dp)                                           :: norm
      REAL(dp), DIMENSION(:), POINTER                    :: V_guess_vec
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_H0, block_ks, block_V_terms, &
                                                            block_X, V_guess
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)

      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)

         ! retrieve potential terms
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=arow, col=acol, block=block_V_terms, found=found)
         CPASSERT(ASSOCIATED(block_V_terms))
         IF (SIZE(block_V_terms) == 0) CYCLE ! protect against corner-case of zero pao parameters

         ! retrieve block_ks and block_H0
         CALL cp_dbcsr_get_block_p(matrix=matrix_ks(1)%matrix, row=arow, col=acol, block=block_ks, found=found)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_H0, row=arow, col=acol, block=block_H0, found=found)
         CPASSERT(ASSOCIATED(block_H0) .AND. ASSOCIATED(block_ks))

         ! guess a potential
         n = SIZE(block_ks, 1)
         ALLOCATE (V_guess_vec(n*n))
         V_guess(1:n, 1:n) => V_guess_vec
         V_guess = block_ks-block_H0

         ! project out anything that the potential terms can not represent
         V_guess_vec = MATMUL(MATMUL(block_V_terms, TRANSPOSE(block_V_terms)), V_guess_vec)

         ! find least-square-fit of guessed potential
         CALL dgelsy_wrapper(block_V_terms, V_guess_vec, block_X(:, 1))

         DEALLOCATE (V_guess_vec)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      norm = cp_dbcsr_frobenius_norm(pao%matrix_X)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Calculated initial guess, norm(X):", norm

      CALL timestop(handle)
   END SUBROUTINE pao_param_initguess_linpot

! **************************************************************************************************
!> \brief Helper routine for calling LAPACK dgelsy
!>  computes the linear least square fit min|| A * X - B ||
!> \param matrix_A possibly non-square matrix
!> \param vec_B ...
!> \param vec_X ...
! **************************************************************************************************
   SUBROUTINE dgelsy_wrapper(matrix_A, vec_B, vec_X)
      REAL(dp), DIMENSION(:, :)                          :: matrix_A
      REAL(dp), DIMENSION(:)                             :: vec_B, vec_X

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dgelsy_wrapper'

      INTEGER                                            :: handle, info, lda, ldb, lwork, m, n, &
                                                            nrhs, rank
      INTEGER, DIMENSION(:), POINTER                     :: jpvt
      REAL(KIND=dp)                                      :: rcon
      REAL(KIND=dp), DIMENSION(:), POINTER               :: B, work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: A

      CALL timeset(routineN, handle)
      NULLIFY (jpvt, work, B)

      IF (SIZE(matrix_A, 1) /= SIZE(vec_B)) CPABORT("DGELSY: matrix_A and vec_B not compatible")
      IF (SIZE(matrix_A, 2) /= SIZE(vec_X)) CPABORT("DGELSY: matrix_A and vec_X not compatible")
      m = SIZE(matrix_A, 1)
      n = SIZE(matrix_A, 2)
      lda = m
      ldb = MAX(m, n)
      nrhs = 1 ! number of right hand sides

      rcon = 0.0_dp ! est. condition number

      ! optional permutations of A
      ALLOCATE (jpvt(n))
      jpvt(:) = 0 ! no permutation

      ! make copies of matrix_A and vec_B
      ALLOCATE (A(m, n), B(MAX(m, n)))
      B(1:m) = vec_B(:)
      A(:, :) = matrix_A(:, :)

      ! query required work space
      lwork = -1
      ALLOCATE (work(1))
      CALL DGELSY(m, n, nrhs, A(1, 1), lda, B(1), ldb, jpvt, rcon, rank, work(1), lwork, info)
      IF (info /= 0) CPABORT("DGELSY work space query failed")

      ! allocate work space
      lwork = INT(REAL(work(1)))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))

      ! perform actual optimization
      CALL DGELSY(m, n, nrhs, A(1, 1), lda, B(1), ldb, jpvt, rcon, rank, work(1), lwork, info)
      IF (info /= 0) CPABORT("DGELSY optimization failed")

      ! copy results
      vec_X(:) = B(1:n)

      DEALLOCATE (A, B, work, jpvt)
      CALL timestop(handle)

   END SUBROUTINE dgelsy_wrapper

END MODULE pao_param_linpot

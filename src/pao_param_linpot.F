!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Common framework for a linear parametrization of the potential.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_param_linpot
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, cp_dbcsr_scale, &
       cp_dbcsr_set
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_sum
  USE pao_input,                       ONLY: pao_fock_param,&
                                             pao_none_precon,&
                                             pao_rotinv_param,&
                                             pao_sinv_precon
  USE pao_linpot_full,                 ONLY: linpot_full_calc_term,&
                                             linpot_full_count_terms
  USE pao_linpot_rotinv,               ONLY: linpot_rotinv_calc_term,&
                                             linpot_rotinv_count_terms
  USE pao_param_fock,                  ONLY: pao_calc_U_fock,&
                                             pao_param_finalize_fock,&
                                             pao_param_init_fock
  USE pao_types,                       ONLY: pao_env_type
  USE pao_utils,                       ONLY: fill_randomly
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_linpot'

  PUBLIC :: pao_param_init_linpot, pao_param_finalize_linpot, pao_calc_U_linpot

CONTAINS

! *****************************************************************************
!> \brief Initialize the linear potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param reuse_matrix_X ...
! *****************************************************************************
  SUBROUTINE pao_param_init_linpot(pao, qs_env, reuse_matrix_X)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: reuse_matrix_X

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_linpot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, iatom, &
                                                kterm, max_terms, natoms, &
                                                nterms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(dp), DIMENSION(:, :), POINTER       :: block_V_term
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)
    CALL pao_param_init_fock(pao, qs_env)

    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    matrix_ks=matrix_ks,&
                    dft_control=dft_control,&
                    natom=natoms)

    IF(dft_control%nspins /= 1) CPABORT("open shell not yet implemented")

    ALLOCATE(row_blk_size(natoms), col_blk_size(natoms))

    ! figure out number of potential terms
    max_terms = 0
    DO iatom=1, natoms
       CALL linpot_count_terms(pao, qs_env, iatom, nterms)
       max_terms = MAX(max_terms, nterms)
       col_blk_size(iatom) = 1 ! a skinny matrix
       row_blk_size(iatom) = nterms
    ENDDO

    IF(reuse_matrix_X) THEN
       IF(pao%iw>0) WRITE(pao%iw,*) "PAO| reusing potential from last MD step"
    ELSE
       ! allocate matrix_X
       CALL cp_dbcsr_init(pao%matrix_X)
       CALL cp_dbcsr_create(pao%matrix_X,&
                            template=matrix_s(1)%matrix,&
                            matrix_type="N",&
                            row_blk_size=row_blk_size,&
                            col_blk_size=col_blk_size)
       CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_X)

       ! init matrix_X
       IF(pao%iw>0) WRITE(pao%iw,*) "PAO| initializing potential with random noise from 0.0 ... +0.1"
       CALL fill_randomly(pao%matrix_X)
       CALL cp_dbcsr_scale(pao%matrix_X, +0.1_dp)
    ENDIF

    ! calculate, normalize, and store potential terms
    ALLOCATE(pao%matrix_V_terms(max_terms))
    DO kterm=1, max_terms
       CALL cp_dbcsr_init(pao%matrix_V_terms(kterm))
       CALL cp_dbcsr_create(pao%matrix_V_terms(kterm), template=matrix_s(1)%matrix)
       CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_V_terms(kterm))
       CALL cp_dbcsr_set(pao%matrix_V_terms(kterm), 0.0_dp)
       CALL cp_dbcsr_iterator_start(iter, pao%matrix_V_terms(kterm))
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_V_term)
          IF(arow /= acol) CPABORT("encountered off-diagonal block")
          iatom = arow
          CALL linpot_count_terms(pao, qs_env, iatom, nterms)
          IF(kterm <= nterms) THEN
            CALL linpot_calc_term(pao, qs_env, iatom, kterm, block_V_term)
            block_V_term = block_V_term / SQRT(SUM(block_V_term*block_V_term))
          ENDIF
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDDO

    ! build preconditioner
    SELECT CASE(pao%preconditioner)
      CASE(pao_none_precon)
        CALL build_preconditioner_none(pao, qs_env)
      CASE(pao_sinv_precon)
        CALL build_preconditioner_Sinv(pao, qs_env)
      CASE DEFAULT
        CPABORT("unkown preconditioner")
    END SELECT

    DEALLOCATE(row_blk_size, col_blk_size)
    CALL timestop(handle)
  END SUBROUTINE pao_param_init_linpot


! *****************************************************************************
!> \brief Finalize the linear potential parametrization
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_param_finalize_linpot(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize_linpot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iterm

    CALL pao_param_finalize_fock(pao)

    DO iterm=1, SIZE(pao%matrix_V_terms)
       CALL cp_dbcsr_release(pao%matrix_V_terms(iterm))
    ENDDO
    DEALLOCATE(pao%matrix_V_terms)

    CALL cp_dbcsr_release(pao%matrix_precon)
    CALL cp_dbcsr_release(pao%matrix_precon_inv)
    ! pao%matrix_X is keept for the next MD step

  END SUBROUTINE pao_param_finalize_linpot


! *****************************************************************************
!> \brief Construct the identity preconditioner
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
SUBROUTINE build_preconditioner_none(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'build_preconditioner_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natoms, nterms
    INTEGER, DIMENSION(:), POINTER           :: blk_size

    CALL get_qs_env(qs_env, natom=natoms)

    ! collecte nterms for all atoms
    ALLOCATE(blk_size(natoms))
    DO iatom=1, natoms
       CALL linpot_count_terms(pao, qs_env, iatom, nterms)
       blk_size(iatom) = nterms
    ENDDO

    ! create preconditioner matrices
    CALL cp_dbcsr_init(pao%matrix_precon)
    CALL cp_dbcsr_create(pao%matrix_precon,&
                         template=pao%matrix_X,&
                         matrix_type="N",&
                         row_blk_size=blk_size,&
                         col_blk_size=blk_size)
    DEALLOCATE(blk_size)
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_precon)
    CALL cp_dbcsr_set(pao%matrix_precon, 0.0_dp)
    CALL cp_dbcsr_add_on_diag(pao%matrix_precon, 1.0_dp)
    CALL cp_dbcsr_init(pao%matrix_precon_inv)
    CALL cp_dbcsr_copy(pao%matrix_precon_inv, pao%matrix_precon)

  END SUBROUTINE build_preconditioner_none


! *****************************************************************************
!> \brief Construct preconditioner based on invers overlap of potential terms.
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE build_preconditioner_Sinv(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'build_preconditioner_Sinv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, i, iatom, j, k, &
                                                natoms, nfiltered, nterms
    INTEGER, DIMENSION(:), POINTER           :: blk_size
    LOGICAL                                  :: found
    REAL(dp)                                 :: s, v, w
    REAL(dp), DIMENSION(:), POINTER          :: evals
    REAL(dp), DIMENSION(:, :), POINTER       :: A, block_precon, &
                                                block_precon_inv, term_i, &
                                                term_j
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL get_qs_env(qs_env, natom=natoms, para_env=para_env)

    ! collecte nterms for all atoms
    ALLOCATE(blk_size(natoms))
    DO iatom=1, natoms
       CALL linpot_count_terms(pao, qs_env, iatom, nterms)
       blk_size(iatom) = nterms
    ENDDO

    ! create preconditioner matrices
    CALL cp_dbcsr_init(pao%matrix_precon)
    CALL cp_dbcsr_create(pao%matrix_precon,&
                         template=pao%matrix_X,&
                         matrix_type="N",&
                         row_blk_size=blk_size,&
                         col_blk_size=blk_size)
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_precon)
    CALL cp_dbcsr_set(pao%matrix_precon, 0.0_dp)
    CALL cp_dbcsr_init(pao%matrix_precon_inv)
    CALL cp_dbcsr_copy(pao%matrix_precon_inv, pao%matrix_precon)

    nfiltered = 0
    CALL cp_dbcsr_iterator_start(iter, pao%matrix_precon)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_precon)
       iatom = arow
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_precon_inv, row=iatom, col=iatom, block=block_precon_inv, found=found)
       CPASSERT(ASSOCIATED(block_precon_inv))

       ! calculate overlap matrix between potential terms
       nterms = SIZE(block_precon, 1)
       ALLOCATE(A(nterms,nterms), evals(nterms))
       DO i=1, nterms
       DO j=1, nterms
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms(i), row=iatom, col=iatom, block=term_i, found=found)
         CPASSERT(ASSOCIATED(term_i))
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms(j), row=iatom, col=iatom, block=term_j, found=found)
         CPASSERT(ASSOCIATED(term_j))
         A(i,j) = SUM(term_i * term_j)
       ENDDO
       ENDDO

       CALL diamat_all(A, evals)
       block_precon = 0.0_dp
       block_precon_inv = 0.0_dp
       DO k=1, nterms
          !IF(ABS(evals(k)) < 1e-12_dp) THEN
          IF(ABS(evals(k)) < 1e-5_dp) THEN ! we ignore tiny eigenvalue
             nfiltered = nfiltered + 1
          ELSE
             v = 1.0_dp/evals(k)
             w = 1.0_dp/v
             DO i=1, nterms
             DO j=1, nterms
                block_precon(i,j)     = block_precon(i,j)     + v*A(i,k)*A(j,k)
                block_precon_inv(i,j) = block_precon_inv(i,j) + w*A(i,k)*A(j,k)
             ENDDO
             ENDDO
          ENDIF
       ENDDO

       ! scale to get line-search steps of order 1
       !TODO: Is this really the best scaling factor?
       s = 0.5_dp * REAL(nterms, dp)
       block_precon = block_precon / s
       block_precon_inv = block_precon_inv * s
       DEALLOCATE(A, evals)

    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    ! print statistics
    CALL mp_sum(nfiltered, para_env%group)
    IF(pao%iw>0) WRITE (pao%iw,*) "PAO| SINV-preconditioner nterm:",SUM(blk_size),&
                                  "nfiltered: ", nfiltered,&
                                  "nterms left:", SUM(blk_size) - nfiltered
    DEALLOCATE(blk_size)

  END SUBROUTINE build_preconditioner_Sinv


! *****************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param iatom ...
!> \param block_X ...
!> \param block_U ...
!> \param block_M1 ...
!> \param block_G ...
! *****************************************************************************
  SUBROUTINE pao_calc_U_linpot(pao, iatom, block_X, block_U, block_M1, block_G)
    TYPE(pao_env_type), POINTER              :: pao
    INTEGER, INTENT(IN)                      :: iatom
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, block_U
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: block_M1, block_G

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_linpot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, kterm, N, nterms
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:), POINTER          :: vec_G, vec_X
    REAL(dp), DIMENSION(:, :), POINTER       :: block_M2, block_V, &
                                                block_V_term

    CALL timeset(routineN,handle)
    N = SIZE(block_U, 1)
    nterms = SIZE(block_X)
    ALLOCATE(block_V(N,N), vec_X(nterms))

    vec_X = block_X(:, 1)

    ! calc block_V
    block_V(:,:) = 0.0_dp
    DO kterm=1, nterms
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms(kterm), row=iatom, col=iatom, block=block_V_term, found=found)
       CPASSERT(ASSOCIATED(block_V_term))
       block_V = block_V + vec_X(kterm) * block_V_term
    ENDDO

    IF(MAXVAL(ABS(block_V-TRANSPOSE(block_V))/MAX(1.0_dp,MAXVAL(ABS(block_V)))) > 1e-12)&
        CPABORT("block_V not symmetric")
    block_V = 0.5_dp*( block_V + TRANSPOSE(block_V) ) ! symmetrize exactly

    IF(.NOT. PRESENT(block_G))THEN
       CALL pao_calc_U_fock(pao, iatom, block_V, block_U)

    ELSE  ! TURNING POINT (if calc grad) --------------------------------------

       CPASSERT(PRESENT(block_M1))
       ALLOCATE(block_M2(N,N), vec_G(nterms))
       CALL pao_calc_U_fock(pao, iatom, block_V, block_U, block_M1, block_M2)

       DO kterm=1, nterms
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms(kterm), row=iatom, col=iatom, block=block_V_term, found=found)
          CPASSERT(ASSOCIATED(block_V_term))
          vec_G(kterm) = SUM(block_V_term * block_M2)
       ENDDO

       block_G(:,1) = vec_G

       DEALLOCATE(block_M2, vec_G)
    ENDIF

    DEALLOCATE(block_V, vec_X)
    CALL timestop(handle)

  END SUBROUTINE pao_calc_U_linpot


! *****************************************************************************
!> \brief Internal routine, counts number of terms in potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param nterms ...
! *****************************************************************************
  SUBROUTINE linpot_count_terms(pao, qs_env, iatom, nterms)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    INTEGER, INTENT(OUT)                     :: nterms

    CHARACTER(len=*), PARAMETER :: routineN = 'linpot_count_terms', &
      routineP = moduleN//':'//routineN

     SELECT CASE(pao%parameterization)
      CASE(pao_fock_param)
        CALL linpot_full_count_terms(qs_env, iatom, nterms)
      CASE(pao_rotinv_param)
        CALL linpot_rotinv_count_terms(qs_env, iatom, nterms)
      CASE DEFAULT
        CPABORT("unkown parameterization")
     END SELECT

  END SUBROUTINE linpot_count_terms


! *****************************************************************************
!> \brief Internal routine, calculates terms in potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param kterm ...
!> \param block_V ...
! *****************************************************************************
  SUBROUTINE linpot_calc_term(pao, qs_env, iatom, kterm, block_V)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom, kterm
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: block_V

    CHARACTER(len=*), PARAMETER :: routineN = 'linpot_calc_term', &
      routineP = moduleN//':'//routineN

     SELECT CASE(pao%parameterization)
      CASE(pao_fock_param)
        CALL linpot_full_calc_term(kterm, block_V)
      CASE(pao_rotinv_param)
        CALL linpot_rotinv_calc_term(qs_env, iatom, kterm, block_V)
      CASE DEFAULT
        CPABORT("unkown parameterization")
     END SELECT

  END SUBROUTINE linpot_calc_term

END MODULE pao_param_linpot

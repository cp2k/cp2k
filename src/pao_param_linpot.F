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
       cp_dbcsr_create, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_reserve_diag_blocks, cp_dbcsr_set
  USE kinds,                           ONLY: dp
  USE pao_input,                       ONLY: pao_fock_param,&
                                             pao_rotinv_param
  USE pao_linpot_full,                 ONLY: linpot_full_calc_term,&
                                             linpot_full_count_terms
  USE pao_linpot_rotinv,               ONLY: linpot_rotinv_calc_term,&
                                             linpot_rotinv_count_terms
  USE pao_param_fock,                  ONLY: pao_calc_U_fock,&
                                             pao_param_finalize_fock,&
                                             pao_param_init_fock
  USE pao_types,                       ONLY: pao_env_type
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
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_H, block_V_term, block_X
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
       CALL cp_dbcsr_set(pao%matrix_X, 0.0_dp)
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
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=arow, col=acol, block=block_X, found=found)
          CALL cp_dbcsr_get_block_p(matrix=matrix_ks(1)%matrix, row=arow, col=acol, block=block_H, found=found)
          CPASSERT(ASSOCIATED(block_X) .AND. ASSOCIATED(block_H))
          IF(arow /= acol) CPABORT("encountered off-diagonal block")
          iatom = arow
          CALL linpot_count_terms(pao, qs_env, iatom, nterms)
          IF(kterm <= nterms) THEN
            CALL linpot_calc_term(pao, qs_env, iatom, kterm, block_V_term)
            block_V_term = block_V_term / SQRT(SUM(block_V_term*block_V_term))
            block_X(kterm,1) = SUM(block_V_term * block_H) ! seems like a good initial guess
          ENDIF
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDDO

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

    ! pao%matrix_X is kept for the next MD step

  END SUBROUTINE pao_param_finalize_linpot


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

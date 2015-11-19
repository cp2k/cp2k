!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Main module for the PAO method
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_main
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_frobenius_norm, cp_dbcsr_get_block_p, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_type
  USE dm_ls_scf_qs,                    ONLY: ls_scf_dm_to_ks,&
                                             ls_scf_qs_atomic_guess
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE linesearch,                      ONLY: linesearch_finalize,&
                                             linesearch_init,&
                                             linesearch_step
  USE machine,                         ONLY: m_walltime
  USE mathlib,                         ONLY: diamat_all
  USE pao_input,                       ONLY: parse_pao_section
  USE pao_io,                          ONLY: pao_read_restart,&
                                             pao_write_restart
  USE pao_param,                       ONLY: pao_calc_grad,&
                                             pao_param_finalize,&
                                             pao_param_init,&
                                             pao_update_AB
  USE pao_preopt,                      ONLY: pao_preopt_run
  USE pao_types,                       ONLY: pao_env_type
  USE pao_utils,                       ONLY: pao_calc_CG_beta,&
                                             pao_print_step_delimiter,&
                                             pao_test_convergence
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_main'

  PUBLIC :: pao_init, pao_scf_init, pao_update, pao_finalize

CONTAINS

! *****************************************************************************
!> \brief Initialize the PAO environment
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_init(qs_env,ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(section_vals_type), POINTER         :: input, section

    IF(.NOT. ls_scf_env%do_pao) RETURN

    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,&
                    input=input)

    pao => ls_scf_env%pao_env
    ls_mstruct => ls_scf_env%ls_mstruct

    CALL parse_pao_section(pao, input)

    ! preconditioner has to be build before pao_param_init()
    CALL pao_build_orthogonalizer(pao, qs_env)

    IF(pao%istep>0) THEN ! was there a previous MD step?
       CALL pao_param_init(pao, qs_env, reuse_matrix_X=.TRUE.)
    ELSE
       CALL pao_param_init(pao, qs_env, reuse_matrix_X=.FALSE.)
       CALL pao_read_restart(pao, qs_env)
    ENDIF

    ! init counters
    pao%step_start_time = m_walltime()
    pao%istep = 0
    pao%last_reset = 0

    CALL cp_dbcsr_init(pao%matrix_X_orig)

    CALL cp_dbcsr_init(pao%matrix_G)
    CALL cp_dbcsr_copy(pao%matrix_G, pao%matrix_X)
    CALL cp_dbcsr_set(pao%matrix_G, 0.0_dp)

    CALL cp_dbcsr_init(pao%matrix_G_prev)
    CALL cp_dbcsr_copy(pao%matrix_G_prev, pao%matrix_G)

    CALL cp_dbcsr_init(pao%matrix_D)
    CALL cp_dbcsr_copy(pao%matrix_D, pao%matrix_G)

    CALL pao_build_selector(qs_env, pao%matrix_Y)

    CALL print_atom_info(pao)

    CALL cp_dbcsr_init(pao%matrix_U)
    CALL cp_dbcsr_create(pao%matrix_U, template=matrix_s(1)%matrix, matrix_type="N", name="PAO matrix_U")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_U)

    CALL cp_dbcsr_init(ls_mstruct%matrix_A)
    CALL cp_dbcsr_create(ls_mstruct%matrix_A, template=pao%matrix_Y)

    CALL cp_dbcsr_init(ls_mstruct%matrix_B)
    CALL cp_dbcsr_create(ls_mstruct%matrix_B, template=pao%matrix_Y)

    CALL pao_update_AB(pao, ls_mstruct)

    section => section_vals_get_subs_vals(input,"DFT%LS_SCF%PAO%LINE_SEARCH")
    CALL linesearch_init(pao%linesearch, section, "PAO|")

    CALL pao_preopt_run(pao, qs_env, ls_scf_env)

    CALL timestop(handle)
  END SUBROUTINE pao_init


! *****************************************************************************
!> \brief Prints a one line summary for each atom.
!> \param pao ...
! *****************************************************************************
  SUBROUTINE print_atom_info(pao)
    TYPE(pao_env_type), POINTER              :: pao

    CHARACTER(len=*), PARAMETER :: routineN = 'print_atom_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natoms
    INTEGER, DIMENSION(:), POINTER           :: pao_basis, param_cols, &
                                                param_rows, pri_basis

    pri_basis => cp_dbcsr_row_block_sizes(pao%matrix_Y)
    pao_basis => cp_dbcsr_col_block_sizes(pao%matrix_Y)

    CPASSERT(SIZE(pao_basis) == SIZE(pri_basis))
    natoms = SIZE(pao_basis)

    param_rows => cp_dbcsr_row_block_sizes(pao%matrix_X)
    param_cols => cp_dbcsr_col_block_sizes(pao%matrix_X)
    CPASSERT(SIZE(param_rows)==natoms.AND.SIZE(param_cols)==natoms)

    DO iatom=1, natoms
      IF(pao%iw>0) WRITE (pao%iw,"(A,I7,T20,A,I3,T45,A,I3,T65,A,I3)")&
                         " PAO| atom: ", iatom,&
                         " prim_basis: ", pri_basis(iatom),&
                         " pao_basis: ", pao_basis(iatom),&
                         " pao_params: ", (param_cols(iatom)*param_rows(iatom))
    ENDDO
  END SUBROUTINE print_atom_info


! *****************************************************************************
!> \brief Signal PAO that a new SCF optimization is about to start
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_scf_init(qs_env,ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_scf_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(cp_dbcsr_type)                      :: matrix_tmp, matrix_tmp2, &
                                                matrix_Z
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao

    IF(.NOT. ls_scf_env%do_pao) RETURN

    CALL timeset(routineN,handle)
    pao => ls_scf_env%pao_env
    ls_mstruct => ls_scf_env%ls_mstruct

    pao%istep = pao%istep + 1
    CALL pao_print_step_delimiter(pao)

    IF(.NOT.ASSOCIATED(pao%matrix_P_guess)) THEN
       CALL ls_scf_qs_atomic_guess(qs_env,ls_scf_env%energy_init)
       IF(pao%iw>0) WRITE(pao%iw,'(A,F20.9)') " PAO| Energy from initial atomic guess:",ls_scf_env%energy_init
    ELSE
       !TODO: use ASPC if run as MD
       pao => ls_scf_env%pao_env
       ls_mstruct => ls_scf_env%ls_mstruct

       CALL cp_dbcsr_init(matrix_tmp)
       CALL cp_dbcsr_create(matrix_tmp, template=ls_scf_env%ls_mstruct%matrix_B,matrix_type="N")

       CALL cp_dbcsr_init(matrix_Z)
       CALL cp_dbcsr_create(matrix_Z, template=ls_scf_env%matrix_p(1),matrix_type="N")

       CALL cp_dbcsr_init(matrix_tmp2)
       CALL cp_dbcsr_create(matrix_tmp2, template=ls_scf_env%matrix_p(1),matrix_type="N")

       DO ispin=1, ls_scf_env%nspins
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_P_guess(ispin), ls_scf_env%ls_mstruct%matrix_B, &
                                 0.0_dp, matrix_tmp)
          CALL cp_dbcsr_multiply("T", "N", 1.0_dp, ls_scf_env%ls_mstruct%matrix_B, matrix_tmp, &
                                 0.0_dp, matrix_Z)

          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_Z, ls_scf_env%matrix_s_inv,&
                                 0.0_dp, matrix_tmp2)
          CALL cp_dbcsr_multiply("T", "N", 1.0_dp, ls_scf_env%matrix_s_inv, matrix_tmp2,&
                                 0.0_dp, ls_scf_env%matrix_p(ispin))
       ENDDO

       CALL cp_dbcsr_release(matrix_Z)
       CALL cp_dbcsr_release(matrix_tmp)
       CALL cp_dbcsr_release(matrix_tmp2)

       ! compute corresponding energy and ks matrix
       CALL ls_scf_dm_to_ks(qs_env,ls_scf_env,ls_scf_env%energy_init,iscf=0)
       IF(pao%iw>0) WRITE (pao%iw,*) "PAO| Initial energy from stored density:",ls_scf_env%energy_init
    ENDIF

    CALL timestop(handle)
  END SUBROUTINE pao_scf_init


! *****************************************************************************
!> \brief Constructs matrix_N and its inverse.
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_build_orthogonalizer(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_orthogonalizer', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, iatom, &
                                                j, k, N, natoms
    LOGICAL                                  :: found
    REAL(dp)                                 :: v, w
    REAL(dp), DIMENSION(:), POINTER          :: evals
    REAL(dp), DIMENSION(:, :), POINTER       :: A, block_N, block_N_inv, &
                                                block_S
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    natom=natoms)

    CALL cp_dbcsr_init(pao%matrix_N)
    CALL cp_dbcsr_create(pao%matrix_N, template=matrix_s(1)%matrix,  name="PAO matrix_N")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N)

    CALL cp_dbcsr_init(pao%matrix_N_inv)
    CALL cp_dbcsr_create(pao%matrix_N_inv, template=matrix_s(1)%matrix, name="PAO matrix_N_inv")
    CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N_inv)

    CALL cp_dbcsr_iterator_start(iter, pao%matrix_N)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_N)
       IF(arow /= acol) CPABORT("got non-diagonal block")
       iatom = arow

       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=arow, col=acol, block=block_N_inv, found=found)
       CPASSERT(ASSOCIATED(block_N_inv))

       CALL cp_dbcsr_get_block_p(matrix=matrix_s(1)%matrix, row=arow, col=acol, block=block_S, found=found)
       CPASSERT(ASSOCIATED(block_S))

       IF(SIZE(block_S,1)/=SIZE(block_S,2)) CPABORT("Expected symmetric block")
       N = SIZE(block_S,1)
       ALLOCATE(A(N,N), evals(N))

       ! take square root of atomic overlap matrix
       A(:,:) = block_S(:,:)
       CALL diamat_all(A, evals) !afterwards A contains the eigenvectors
       DO k=1, N
          ! NOTE: To maintain a consistent notation with the Berghold paper,
          ! the "_inv" is swapped: N^{-1}=sqrt(S); N=sqrt(S)^{-1}
          w = 1.0_dp / SQRT(evals(k))
          v = SQRT(evals(k))
          DO i=1, N
             DO j=1, N
                block_N(i,j)     = block_N(i,j)     + w*A(i,k)*A(j,k)
                block_N_inv(i,j) = block_N_inv(i,j) + v*A(i,k)*A(j,k)
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE(A, evals)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE pao_build_orthogonalizer


! *****************************************************************************
!> \brief Finalize the PAO environment
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_finalize(ls_scf_env)
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao

    IF(.NOT. ls_scf_env%do_pao) RETURN

    CALL timeset(routineN,handle)
    pao => ls_scf_env%pao_env
    ls_mstruct => ls_scf_env%ls_mstruct

    CALL pao_param_finalize(pao)

    ! We keep pao%matrix_X for next scf-run, e.g. during MD or GEO-OPT
    CALL cp_dbcsr_release(pao%matrix_X_orig)
    CALL cp_dbcsr_release(pao%matrix_G)
    CALL cp_dbcsr_release(pao%matrix_G_prev)
    CALL cp_dbcsr_release(pao%matrix_D)
    CALL cp_dbcsr_release(pao%matrix_U)
    CALL cp_dbcsr_release(pao%matrix_Y)
    CALL cp_dbcsr_release(pao%matrix_N)
    CALL cp_dbcsr_release(pao%matrix_N_inv)

    CALL cp_dbcsr_release(ls_mstruct%matrix_A)
    CALL cp_dbcsr_release(ls_mstruct%matrix_B)

    IF(ASSOCIATED(pao%matrix_P_guess)) THEN
       DO ispin=1, SIZE(pao%matrix_P_guess)
          CALL cp_dbcsr_release(pao%matrix_P_guess(ispin))
       ENDDO
       DEALLOCATE(pao%matrix_P_guess)
    ENDIF

    CALL linesearch_finalize(pao%linesearch)
    !TODO: should finish printkey
    CALL timestop(handle)
  END SUBROUTINE pao_finalize


! *****************************************************************************
!> \brief Called after the SCF optimization, updates the PAO basis.
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param pao_is_done ...
! *****************************************************************************
  SUBROUTINE pao_update(qs_env,ls_scf_env,pao_is_done)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    LOGICAL, INTENT(OUT)                     :: pao_is_done

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: beta
    TYPE(cp_dbcsr_type)                      :: matrix_M
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_energy_type), POINTER            :: energy

    IF(.NOT. ls_scf_env%do_pao) THEN
       pao_is_done = .TRUE.
       RETURN
    ENDIF

    CALL timeset(routineN,handle)
    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env
    CALL get_qs_env(qs_env, energy=energy)
    pao_is_done = .FALSE.

    IF(pao%istep >= pao%max_pao) THEN
        CPWARN("PAO not converged!")
        pao_is_done = .TRUE.
    ELSE

       IF(pao%linesearch%starts) THEN
          ! store prim. matrix_P before linesearch starts to wander off
          CALL pao_store_P(qs_env, ls_scf_env)
          CALL pao_write_restart(pao, qs_env, energy%total)

          ! calc new gradient
          CALL pao_calc_outer_grad_orig(qs_env, ls_scf_env, matrix_M)
          CALL pao_calc_grad(pao, matrix_M)
          CALL cp_dbcsr_release(matrix_M)

          ! determine CG mixing factor
          CALL pao_calc_CG_beta(pao, beta)

          ! calculate new CG direction matrix_D
          CALL cp_dbcsr_add(pao%matrix_D, pao%matrix_G, beta, 1.0_dp)

          ! backup gradient for calc_CG_beta
          CALL cp_dbcsr_copy(pao%matrix_G_prev, pao%matrix_G)

          !backup X
          CALL cp_dbcsr_copy(pao%matrix_X_orig, pao%matrix_X)

          ! print info and convergence test
          CALL pao_test_convergence(pao, ls_scf_env, energy%total, "SCF", pao_is_done)
          IF(pao_is_done .AND. pao%iw>0) WRITE(pao%iw,*) "PAO| converged after ", pao%istep, " steps :-)"
       ENDIF

       IF(.NOT.pao_is_done) THEN
          CALL linesearch_step(pao%linesearch,&
                               energy=energy%total,&
                               slope=cp_dbcsr_frobenius_norm(pao%matrix_G))

          CALL cp_dbcsr_copy(pao%matrix_X, pao%matrix_X_orig) !restore X
          CALL cp_dbcsr_add(pao%matrix_X, pao%matrix_D, 1.0_dp, -pao%linesearch%step_size)
          CALL pao_update_AB(pao, ls_mstruct)
       ENDIF
    ENDIF

    CALL timestop(handle)
  END SUBROUTINE pao_update


! *****************************************************************************
!> \brief Stores density matrix as initial guess for next SCF optimization.
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_store_P(qs_env, ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    INTEGER                                  :: ispin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type)                      :: matrix_SB, matrix_tmp
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao

    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env

    IF(pao%iw>0) WRITE (pao%iw,*) "PAO| Storing primary density matrix for next PAO interation."

    CALL get_qs_env(qs_env, matrix_s=matrix_s)

    !lazy allocation
    IF(.NOT.ASSOCIATED(pao%matrix_P_guess)) THEN
       ALLOCATE(pao%matrix_P_guess(ls_scf_env%nspins))
       DO ispin=1, ls_scf_env%nspins
          CALL cp_dbcsr_init(pao%matrix_P_guess(ispin))
          CALL cp_dbcsr_create(pao%matrix_P_guess(ispin), template=matrix_s(1)%matrix, matrix_type="N")
       ENDDO
    ENDIF

    CALL cp_dbcsr_init(matrix_SB)
    CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B, &
                           0.0_dp, matrix_SB)


    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp, template=ls_mstruct%matrix_B)
    DO ispin=1, ls_scf_env%nspins
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_SB, ls_scf_env%matrix_p(ispin), &
                              0.0_dp, matrix_tmp)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_SB, &
                              0.0_dp, pao%matrix_P_guess(ispin))
    ENDDO

    CALL cp_dbcsr_release(matrix_SB)
    CALL cp_dbcsr_release(matrix_tmp)
  END SUBROUTINE pao_store_P


! *****************************************************************************
!> \brief Build rectangular matrix to converert between primary and PAO basis.
!> \param qs_env ...
!> \param matrix_Y ...
! *****************************************************************************
  SUBROUTINE pao_build_selector(qs_env, matrix_Y)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_type)                      :: matrix_Y

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_selector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, ikind, &
                                                nblocks, pao_basis_size
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes_aux, blk_sizes_pri
    REAL(dp), DIMENSION(:, :), POINTER       :: block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    NULLIFY(ks_env, matrix_s)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set)

    blk_sizes_pri => cp_dbcsr_col_block_sizes(matrix_s(1)%matrix)
    nblocks = SIZE(blk_sizes_pri)
    IF(nblocks /= SIZE(particle_set, 1)) CPABORT("nblocks wrong")

    ALLOCATE(blk_sizes_aux(nblocks))
    DO i=1, nblocks
      CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), pao_basis_size=pao_basis_size)
      IF(pao_basis_size < 1) CPABORT("PAO_BASIS_SIZE too small.")
      blk_sizes_aux(i) = pao_basis_size
    ENDDO

    CALL cp_dbcsr_init(matrix_Y)
    CALL cp_dbcsr_create(matrix_Y,&
                         template=matrix_s(1)%matrix,&
                         matrix_type="N",&
                         row_blk_size=blk_sizes_pri,&
                         col_blk_size=blk_sizes_aux,&
                         name="PAO matrix_Y")
    DEALLOCATE(blk_sizes_aux)

    CALL cp_dbcsr_reserve_diag_blocks(matrix_Y)

    CALL cp_dbcsr_iterator_start(iter, matrix_Y)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block)
       IF(SIZE(block,2)>SIZE(block,1)) CPABORT("Unexpected block sizes")
       block = 0.0_dp
       DO i=1, SIZE(block, 2)
          block(i,i) = 1.0_dp
       ENDDO
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE pao_build_selector


! *****************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M ...
! *****************************************************************************
   SUBROUTINE pao_calc_outer_grad_orig(qs_env,ls_scf_env,matrix_M)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    TYPE(cp_dbcsr_type)                      :: matrix_M

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_outer_grad_orig', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_dbcsr_type) :: matrix_HB, matrix_HBP, matrix_HP, matrix_M_pre, &
      matrix_PHP, matrix_SB, matrix_SBPHP, matrix_tmp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env

    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control)
    nspin = dft_control%nspins

    !TODO: handle openshell case properly
    IF(nspin /= 1) CPABORT("open shell not yet implemented")

    CALL cp_dbcsr_init(matrix_HB)
    CALL cp_dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_HB)

    CALL cp_dbcsr_init(matrix_HBP)
    CALL cp_dbcsr_create(matrix_HBP, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_HB, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_HBP)

    CALL cp_dbcsr_init(matrix_HP)
    CALL cp_dbcsr_create(matrix_HP, template=ls_scf_env%matrix_p(1), matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_ks(1), ls_scf_env%matrix_p(1),&
                          0.0_dp, matrix_HP)

    CALL cp_dbcsr_init(matrix_PHP)
    CALL cp_dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_p(1), matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), matrix_HP,&
                           0.0_dp, matrix_PHP)

    CALL cp_dbcsr_init(matrix_SB)
    CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_SB)

    CALL cp_dbcsr_init(matrix_SBPHP)
    CALL cp_dbcsr_create(matrix_SBPHP, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_SB, matrix_PHP,&
                           0.0_dp, matrix_SBPHP)

    IF(SIZE(ls_scf_env%matrix_p)==1) CALL cp_dbcsr_scale(matrix_SBPHP, 0.5_dp)

    CALL cp_dbcsr_init(matrix_M_pre)
    CALL cp_dbcsr_copy(matrix_M_pre, matrix_HBP)
    CALL cp_dbcsr_add(matrix_M_pre, matrix_SBPHP, 2.0_dp, -2.0_dp)

    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp, template=matrix_M_pre, matrix_type="N")
    CALL cp_dbcsr_multiply("T", "N", 1.0_dp, pao%matrix_N, matrix_M_pre,&
                           0.0_dp, matrix_tmp)

    CALL cp_dbcsr_init(matrix_M)
    CALL cp_dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, pao%matrix_Y,&
                           0.0_dp, matrix_M)

    CALL cp_dbcsr_release(matrix_HB)
    CALL cp_dbcsr_release(matrix_HBP)
    CALL cp_dbcsr_release(matrix_HP)
    CALL cp_dbcsr_release(matrix_PHP)
    CALL cp_dbcsr_release(matrix_SB)
    CALL cp_dbcsr_release(matrix_SBPHP)
    CALL cp_dbcsr_release(matrix_M_pre)
    CALL cp_dbcsr_release(matrix_tmp)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_outer_grad_orig

END MODULE pao_main

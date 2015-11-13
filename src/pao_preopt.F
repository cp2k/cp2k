!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief PAO pre-optimization
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_preopt
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_binary_read, cp_dbcsr_checksum, &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_copy_into_existing, &
       cp_dbcsr_create, cp_dbcsr_desymmetrize, cp_dbcsr_distribution, &
       cp_dbcsr_frobenius_norm, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type
  USE cp_log_handling,                 ONLY: cp_to_string
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dm_ls_scf_qs,                    ONLY: matrix_ls_create
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE kinds,                           ONLY: dp
  USE linesearch,                      ONLY: linesearch_reset,&
                                             linesearch_step
  USE message_passing,                 ONLY: mp_max
  USE pao_input,                       ONLY: pao_exp_param
  USE pao_io,                          ONLY: pao_write_restart
  USE pao_param,                       ONLY: pao_calc_grad,&
                                             pao_update_AB
  USE pao_preopt_mcweeny,              ONLY: pao_calc_outer_grad_mcweeny,&
                                             pao_preopt_dm_mcweeny,&
                                             pao_preopt_finalize_mcweeny
  USE pao_preopt_trs4,                 ONLY: pao_calc_outer_grad_lnv,&
                                             pao_preopt_dm_trs4,&
                                             pao_preopt_finalize_trs4
  USE pao_types,                       ONLY: pao_env_type
  USE pao_utils,                       ONLY: pao_calc_CG_beta,&
                                             pao_print_oneliner,&
                                             pao_print_step_delimiter,&
                                             pao_test_convegence
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_did_change
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_preopt'

  PUBLIC :: pao_preopt_run

CONTAINS

! *****************************************************************************
!> \brief Main routine which performs the pre-optimization
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_preopt_run(pao, qs_env, ls_scf_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_preopt_run', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: mode
    INTEGER                                  :: handle
    LOGICAL                                  :: preopt_is_done
    REAL(KIND=dp)                            :: beta, energy
    TYPE(cp_dbcsr_type)                      :: matrix_M

    IF(LEN_TRIM(pao%preopt_dm_file)==0) RETURN ! preopt disabled

    CALL timeset(routineN,handle)

    CALL pao_preopt_init(pao, qs_env, ls_scf_env)

    preopt_is_done = .FALSE.
    mode = "TRS4"

    DO WHILE(.NOT.preopt_is_done)
       IF(pao%istep+1 >= pao%max_pao) EXIT ! leave the last step for SCF-PAO
       pao%istep = pao%istep + 1
       CALL pao_print_step_delimiter(pao)

       ! calculate new density matrix
       IF(TRIM(mode) == "TRS4") THEN
          CALL pao_preopt_dm_trs4(qs_env, ls_scf_env)
       ELSE
          CALL pao_preopt_dm_mcweeny(pao, qs_env, ls_scf_env)
       ENDIF

       ! calc energy and check trace_PS
       CALL pao_calc_trace_PH(ls_scf_env, energy)
       CALL pao_check_trace_PS(ls_scf_env)

       IF(pao%linesearch%starts) THEN
          ! calculate new gradient
          IF(TRIM(mode) == "TRS4") THEN
             CALL pao_calc_outer_grad_lnv(qs_env, ls_scf_env, matrix_M)
          ELSE
             CALL pao_calc_outer_grad_mcweeny(qs_env, ls_scf_env, matrix_M)
          ENDIF
          CALL pao_calc_grad(pao, matrix_M)
          CALL cp_dbcsr_release(matrix_M)

          ! print info
          CALL pao_print_oneliner(pao, ls_scf_env, energy, mode)

          ! determine CG mixing factor
          CALL pao_calc_CG_beta(pao, beta)

          ! calculate new CG direction matrix_D
          CALL cp_dbcsr_add(pao%matrix_D, pao%matrix_G, beta, 1.0_dp)

          ! backup gradient for calc_CG_beta
          CALL cp_dbcsr_copy(pao%matrix_G_prev, pao%matrix_G)

          !backup X
          CALL cp_dbcsr_copy(pao%matrix_X_orig, pao%matrix_X)

          IF(TRIM(mode) == "McWeeny") THEN
             CALL pao_check_grad_mcweeny(pao, qs_env, ls_scf_env)
             CALL pao_test_convegence(pao, ls_scf_env, pao%eps_preopt, energy, preopt_is_done)
          ENDIF
       ENDIF

       IF(.NOT.preopt_is_done) THEN
          CALL linesearch_step(pao%linesearch,&
                               energy=energy,&
                               slope=cp_dbcsr_frobenius_norm(pao%matrix_G))
          CALL cp_dbcsr_copy(pao%matrix_X, pao%matrix_X_orig) !restore X
          CALL cp_dbcsr_add(pao%matrix_X, pao%matrix_D, 1.0_dp, -pao%linesearch%step_size)
          CALL pao_update_AB(pao, ls_scf_env%ls_mstruct)
       ENDIF

       ! switch from TRS4 to McWeeny
       IF(TRIM(mode) == "TRS4" .AND. pao%linesearch%starts .AND. pao%istep>pao%max_trs4_preopt) THEN
          CALL pao_preopt_finalize_trs4(ls_scf_env)
          CALL linesearch_reset(pao%linesearch)
          CALL cp_dbcsr_set(pao%matrix_D, 0.0_dp) ! reset CG
          CALL cp_dbcsr_set(pao%matrix_G_prev, 0.0_dp)
          mode = "McWeeny"
       END IF
    ENDDO

    CALL pao_write_restart(pao, qs_env)

    CALL pao_preopt_finalize(pao, ls_scf_env, mode)

    CALL timestop(handle)
  END SUBROUTINE pao_preopt_run


! *****************************************************************************
!> \brief Initialize the pre-optimization
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_preopt_init(pao, qs_env, ls_scf_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_preopt_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    REAL(KIND=dp)                            :: cs_pos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_dbcsr_type)                      :: matrix_tmp
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    nspin = ls_scf_env%nspins

    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    matrix_ks=matrix_ks,&
                    rho=rho,&
                    energy=energy)

    CALL qs_rho_get(rho, rho_ao=rho_ao)

    IF(nspin/=1) CPABORT("open shell not yet implemented")

    DO ispin=1, nspin
       CALL cp_dbcsr_init(matrix_tmp) ! read into tmp-matrix to preserve sparsity of rho_ao
       CALL cp_dbcsr_binary_read(pao%preopt_dm_file, matrix_new=matrix_tmp,& 
                                 distribution=cp_dbcsr_distribution(matrix_s(1)%matrix))
       cs_pos = cp_dbcsr_checksum(matrix_tmp, pos=.TRUE.)
       IF (pao%iw>0) WRITE(pao%iw,'(T2,A,E20.8)') "Read restart DM "//TRIM(pao%preopt_dm_file)//&
                                                  " with checksum: ",cs_pos
       CALL cp_dbcsr_copy_into_existing(rho_ao(ispin)%matrix, matrix_tmp)
       CALL cp_dbcsr_release(matrix_tmp)
    ENDDO

    ! calculate corresponding ks matrix
    CALL qs_rho_update_rho(rho,qs_env=qs_env)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.)
    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE.,&
                             just_energy=.FALSE., print_active=.TRUE.)
    IF (pao%iw>0) WRITE(pao%iw,*) "PAO| Quickstep energy from restart density:",energy%total

    ! allocate LS matrices
    CALL matrix_ls_create(matrix_ls=ls_scf_env%matrix_s, matrix_qs=matrix_s(1)%matrix,&
                          ls_mstruct=ls_scf_env%ls_mstruct)

    ALLOCATE(ls_scf_env%matrix_p(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(ls_scf_env%matrix_p(ispin))
       CALL cp_dbcsr_create(ls_scf_env%matrix_p(ispin),template=ls_scf_env%matrix_s,&
                            matrix_type="N")
    ENDDO

    ALLOCATE(ls_scf_env%matrix_ks(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(ls_scf_env%matrix_ks(ispin))
       CALL cp_dbcsr_create(ls_scf_env%matrix_ks(ispin),template=ls_scf_env%matrix_s,&
                            matrix_type="N")
    ENDDO

    ! prevent warnings from cp_dbcsr_release()
    CALL cp_dbcsr_init(ls_scf_env%matrix_s_inv)
    CALL cp_dbcsr_init(ls_scf_env%matrix_s_sqrt)
    CALL cp_dbcsr_init(ls_scf_env%matrix_s_sqrt_inv)

    CALL timestop(handle)
  END SUBROUTINE pao_preopt_init


! *****************************************************************************
!> \brief Finalize the pre-optimization
!> \param pao ...
!> \param ls_scf_env ...
!> \param mode ...
! *****************************************************************************
  SUBROUTINE pao_preopt_finalize(pao, ls_scf_env, mode)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    CHARACTER(LEN=10)                        :: mode

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_preopt_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin

    CALL timeset(routineN,handle)

    IF(TRIM(mode) == "TRS4")&
       CALL pao_preopt_finalize_trs4(ls_scf_env)

    IF(TRIM(mode) == "McWeeny")&
       CALL pao_preopt_finalize_mcweeny(pao)

    CALL cp_dbcsr_release(ls_scf_env%matrix_s)

    DO ispin=1, ls_scf_env%nspins
       CALL cp_dbcsr_release(ls_scf_env%matrix_p(ispin))
       CALL cp_dbcsr_release(ls_scf_env%matrix_ks(ispin))
    ENDDO

    DEALLOCATE(ls_scf_env%matrix_ks, ls_scf_env%matrix_p)

    CALL linesearch_reset(pao%linesearch)
    CALL cp_dbcsr_set(pao%matrix_D, 0.0_dp)
    CALL cp_dbcsr_set(pao%matrix_G_prev, 0.0_dp)

    CALL timestop(handle)
  END SUBROUTINE pao_preopt_finalize


! *****************************************************************************
!> \brief Calculate the pre-optimization energy
!> \param ls_scf_env ...
!> \param energy ...
! *****************************************************************************
  SUBROUTINE pao_calc_trace_PH(ls_scf_env, energy)
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    REAL(KIND=dp), INTENT(OUT)               :: energy

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_trace_PH', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: trace_PH

    CALL timeset(routineN,handle)
    energy = 0.0_dp
    DO ispin=1, ls_scf_env%nspins
       CALL cp_dbcsr_trace(ls_scf_env%matrix_p(ispin), ls_scf_env%matrix_ks(ispin), trace_PH)
       energy = energy + trace_PH
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE pao_calc_trace_PH


! *****************************************************************************
!> \brief Ensure that the number of electrons is correct.
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_check_trace_PS(ls_scf_env)
    TYPE(ls_scf_env_type)                    :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_check_trace_PS', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: tmp, trace_PS
    TYPE(cp_dbcsr_type)                      :: matrix_S_desym

    CALL timeset(routineN,handle)
    CALL cp_dbcsr_init(matrix_S_desym)
    CALL cp_dbcsr_create(matrix_S_desym, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_desymmetrize(ls_scf_env%matrix_s, matrix_S_desym)

    trace_PS = 0.0_dp
    DO ispin=1, ls_scf_env%nspins
       CALL cp_dbcsr_trace(ls_scf_env%matrix_p(ispin), matrix_S_desym, tmp)
       trace_PS = trace_PS + tmp
    ENDDO

    CALL cp_dbcsr_release(matrix_S_desym)

    IF(ABS(ls_scf_env%nelectron_total - trace_PS) > 0.5)&
         CPABORT("Number of electrons wrong. Trace(PS) ="//cp_to_string(trace_PS))

    CALL timestop(handle)
  END SUBROUTINE pao_check_trace_PS


! *****************************************************************************
!> \brief Debugging routine for checking the analytic gradient.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_check_grad_mcweeny(pao, qs_env, ls_scf_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_check_grad_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, j, natoms
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes_col, blk_sizes_row
    LOGICAL                                  :: found
    REAL(dp)                                 :: delta_max, energy1, energy2, &
                                                eps, Gij_num, symm
    REAL(dp), DIMENSION(:, :), POINTER       :: block_G, block_X, block_X0
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct

    IF(pao%check_grad_mcweeny_tol<0.0_dp) RETURN ! no checking

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct

    SELECT CASE(pao%parameterization)
      CASE(pao_exp_param)
        symm = -1.0_dp ! anti-symmetric
      CASE DEFAULT
        symm = 0.0_dp  ! no symmetry
    END SELECT

    CALL get_qs_env(qs_env, para_env=para_env, natom=natoms)

    eps = 1.0e-5_dp
    delta_max = 0.0_dp

    blk_sizes_col => cp_dbcsr_col_block_sizes(pao%matrix_X)
    blk_sizes_row => cp_dbcsr_row_block_sizes(pao%matrix_X)

    ! can not use an iterator here, because other DBCSR routines are called within loop.
    DO iatom=1, natoms
       IF(pao%iw>0) WRITE(pao%iw,*) 'PAO| checking McWeeny gradient of atom ',iatom
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)

       IF(ASSOCIATED(block_X)) THEN !only one node actually has the block
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
          CPASSERT(ASSOCIATED(block_G))
          ALLOCATE(block_X0(blk_sizes_row(iatom),blk_sizes_col(iatom)))
          block_X0 = block_X
       ENDIF

       DO i=1, blk_sizes_row(iatom)
          DO j=1, blk_sizes_col(iatom)
             IF(ASSOCIATED(block_X)) THEN
                block_X = block_X0
                block_X(i,j) = block_X(i,j) + eps
                IF(symm/=0.0) block_X(j,i) = block_X(j,i) + symm*eps
             ENDIF

             CALL pao_update_AB(pao, ls_mstruct)
             CALL pao_preopt_dm_mcweeny(pao, qs_env, ls_scf_env)
             CALL pao_calc_trace_PH(ls_scf_env, energy1)

             IF(ASSOCIATED(block_X)) THEN
                block_X = block_X0
                block_X(i,j) = block_X(i,j) - eps
                IF(symm/=0.0) block_X(j,i) = block_X(j,i) - symm*eps
             ENDIF

             CALL pao_update_AB(pao, ls_mstruct)
             CALL pao_preopt_dm_mcweeny(pao, qs_env, ls_scf_env)
             CALL pao_calc_trace_PH(ls_scf_env, energy2)

             IF(ASSOCIATED(block_X)) THEN
                Gij_num = (energy1 - energy2) / (2.0_dp*eps)
                delta_max = MAX(delta_max, ABS(Gij_num - block_G(i,j)))
             ENDIF
          ENDDO
       ENDDO

       IF(ASSOCIATED(block_X)) THEN
          block_X = block_X0
          DEALLOCATE(block_X0)
       ENDIF
    END DO

    CALL mp_max(delta_max, para_env%group)
    IF(pao%iw>0) WRITE(pao%iw,*) 'PAO| checked McWeeny gradient, max delta:', delta_max
    IF(delta_max > pao%check_grad_mcweeny_tol)  CALL cp_abort(__LOCATION__,&
       "Analytic and numeric gradients of McWeeny differ too much:"//cp_to_string(delta_max))

    CALL timestop(handle)
  END SUBROUTINE pao_check_grad_mcweeny

END MODULE pao_preopt

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Main module for the PAO method
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_main
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_init, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_reserve_diag_blocks, cp_dbcsr_set, cp_dbcsr_type
  USE cp_external_control,             ONLY: external_control
  USE dm_ls_scf_methods,               ONLY: ls_scf_init_matrix_S
  USE dm_ls_scf_qs,                    ONLY: ls_scf_qs_atomic_guess
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE linesearch,                      ONLY: linesearch_finalize,&
                                             linesearch_init,&
                                             linesearch_reset,&
                                             linesearch_step
  USE machine,                         ONLY: m_walltime
  USE pao_input,                       ONLY: parse_pao_section
  USE pao_io,                          ONLY: pao_read_restart,&
                                             pao_write_restart
  USE pao_methods,                     ONLY: &
       pao_build_orthogonalizer, pao_build_selector, pao_calc_CG_beta, &
       pao_calc_outer_grad_lnv, pao_calc_trace_ph, pao_check_grad_full, &
       pao_check_trace_ps, pao_dm_trs4, pao_print_atom_info, &
       pao_read_preopt_dm, pao_test_convergence
  USE pao_param,                       ONLY: pao_calc_grad,&
                                             pao_param_finalize,&
                                             pao_param_init,&
                                             pao_update_AB
  USE pao_types,                       ONLY: pao_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_main'

  PUBLIC :: pao_init, pao_update, pao_finalize

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

    ! get initial P and KS matrices in primary basis
    IF(LEN_TRIM(pao%preopt_dm_file)>0) THEN
       CALL pao_read_preopt_dm(pao, qs_env)
    ELSE
       CALL ls_scf_qs_atomic_guess(qs_env, ls_scf_env%energy_init)
       IF(pao%iw>0) WRITE(pao%iw,'(A,F20.9)') " PAO| Energy from initial atomic guess:", ls_scf_env%energy_init
    ENDIF

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

    CALL pao_print_atom_info(pao)

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

    CALL timestop(handle)
  END SUBROUTINE pao_init


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

    INTEGER                                  :: handle, icycle
    LOGICAL                                  :: cycle_converged, should_stop
    REAL(KIND=dp)                            :: beta, energy
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type)                      :: matrix_M, matrix_X_start
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao

    IF(.NOT. ls_scf_env%do_pao) THEN
       pao_is_done = .TRUE.
       RETURN
    ENDIF

    CALL timeset(routineN,handle)
    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env

    ! reset CG
    CALL linesearch_reset(pao%linesearch)
    CALL cp_dbcsr_set(pao%matrix_D, 0.0_dp)
    CALL cp_dbcsr_set(pao%matrix_G_prev, 0.0_dp)

    IF (pao%mixing /= 1.0_dp) THEN
       CALL cp_dbcsr_init(matrix_X_start)
       CALL cp_dbcsr_copy(matrix_X_start, pao%matrix_X)
    ENDIF

    cycle_converged = .FALSE.
    icycle = 0

    DO WHILE(.TRUE.)
       pao%istep = pao%istep + 1

       IF(pao%iw>0) WRITE(pao%iw,"(A,I9,A)") " PAO| ======================= Iteration: ",&
                 pao%istep," ============================="

       ! calculate new density matrix
       CALL pao_dm_trs4(qs_env, ls_scf_env)

       ! calc energy and check trace_PS
       CALL pao_calc_trace_PH(ls_scf_env, energy)
       CALL pao_check_trace_PS(ls_scf_env)

       IF(pao%linesearch%starts) THEN
          icycle = icycle + 1
          ! calc new gradient
          CALL pao_calc_outer_grad_lnv(qs_env, ls_scf_env, matrix_M)
          CALL pao_calc_grad(pao, matrix_M)
          CALL cp_dbcsr_release(matrix_M)
          CALL pao_check_grad_full(pao, qs_env, ls_scf_env)

          ! determine CG mixing factor
          CALL pao_calc_CG_beta(pao, beta)

          ! calculate new CG direction matrix_D
          CALL cp_dbcsr_add(pao%matrix_D, pao%matrix_G, beta, 1.0_dp)

          ! backup gradient for calc_CG_beta
          CALL cp_dbcsr_copy(pao%matrix_G_prev, pao%matrix_G)

          !backup X
          CALL cp_dbcsr_copy(pao%matrix_X_orig, pao%matrix_X)

          ! print info and convergence test
          CALL pao_test_convergence(pao, ls_scf_env, energy, icycle, cycle_converged)
          IF(cycle_converged) THEN
             pao_is_done = icycle < 3
             IF(pao_is_done .AND. pao%iw>0) WRITE(pao%iw,*) "PAO| converged after ", pao%istep, " steps :-)"
             EXIT
          ENDIF
         
          ! if we have reached the maximum number of cycles exit in order 
          ! to restart with a fresh hamiltonian
          IF (icycle >= pao%max_cycles) THEN
             IF (pao%iw>0) WRITE(pao%iw,*) "PAO| CG not yet converged after ", icycle , " cylces."
             EXIT
          ENDIF
       ENDIF

       ! check for early abort without convergence?
       CALL external_control(should_stop,"PAO",start_time=qs_env%start_time,target_time=qs_env%target_time)
       IF (should_stop .OR. pao%istep >= pao%max_pao) THEN
          CPWARN("PAO not converged!")
          pao_is_done = .TRUE.
          EXIT
       ENDIF

       ! perform line-search step
       CALL linesearch_step(pao%linesearch,&
                            energy=energy,&
                            slope=cp_dbcsr_frobenius_norm(pao%matrix_G))
       CALL cp_dbcsr_copy(pao%matrix_X, pao%matrix_X_orig) !restore X
       CALL cp_dbcsr_add(pao%matrix_X, pao%matrix_D, 1.0_dp, -pao%linesearch%step_size)
       CALL pao_update_AB(pao, ls_mstruct)
    ENDDO
 
    IF (pao%mixing /= 1.0_dp) THEN
       CALL cp_dbcsr_add(pao%matrix_X,matrix_X_start,pao%mixing,1.0_dp-pao%mixing)

       ! update matrices to match the new X
       CALL pao_update_AB(pao, ls_mstruct)
       CALL cp_dbcsr_release(ls_scf_env%matrix_s_inv)
       CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt)
       CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt_inv)
       CALL get_qs_env(qs_env, matrix_s=matrix_s)
       CALL ls_scf_init_matrix_s(matrix_s(1)%matrix, ls_scf_env)
       CALL cp_dbcsr_release(matrix_X_start)
    ENDIF

    CALL pao_write_restart(pao, qs_env, energy)

    CALL timestop(handle)
  END SUBROUTINE pao_update


! *****************************************************************************
!> \brief Finalize the PAO environment
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_finalize(ls_scf_env)
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
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

    CALL linesearch_finalize(pao%linesearch)
    !TODO: should finish printkey
    CALL timestop(handle)
  END SUBROUTINE pao_finalize

END MODULE pao_main

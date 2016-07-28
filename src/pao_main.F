!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main module for the PAO method
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_main
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_init, cp_dbcsr_p_type, &
        cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, cp_dbcsr_set, cp_dbcsr_type
   USE cp_external_control,             ONLY: external_control
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
        pao_add_forces, pao_build_core_hamiltonian, pao_build_diag_distribution, &
        pao_build_matrix_X, pao_build_orthogonalizer, pao_build_selector, pao_calc_energy, &
        pao_calc_outer_grad_lnv, pao_check_grad, pao_check_trace_ps, pao_guess_initial_P, &
        pao_init_gaussian_radii, pao_print_atom_info, pao_store_P, pao_test_convergence
   USE pao_ml,                          ONLY: pao_ml_init,&
                                              pao_ml_predict
   USE pao_optimizer,                   ONLY: pao_opt_finalize,&
                                              pao_opt_init,&
                                              pao_opt_new_dir
   USE pao_param,                       ONLY: pao_calc_U,&
                                              pao_param_finalize,&
                                              pao_param_init,&
                                              pao_param_initial_guess,&
                                              pao_update_AB
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_main'

   PUBLIC :: pao_init, pao_update, pao_post_scf, pao_optimization_start, pao_optimization_end

CONTAINS

! **************************************************************************************************
!> \brief Initialize the PAO environment
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_init(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_init'

      INTEGER                                            :: handle
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(section_vals_type), POINTER                   :: input

      IF (.NOT. ls_scf_env%do_pao) RETURN

      CALL timeset(routineN, handle)
      pao => ls_scf_env%pao_env
      CALL get_qs_env(qs_env=qs_env, input=input, matrix_s=matrix_s)

      ! parse input
      CALL parse_pao_section(pao, input)

      CALL pao_init_gaussian_radii(pao, qs_env)

      ! train machine learning
      CALL pao_ml_init(pao, qs_env)

      CALL timestop(handle)
   END SUBROUTINE pao_init

! **************************************************************************************************
!> \brief Start a PAO optimization run.
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_optimization_start(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_optimization_start'

      INTEGER                                            :: handle
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(section_vals_type), POINTER                   :: input, section

      IF (.NOT. ls_scf_env%do_pao) RETURN

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, &
                      input=input)

      pao => ls_scf_env%pao_env
      ls_mstruct => ls_scf_env%ls_mstruct

      ! reset state
      pao%step_start_time = m_walltime()
      pao%istep = 0
      pao%matrix_P_ready = .FALSE.

      ! ready stuff that does not depend on atom positions
      IF (.NOT. pao%constants_ready) THEN
         CALL pao_build_diag_distribution(pao, qs_env)
         CALL pao_build_orthogonalizer(pao, qs_env)
         CALL pao_build_selector(pao, qs_env)
         CALL pao_build_core_hamiltonian(pao, qs_env)
         pao%constants_ready = .TRUE.
      ENDIF

      CALL pao_param_init(pao, qs_env)

      ! ready PAO parameter matrix_X
      IF (.NOT. pao%matrix_X_ready) THEN
         CALL pao_build_matrix_X(pao, qs_env)
         CALL pao_print_atom_info(pao)
         IF (LEN_TRIM(pao%restart_file) > 0) THEN
            CALL pao_read_restart(pao, qs_env)
         ELSE IF (SIZE(pao%ml_training_set) > 0) THEN
            CALL pao_ml_predict(pao, qs_env)
         ELSE
            CALL pao_param_initial_guess(pao, qs_env)
         ENDIF
         pao%matrix_X_ready = .TRUE.
      ELSE IF (SIZE(pao%ml_training_set) > 0) THEN
         CALL pao_ml_predict(pao, qs_env)
      ELSE
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| reusing matrix_X from previous optimization"
      ENDIF

      ! init line-search
      section => section_vals_get_subs_vals(input, "DFT%LS_SCF%PAO%LINE_SEARCH")
      CALL linesearch_init(pao%linesearch, section, "PAO|")

      ! create some more matrices
      CALL cp_dbcsr_init(pao%matrix_X_orig)
      CALL cp_dbcsr_init(pao%matrix_G)
      CALL cp_dbcsr_copy(pao%matrix_G, pao%matrix_X)
      CALL cp_dbcsr_set(pao%matrix_G, 0.0_dp)

      CALL cp_dbcsr_init(pao%matrix_U)
      CALL cp_dbcsr_create(pao%matrix_U, &
                           name="PAO matrix_U", &
                           matrix_type="N", &
                           dist=pao%diag_distribution, &
                           template=matrix_s(1)%matrix)
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_U)

      CALL cp_dbcsr_init(ls_mstruct%matrix_A)
      CALL cp_dbcsr_create(ls_mstruct%matrix_A, template=pao%matrix_Y)
      CALL cp_dbcsr_init(ls_mstruct%matrix_B)
      CALL cp_dbcsr_create(ls_mstruct%matrix_B, template=pao%matrix_Y)

      ! fill PAO transformation matrices
      CALL pao_update_AB(pao, qs_env, ls_mstruct)

      CALL timestop(handle)
   END SUBROUTINE pao_optimization_start

! **************************************************************************************************
!> \brief Called after the SCF optimization, updates the PAO basis.
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param pao_is_done ...
! **************************************************************************************************
   SUBROUTINE pao_update(qs_env, ls_scf_env, pao_is_done)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(OUT)                               :: pao_is_done

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_update'

      INTEGER                                            :: handle, icycle
      LOGICAL                                            :: cycle_converged, do_mixing, should_stop
      REAL(KIND=dp)                                      :: energy, penalty
      TYPE(cp_dbcsr_type)                                :: matrix_M, matrix_X_mixing
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao

      IF (.NOT. ls_scf_env%do_pao) THEN
         pao_is_done = .TRUE.
         RETURN
      ENDIF

      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env

      IF (.NOT. pao%matrix_P_ready) THEN
         CALL pao_guess_initial_P(pao, qs_env, ls_scf_env)
         pao%matrix_P_ready = .TRUE.
      ENDIF

      IF (pao%max_pao == 0) THEN
         pao_is_done = .TRUE.
         RETURN
      ENDIF

      IF (pao%need_initial_scf) THEN
         pao_is_done = .FALSE.
         pao%need_initial_scf = .FALSE.
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Performing initial SCF optimization."
         RETURN
      ENDIF

      CALL timeset(routineN, handle)

      ! perform mixing once we are well into the optimization
      do_mixing = pao%mixing /= 1.0_dp .AND. pao%istep > 1
      IF (do_mixing) THEN
         CALL cp_dbcsr_init(matrix_X_mixing)
         CALL cp_dbcsr_copy(matrix_X_mixing, pao%matrix_X)
      ENDIF

      cycle_converged = .FALSE.
      icycle = 0
      CALL linesearch_reset(pao%linesearch)
      CALL pao_opt_init(pao)

      DO WHILE (.TRUE.)
         pao%istep = pao%istep+1

         IF (pao%iw > 0) WRITE (pao%iw, "(A,I9,A)") " PAO| ======================= Iteration: ", &
            pao%istep, " ============================="

         ! calc energy and check trace_PS
         CALL pao_calc_energy(pao, qs_env, ls_scf_env, energy)
         CALL pao_check_trace_PS(ls_scf_env)

         IF (pao%linesearch%starts) THEN
            icycle = icycle+1
            ! calc new gradient including penalty terms
            CALL pao_calc_outer_grad_lnv(qs_env, ls_scf_env, matrix_M)
            CALL pao_calc_U(pao, qs_env, matrix_M=matrix_M, matrix_G=pao%matrix_G, penalty=penalty)
            CALL cp_dbcsr_release(matrix_M)
            CALL pao_check_grad(pao, qs_env, ls_scf_env)

            ! calculate new direction for line-search
            CALL pao_opt_new_dir(pao, icycle)

            !backup X
            CALL cp_dbcsr_copy(pao%matrix_X_orig, pao%matrix_X)

            ! print info and convergence test
            CALL pao_test_convergence(pao, ls_scf_env, energy, cycle_converged)
            IF (cycle_converged) THEN
               pao_is_done = icycle < 3
               IF (pao_is_done .AND. pao%iw > 0) WRITE (pao%iw, *) "PAO| converged after ", pao%istep, " steps :-)"
               EXIT
            ENDIF

            ! if we have reached the maximum number of cycles exit in order
            ! to restart with a fresh hamiltonian
            IF (icycle >= pao%max_cycles) THEN
               IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| CG not yet converged after ", icycle, " cylces."
               pao_is_done = .FALSE.
               EXIT
            ENDIF

            IF (MOD(icycle, pao%write_cycles) == 0) &
               CALL pao_write_restart(pao, qs_env, energy) ! write an intermediate restart file
         ENDIF

         ! check for early abort without convergence?
         CALL external_control(should_stop, "PAO", start_time=qs_env%start_time, target_time=qs_env%target_time)
         IF (should_stop .OR. pao%istep >= pao%max_pao) THEN
            CPWARN("PAO not converged!")
            pao_is_done = .TRUE.
            EXIT
         ENDIF

         ! perform line-search step
         CALL linesearch_step(pao%linesearch, energy=energy, slope=pao%norm_G**2)

         IF (pao%linesearch%step_size < 1e-10_dp) CPABORT("PAO gradient is wrong.")

         CALL cp_dbcsr_copy(pao%matrix_X, pao%matrix_X_orig) !restore X
         CALL cp_dbcsr_add(pao%matrix_X, pao%matrix_D, 1.0_dp, pao%linesearch%step_size)
      ENDDO

      ! perform mixing of matrix_X
      IF (do_mixing) THEN
         CALL cp_dbcsr_add(pao%matrix_X, matrix_X_mixing, pao%mixing, 1.0_dp-pao%mixing)
         CALL cp_dbcsr_release(matrix_X_mixing)
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Recalculating energy after mixing."
         CALL pao_calc_energy(pao, qs_env, ls_scf_env, energy)
      ENDIF

      CALL pao_write_restart(pao, qs_env, energy)
      CALL pao_opt_finalize(pao)

      CALL timestop(handle)
   END SUBROUTINE pao_update

! **************************************************************************************************
!> \brief Calculate PAO forces and store density matrix for future ASPC extrapolations
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param pao_is_done ...
! **************************************************************************************************
   SUBROUTINE pao_post_scf(qs_env, ls_scf_env, pao_is_done)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: pao_is_done

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_post_scf'

      INTEGER                                            :: handle

      IF (.NOT. ls_scf_env%do_pao) RETURN
      IF (.NOT. pao_is_done) RETURN

      CALL timeset(routineN, handle)

      CALL pao_store_P(qs_env, ls_scf_env)
      IF (ls_scf_env%calculate_forces) CALL pao_add_forces(qs_env, ls_scf_env)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Finish a PAO optimization run.
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_optimization_end(ls_scf_env)
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_optimization_end'

      INTEGER                                            :: handle
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao

      IF (.NOT. ls_scf_env%do_pao) RETURN

      pao => ls_scf_env%pao_env
      ls_mstruct => ls_scf_env%ls_mstruct

      CALL timeset(routineN, handle)

      CALL pao_param_finalize(pao)

      ! We keep pao%matrix_X for next scf-run, e.g. during MD or GEO-OPT
      CALL cp_dbcsr_release(pao%matrix_X_orig)
      CALL cp_dbcsr_release(pao%matrix_G)
      CALL cp_dbcsr_release(pao%matrix_U)
      CALL cp_dbcsr_release(ls_mstruct%matrix_A)
      CALL cp_dbcsr_release(ls_mstruct%matrix_B)

      CALL linesearch_finalize(pao%linesearch)

      CALL timestop(handle)
   END SUBROUTINE pao_optimization_end

END MODULE pao_main

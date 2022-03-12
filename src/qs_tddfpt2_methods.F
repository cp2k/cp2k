!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_methods
   USE admm_types,                      ONLY: admm_type
   USE bibliography,                    ONLY: Grimme2013,&
                                              Grimme2016,&
                                              Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE exstates_types,                  ONLY: excited_energy_type
   USE header,                          ONLY: tddfpt_header
   USE input_constants,                 ONLY: tddfpt_dipole_velocity,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE lri_environment_methods,         ONLY: lri_print_stat
   USE lri_environment_types,           ONLY: lri_density_release,&
                                              lri_env_release
   USE machine,                         ONLY: m_flush
   USE mulliken,                        ONLY: compute_charges
   USE physcon,                         ONLY: evolt
   USE qs_2nd_kernel_ao,                ONLY: build_dm_response
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kernel_methods,               ONLY: create_kernel_env
   USE qs_kernel_types,                 ONLY: full_kernel_env_type,&
                                              kernel_env_type,&
                                              release_kernel_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_eigensolver,          ONLY: tddfpt_davidson_solver,&
                                              tddfpt_orthogonalize_psi1_psi0,&
                                              tddfpt_orthonormalize_psi1_psi1
   USE qs_tddfpt2_forces,               ONLY: tddfpt_forces,&
                                              tddfpt_resvec1,&
                                              tddfpt_resvec1_admm,&
                                              tddfpt_resvec2,&
                                              tddfpt_resvec2_xtb,&
                                              tddfpt_resvec3
   USE qs_tddfpt2_lri_utils,            ONLY: tddfpt2_lri_init
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_excitation_analysis,&
                                              tddfpt_print_nto_analysis,&
                                              tddfpt_print_summary
   USE qs_tddfpt2_restart,              ONLY: tddfpt_read_restart,&
                                              tddfpt_write_newtonx_output,&
                                              tddfpt_write_restart
   USE qs_tddfpt2_stda_types,           ONLY: allocate_stda_env,&
                                              deallocate_stda_env,&
                                              stda_env_type,&
                                              stda_init_param
   USE qs_tddfpt2_stda_utils,           ONLY: stda_init_matrices
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_sub_env_init,&
                                              tddfpt_sub_env_release,&
                                              tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: stda_create_work_matrices,&
                                              tddfpt_create_work_matrices,&
                                              tddfpt_ground_state_mos,&
                                              tddfpt_release_work_matrices,&
                                              tddfpt_work_matrices
   USE qs_tddfpt2_utils,                ONLY: tddfpt_guess_vectors,&
                                              tddfpt_init_mos,&
                                              tddfpt_oecorr,&
                                              tddfpt_release_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
   USE xc_write_output,                 ONLY: xc_write
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_methods'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT calculation.
!> \param qs_env  Quickstep environment
!> \param calc_forces ...
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 refactored to be used with Davidson eigensolver [Sergey Chulkov]
!>    * 03.2017 cleaned and refactored [Sergey Chulkov]
!> \note Based on the subroutines tddfpt_env_init(), and tddfpt_env_deallocate().
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt'

      CHARACTER(len=20)                                  :: nstates_str
      INTEGER                                            :: energy_unit, handle, ispin, istate, &
                                                            iter, log_unit, mult, my_state, nao, &
                                                            niters, nspins, nstates, nstates_read, &
                                                            old_state
      INTEGER, DIMENSION(maxspins)                       :: nactive
      LOGICAL                                            :: do_admm, do_exck, do_hfx, do_hfxlr, &
                                                            is_restarted, state_change
      REAL(kind=dp)                                      :: conv
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dipole_op_mos_occ, evects, S_evects
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_oep, matrix_s, &
                                                            matrix_s_aux_fit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env, kernel_env_admm_aux
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: lri_section, tddfpt_print_section, &
                                                            tddfpt_section, xc_section
      TYPE(stda_env_type), TARGET                        :: stda_kernel
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! input section print/xc
      NULLIFY (tddfpt_section)
      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
      CALL tddfpt_input(qs_env, do_hfx, do_admm, do_exck, do_hfxlr, &
                        xc_section, tddfpt_print_section, lri_section)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, dft_control=dft_control, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, scf_env=scf_env)
      tddfpt_control => dft_control%tddfpt2_control
      tddfpt_control%do_hfx = do_hfx
      tddfpt_control%do_admm = do_admm
      tddfpt_control%do_hfxlr = do_hfxlr
      tddfpt_control%do_exck = do_exck

      CALL cite_reference(Iannuzzi2005)
      IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         CALL cite_reference(Grimme2013)
         CALL cite_reference(Grimme2016)
      END IF

      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "PROGRAM_BANNER", extension=".tddfptLog")
      CALL tddfpt_header(log_unit)
      CALL kernel_info(log_unit, dft_control, tddfpt_control, xc_section)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "PROGRAM_BANNER")

      ! obtain occupied and virtual (unoccupied) ground-state Kohn-Sham orbitals
      NULLIFY (gs_mos)
      CALL tddfpt_init_mos(qs_env, gs_mos)

      ! obtain corrected KS-matrix
      CALL tddfpt_oecorr(qs_env, gs_mos, matrix_ks_oep)

      IF ((tddfpt_control%do_lrigpw) .AND. (tddfpt_control%kernel /= tddfpt_kernel_full)) THEN
         CALL cp_abort(__LOCATION__, "LRI only implemented for full kernel")
      END IF

      NULLIFY (admm_env)
      IF (do_admm) CALL get_qs_env(qs_env, admm_env=admm_env)

      IF (ASSOCIATED(matrix_ks_oep)) matrix_ks => matrix_ks_oep

      ! components of the dipole operator
      CALL tddfpt_dipole_operator(dipole_op_mos_occ, tddfpt_control, gs_mos, qs_env)

      nspins = SIZE(gs_mos)
      ! multiplicity of molecular system
      IF (nspins > 1) THEN
         mult = ABS(SIZE(gs_mos(1)%evals_occ) - SIZE(gs_mos(2)%evals_occ)) + 1
         IF (mult > 2) &
            CALL cp_warn(__LOCATION__, "There is a convergence issue for multiplicity >= 3")
      ELSE
         IF (tddfpt_control%rks_triplets) THEN
            mult = 3
         ELSE
            mult = 1
         END IF
      END IF

      ! split mpi communicator
      ALLOCATE (evects(nspins, 1))
      DO ispin = 1, nspins
         evects(ispin, 1)%matrix => gs_mos(ispin)%mos_occ
      END DO
      CALL tddfpt_sub_env_init(sub_env, qs_env, mos_occ=evects(:, 1), kernel=tddfpt_control%kernel)
      DEALLOCATE (evects)

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         ! create environment for Full Kernel
         IF (dft_control%qs_control%xtb) THEN
            CPABORT("TDDFPT: xTB only works with sTDA Kernel")
         END IF

         ! allocate pools and work matrices
         nstates = tddfpt_control%nstates
         CALL tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, do_admm, qs_env, sub_env)

         CALL tddfpt_construct_ground_state_orb_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                        is_rks_triplets=tddfpt_control%rks_triplets, &
                                                        qs_env=qs_env, sub_env=sub_env, &
                                                        wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub)

         tddfpt_control%do_response_now = .TRUE. ! needs to be after tddfpt_construct_ground_state_orb_density
      END IF

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         IF (do_admm) THEN
            ! Full kernel with ADMM
            CALL create_kernel_env(kernel_env=full_kernel_env, &
                                   rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                   xc_section=admm_env%xc_section_primary, &
                                   is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)

            CALL tddfpt_construct_aux_fit_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                  rho_aux_fit_struct=work_matrices%rho_aux_fit_struct_sub, &
                                                  qs_env=qs_env, sub_env=sub_env, &
                                                  wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub, &
                                                  wfm_rho_aux_fit=work_matrices%rho_ao_aux_fit_fm_sub, &
                                                  wfm_aux_orb=work_matrices%wfm_aux_orb_sub)

            CALL create_kernel_env(kernel_env=kernel_env_admm_aux, &
                                   rho_struct_sub=work_matrices%rho_aux_fit_struct_sub, &
                                   xc_section=admm_env%xc_section_aux, &
                                   is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)
            kernel_env%full_kernel => full_kernel_env
            kernel_env%admm_kernel => kernel_env_admm_aux
         ELSE
            ! Full kernel
            CALL create_kernel_env(kernel_env=full_kernel_env, &
                                   rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                   xc_section=xc_section, &
                                   is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)
            kernel_env%full_kernel => full_kernel_env
            NULLIFY (kernel_env%admm_kernel)
         END IF
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         ! sTDA kernel
         nactive = 0
         CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)
         DO ispin = 1, SIZE(gs_mos)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nactive(ispin))
         END DO
         CALL allocate_stda_env(qs_env, stda_kernel, nao, nactive)
         ! sTDA parameters
         CALL stda_init_param(qs_env, stda_kernel, tddfpt_control%stda_control)
         ! allocate pools and work matrices
         nstates = tddfpt_control%nstates
         CALL stda_create_work_matrices(work_matrices, gs_mos, nstates, qs_env, sub_env)
         !
         CALL stda_init_matrices(qs_env, stda_kernel, sub_env, work_matrices, tddfpt_control)
         !
         kernel_env%stda_kernel => stda_kernel
         NULLIFY (kernel_env%full_kernel)
         NULLIFY (kernel_env%admm_kernel)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_none) THEN
         ! allocate pools and work matrices
         nstates = tddfpt_control%nstates
         CALL stda_create_work_matrices(work_matrices, gs_mos, nstates, qs_env, sub_env)
         NULLIFY (kernel_env%full_kernel)
         NULLIFY (kernel_env%admm_kernel)
         NULLIFY (kernel_env%stda_kernel)
      ELSE
         CPABORT('Unknown kernel type')
      END IF

      IF (tddfpt_control%do_lrigpw) THEN
         CALL tddfpt2_lri_init(qs_env, kernel_env, lri_section, tddfpt_print_section)
      END IF

      ALLOCATE (evals(nstates))
      ALLOCATE (evects(nspins, nstates), S_evects(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (evects(ispin, istate)%matrix, S_evects(ispin, istate)%matrix)
            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, S_evects(ispin, istate)%matrix)
         END DO
      END DO

      ! reuse Ritz vectors from the previous calculation if available
      IF (tddfpt_control%is_restart) THEN
         nstates_read = tddfpt_read_restart(evects=evects, evals=evals, gs_mos=gs_mos, &
                                            logger=logger, tddfpt_section=tddfpt_section, &
                                            tddfpt_print_section=tddfpt_print_section, &
                                            fm_pool_ao_mo_occ=work_matrices%fm_pool_ao_mo_occ, &
                                            blacs_env_global=blacs_env)
      ELSE
         nstates_read = 0
      END IF

      is_restarted = nstates_read >= nstates

      ! build the list of missed singly excited states and sort them in ascending order
      ! according to their excitation energies
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "GUESS_VECTORS", extension=".tddfptLog")
      CALL tddfpt_guess_vectors(evects=evects, evals=evals, gs_mos=gs_mos, log_unit=log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "GUESS_VECTORS")

      CALL tddfpt_orthogonalize_psi1_psi0(evects, work_matrices%S_C0_C0T)
      CALL tddfpt_orthonormalize_psi1_psi1(evects, nstates, S_evects, matrix_s(1)%matrix)

      niters = tddfpt_control%niters
      IF (niters > 0) THEN
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "ITERATION_INFO", extension=".tddfptLog")
         energy_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "DETAILED_ENERGY", extension=".tddfptLog")

         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                      TDDFPT WAVEFUNCTION OPTIMIZATION                       -", &
               "-------------------------------------------------------------------------------"

            WRITE (log_unit, '(/,T11,A,T27,A,T40,A,T62,A)') "Step", "Time", "Convergence", "Conv. states"
            WRITE (log_unit, '(1X,79("-"))')
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         DO
            ! *** perform Davidson iterations ***
            conv = tddfpt_davidson_solver(evects=evects, evals=evals, S_evects=S_evects, gs_mos=gs_mos, &
                                          tddfpt_control=tddfpt_control, &
                                          matrix_ks=matrix_ks, qs_env=qs_env, &
                                          kernel_env=kernel_env, &
                                          sub_env=sub_env, logger=logger, &
                                          iter_unit=log_unit, energy_unit=energy_unit, &
                                          tddfpt_print_section=tddfpt_print_section, &
                                          work_matrices=work_matrices)

            ! at this point at least one of the following conditions are met:
            ! a) convergence criteria has been achieved;
            ! b) maximum number of iterations has been reached;
            ! c) Davidson iterations must be restarted due to lack of Krylov vectors or numerical instability

            CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)
            ! terminate the loop if either (a) or (b) is true ...
            IF ((conv <= tddfpt_control%conv .AND. is_restarted) .OR. iter >= niters) EXIT

            ! ... otherwise restart Davidson iterations
            is_restarted = .TRUE.
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(1X,25("-"),1X,A,1X,25("-"))') "Restart Davidson iterations"
               CALL m_flush(log_unit)
            END IF
         END DO

         ! write TDDFPT restart file at the last iteration if requested to do so
         CALL cp_iterate(logger%iter_info, increment=0, last=.TRUE.)
         CALL tddfpt_write_restart(evects=evects, evals=evals, gs_mos=gs_mos, &
                                   logger=logger, tddfpt_print_section=tddfpt_print_section)

         CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

         ! print convergence summary
         IF (log_unit > 0) THEN
            CALL integer_to_string(iter, nstates_str)
            IF (conv <= tddfpt_control%conv) THEN
               WRITE (log_unit, "(1X,A)") "", &
                  "-------------------------------------------------------------------------------", &
                  "-  TDDFPT run converged in "//TRIM(nstates_str)//" iteration(s) ", &
                  "-------------------------------------------------------------------------------"
            ELSE
               WRITE (log_unit, "(1X,A)") "", &
                  "-------------------------------------------------------------------------------", &
                  "-  TDDFPT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ", &
                  "-------------------------------------------------------------------------------"
            END IF
         END IF

         CALL cp_print_key_finished_output(energy_unit, logger, tddfpt_print_section, "DETAILED_ENERGY")
         CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "ITERATION_INFO")
      ELSE
         CALL cp_warn(__LOCATION__, "Skipping TDDFPT wavefunction optimization")
      END IF

      CALL tddfpt_write_newtonx_output(evects, evals, gs_mos, logger, tddfpt_print_section, qs_env)

      IF (ASSOCIATED(matrix_ks_oep) .AND. tddfpt_control%dipole_form == tddfpt_dipole_velocity) THEN
         CALL cp_warn(__LOCATION__, &
                      "Transition dipole moments and oscillator strengths are likely to be incorrect "// &
                      "when computed using an orbital energy correction XC-potential together with "// &
                      "the velocity form of dipole transition integrals")
      END IF

      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit()
      CALL tddfpt_print_summary(log_unit, evects, evals, mult, dipole_op_mos_occ, tddfpt_control%dipole_form)
      CALL tddfpt_print_excitation_analysis(log_unit, evects, evals, gs_mos, matrix_s(1)%matrix, &
                                            min_amplitude=tddfpt_control%min_excitation_amplitude)
      CALL tddfpt_print_nto_analysis(qs_env, evects, evals, gs_mos, matrix_s(1)%matrix, &
                                     tddfpt_print_section)

      IF (tddfpt_control%do_lrigpw) THEN
         CALL lri_print_stat(qs_env, ltddfpt=.TRUE., tddfpt_lri_env=kernel_env%full_kernel%lri_env)
      END IF

      ! excited state potential energy surface
      IF (qs_env%excited_state) THEN
         IF (ASSOCIATED(matrix_ks_oep)) CALL get_qs_env(qs_env, matrix_ks=matrix_ks) ! for gradients unshifted KS matrix
         CALL get_qs_env(qs_env, exstate_env=ex_env)
         state_change = .FALSE.
         IF (ex_env%state > 0) THEN
            my_state = ex_env%state
         ELSEIF (ex_env%state < 0) THEN
            ! state following
            my_state = ABS(ex_env%state)
            CALL assign_state(qs_env, matrix_s, evects, mos, ex_env%fingerprint, my_state)
            IF (my_state /= ABS(ex_env%state)) THEN
               state_change = .TRUE.
               old_state = ABS(ex_env%state)
            END IF
            ex_env%state = -my_state
         ELSE
            CALL cp_warn(__LOCATION__, "Active excited state not assigned. Use the first state.")
            my_state = 1
         END IF
         CPASSERT(my_state > 0)
         IF (my_state > nstates) THEN
            CALL cp_warn(__LOCATION__, "There were not enough excited states calculated.")
            CPABORT("excited state potential energy surface")
         END IF
         !
         ! energy
         ex_env%evalue = evals(my_state)
         ! excitation vector
         IF (ASSOCIATED(ex_env%evect)) THEN
            DO ispin = 1, SIZE(ex_env%evect)
               CALL cp_fm_release(ex_env%evect(ispin)%matrix)
            END DO
            DEALLOCATE (ex_env%evect)
         END IF
         ALLOCATE (ex_env%evect(nspins))
         DO ispin = 1, nspins
            CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, matrix_struct=matrix_struct)
            NULLIFY (ex_env%evect(ispin)%matrix)
            CALL cp_fm_create(ex_env%evect(ispin)%matrix, matrix_struct)
            CALL cp_fm_to_fm(evects(ispin, my_state)%matrix, ex_env%evect(ispin)%matrix)
         END DO

         IF (log_unit > 0) THEN
            IF (state_change) THEN
               WRITE (log_unit, "(1X,A,I4,A,I4)") "Target state has been changed from ", &
                  old_state, " to new state ", my_state
            END IF
            WRITE (log_unit, "(1X,A,I4,A,F12.5,A)") "Calculate properties for state:", &
               my_state, "      with excitation energy ", ex_env%evalue*evolt, " eV"
         END IF

         IF (calc_forces) THEN
            ! rhs of linres equation
            IF (ASSOCIATED(ex_env%cpmos)) THEN
               DO ispin = 1, SIZE(ex_env%cpmos)
                  CALL cp_fm_release(ex_env%cpmos(ispin)%matrix)
               END DO
               DEALLOCATE (ex_env%cpmos)
            END IF
            ALLOCATE (ex_env%cpmos(nspins))
            DO ispin = 1, nspins
               CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, matrix_struct=matrix_struct)
               NULLIFY (ex_env%cpmos(ispin)%matrix)
               CALL cp_fm_create(ex_env%cpmos(ispin)%matrix, matrix_struct)
               CALL cp_fm_set_all(ex_env%cpmos(ispin)%matrix, 0.0_dp)
            END DO
            CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
            CALL dbcsr_allocate_matrix_set(ex_env%matrix_pe, nspins)
            DO ispin = 1, nspins
               ALLOCATE (ex_env%matrix_pe(ispin)%matrix)
               CALL dbcsr_create(ex_env%matrix_pe(ispin)%matrix, template=matrix_s(1)%matrix)
               CALL dbcsr_copy(ex_env%matrix_pe(ispin)%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(ex_env%matrix_pe(ispin)%matrix, 0.0_dp)
               CALL tddfpt_resvec1(ex_env%evect(ispin)%matrix, gs_mos(ispin)%mos_occ, &
                                   matrix_s(1)%matrix, ex_env%matrix_pe(ispin)%matrix)
            END DO
            !
            ! ground state ADMM!
            IF (dft_control%do_admm) THEN
               CALL get_qs_env(qs_env, admm_env=admm_env, matrix_s_aux_fit=matrix_s_aux_fit)
               CALL dbcsr_allocate_matrix_set(ex_env%matrix_pe_admm, nspins)
               DO ispin = 1, nspins
                  ALLOCATE (ex_env%matrix_pe_admm(ispin)%matrix)
                  CALL dbcsr_create(ex_env%matrix_pe_admm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_copy(ex_env%matrix_pe_admm(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_set(ex_env%matrix_pe_admm(ispin)%matrix, 0.0_dp)
                  CALL tddfpt_resvec1_admm(ex_env%matrix_pe(ispin)%matrix, &
                                           admm_env, ex_env%matrix_pe_admm(ispin)%matrix)
               END DO
            END IF
            !
            CALL dbcsr_allocate_matrix_set(ex_env%matrix_hz, nspins)
            DO ispin = 1, nspins
               ALLOCATE (ex_env%matrix_hz(ispin)%matrix)
               CALL dbcsr_create(ex_env%matrix_hz(ispin)%matrix, template=matrix_s(1)%matrix)
               CALL dbcsr_copy(ex_env%matrix_hz(ispin)%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(ex_env%matrix_hz(ispin)%matrix, 0.0_dp)
            END DO
            IF (dft_control%qs_control%xtb) THEN
               CALL tddfpt_resvec2_xtb(qs_env, ex_env%matrix_pe, gs_mos, ex_env%matrix_hz, ex_env%cpmos)
            ELSE
               CALL tddfpt_resvec2(qs_env, ex_env%matrix_pe, ex_env%matrix_pe_admm, &
                                   gs_mos, ex_env%matrix_hz, ex_env%cpmos)
            END IF
            !
            CALL dbcsr_allocate_matrix_set(ex_env%matrix_px1, nspins)
            DO ispin = 1, nspins
               ALLOCATE (ex_env%matrix_px1(ispin)%matrix)
               CALL dbcsr_create(ex_env%matrix_px1(ispin)%matrix, template=matrix_s(1)%matrix)
               CALL dbcsr_copy(ex_env%matrix_px1(ispin)%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(ex_env%matrix_px1(ispin)%matrix, 0.0_dp)
            END DO
            ! Kernel ADMM
            IF (tddfpt_control%do_admm) THEN
               CALL get_qs_env(qs_env, admm_env=admm_env, matrix_s_aux_fit=matrix_s_aux_fit)
               CALL dbcsr_allocate_matrix_set(ex_env%matrix_px1_admm, nspins)
               DO ispin = 1, nspins
                  ALLOCATE (ex_env%matrix_px1_admm(ispin)%matrix)
                  CALL dbcsr_create(ex_env%matrix_px1_admm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_copy(ex_env%matrix_px1_admm(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_set(ex_env%matrix_px1_admm(ispin)%matrix, 0.0_dp)
               END DO
            END IF
            ! TDA forces
            CALL tddfpt_forces(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices)
            ! Rotate res vector cpmos into original frame of occupied orbitals
            CALL tddfpt_resvec3(qs_env, ex_env%cpmos, work_matrices)
         END IF
      END IF

      ! clean up
      DO istate = 1, SIZE(evects, 2)
         DO ispin = 1, nspins
            CALL cp_fm_release(evects(ispin, istate)%matrix)
            CALL cp_fm_release(S_evects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (evects, S_evects, evals)

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         IF (do_admm) CALL release_kernel_env(kernel_env%admm_kernel)
         IF (tddfpt_control%do_lrigpw) THEN
            CALL lri_env_release(kernel_env%full_kernel%lri_env)
            CALL lri_density_release(kernel_env%full_kernel%lri_density)
         END IF
         CALL release_kernel_env(kernel_env%full_kernel)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         CALL deallocate_stda_env(stda_kernel)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_none) THEN
         !
      ELSE
         CPABORT('Unknown kernel type')
      END IF
      CALL tddfpt_release_work_matrices(work_matrices, sub_env)
      CALL tddfpt_sub_env_release(sub_env)

      IF (ALLOCATED(dipole_op_mos_occ)) THEN
         DO ispin = nspins, 1, -1
            DO istate = SIZE(dipole_op_mos_occ, 1), 1, -1
               CALL cp_fm_release(dipole_op_mos_occ(istate, ispin)%matrix)
            END DO
         END DO
         DEALLOCATE (dipole_op_mos_occ)
      END IF

      DO ispin = nspins, 1, -1
         CALL tddfpt_release_ground_state_mos(gs_mos(ispin))
      END DO
      DEALLOCATE (gs_mos)

      IF (ASSOCIATED(matrix_ks_oep)) &
         CALL dbcsr_deallocate_matrix_set(matrix_ks_oep)

      CALL timestop(handle)

   END SUBROUTINE tddfpt

! **************************************************************************************************
!> \brief TDDFPT input
!> \param qs_env  Quickstep environment
!> \param do_hfx ...
!> \param do_admm ...
!> \param do_exck ...
!> \param do_hfxlr ...
!> \param xc_section ...
!> \param tddfpt_print_section ...
!> \param lri_section ...
! **************************************************************************************************
   SUBROUTINE tddfpt_input(qs_env, do_hfx, do_admm, do_exck, do_hfxlr, &
                           xc_section, tddfpt_print_section, lri_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(INOUT)                             :: do_hfx, do_admm, do_exck, do_hfxlr
      TYPE(section_vals_type), POINTER                   :: xc_section, tddfpt_print_section, &
                                                            lri_section

      CHARACTER(len=20)                                  :: nstates_str
      LOGICAL                                            :: exar, exf, exgcp, exhf, exhfxk, exk, &
                                                            explicit_root, expot, exvdw, exwfn
      REAL(kind=dp)                                      :: C_hf
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, tddfpt_section, &
                                                            xc_root, xc_sub
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      NULLIFY (dft_control, input)
      CALL get_qs_env(qs_env, dft_control=dft_control, input=input)
      tddfpt_control => dft_control%tddfpt2_control

      ! no k-points
      CPASSERT(dft_control%nimages <= 1)

      IF (tddfpt_control%nstates <= 0) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL cp_warn(__LOCATION__, "TDDFPT calculation was requested for "// &
                      TRIM(nstates_str)//" excited states: nothing to do.")
         RETURN
      END IF

      NULLIFY (tddfpt_section, tddfpt_print_section)
      tddfpt_section => section_vals_get_subs_vals(input, "PROPERTIES%TDDFPT")
      tddfpt_print_section => section_vals_get_subs_vals(tddfpt_section, "PRINT")

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         NULLIFY (xc_root)
         xc_root => section_vals_get_subs_vals(tddfpt_section, "XC")
         CALL section_vals_get(xc_root, explicit=explicit_root)
         NULLIFY (xc_section)
         IF (explicit_root) THEN
            ! No ADIABATIC_RESCALING option possible
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "ADIABATIC_RESCALING")
            CALL section_vals_get(xc_sub, explicit=exar)
            IF (exar) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel with ADIABATIC_RESCALING not possible.")
               CPABORT("TDDFPT Input")
            END IF
            ! No GCP_POTENTIAL option possible
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "GCP_POTENTIAL")
            CALL section_vals_get(xc_sub, explicit=exgcp)
            IF (exgcp) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel with GCP_POTENTIAL not possible.")
               CPABORT("TDDFPT Input")
            END IF
            ! No VDW_POTENTIAL option possible
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "VDW_POTENTIAL")
            CALL section_vals_get(xc_sub, explicit=exvdw)
            IF (exvdw) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel with VDW_POTENTIAL not possible.")
               CPABORT("TDDFPT Input")
            END IF
            ! No WF_CORRELATION option possible
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "WF_CORRELATION")
            CALL section_vals_get(xc_sub, explicit=exwfn)
            IF (exwfn) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel with WF_CORRELATION not possible.")
               CPABORT("TDDFPT Input")
            END IF
            ! No XC_POTENTIAL option possible
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "XC_POTENTIAL")
            CALL section_vals_get(xc_sub, explicit=expot)
            IF (expot) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel with XC_POTENTIAL not possible.")
               CPABORT("TDDFPT Input")
            END IF
            !
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "XC_FUNCTIONAL")
            CALL section_vals_get(xc_sub, explicit=exf)
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "XC_KERNEL")
            CALL section_vals_get(xc_sub, explicit=exk)
            IF ((exf .AND. exk) .OR. .NOT. (exf .OR. exk)) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel needs XC_FUNCTIONAL or XC_KERNEL section.")
               CPABORT("TDDFPT Input")
            END IF
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "HF")
            CALL section_vals_get(xc_sub, explicit=exhf)
            NULLIFY (xc_sub)
            xc_sub => section_vals_get_subs_vals(xc_root, "HFX_KERNEL")
            CALL section_vals_get(xc_sub, explicit=exhfxk)
            IF (exhf .AND. exhfxk) THEN
               CALL cp_warn(__LOCATION__, "TDDFPT Kernel: HF and HFX_KERNEL are exclusive.")
               CPABORT("TDDFPT Input")
            END IF
            !
            xc_section => xc_root
            hfx_section => section_vals_get_subs_vals(xc_section, "HF")
            CALL section_vals_get(hfx_section, explicit=do_hfx)
            IF (do_hfx) THEN
               CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
               do_hfx = (C_hf /= 0.0_dp)
            END IF
            do_admm = do_hfx .AND. dft_control%do_admm
            IF (do_admm) THEN
               ! 'admm_env%xc_section_primary' and 'admm_env%xc_section_aux' need to be redefined
               CALL cp_abort(__LOCATION__, &
                             "ADMM is not implemented for a TDDFT kernel XC-functional which is different from "// &
                             "the one used for the ground-state calculation. A ground-state 'admm_env' cannot be reused.")
            END IF
            ! SET HFX_KERNEL and/or XC_KERNEL
            IF (exhfxk) THEN
               do_exck = .TRUE.
            ELSE
               do_exck = .FALSE.
            END IF
            IF (exk) THEN
               do_hfxlr = .TRUE.
            ELSE
               do_hfxlr = .FALSE.
            END IF
         ELSE
            xc_section => section_vals_get_subs_vals(input, "DFT%XC")
            hfx_section => section_vals_get_subs_vals(xc_section, "HF")
            CALL section_vals_get(hfx_section, explicit=do_hfx)
            IF (do_hfx) THEN
               CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
               do_hfx = (C_hf /= 0.0_dp)
            END IF
            do_admm = do_hfx .AND. dft_control%do_admm
            do_exck = .FALSE.
            do_hfxlr = .FALSE.
         END IF
      ELSE
         do_hfx = .FALSE.
         do_admm = .FALSE.
         do_exck = .FALSE.
         do_hfxlr = .FALSE.
      END IF

      ! reset rks_triplets if UKS is in use
      IF (tddfpt_control%rks_triplets .AND. dft_control%nspins > 1) THEN
         tddfpt_control%rks_triplets = .FALSE.
         CALL cp_warn(__LOCATION__, "Keyword RKS_TRIPLETS has been ignored for spin-polarised calculations")
      END IF

      ! lri input
      IF (tddfpt_control%do_lrigpw) THEN
         lri_section => section_vals_get_subs_vals(tddfpt_section, "LRIGPW")
      END IF

   END SUBROUTINE tddfpt_input

! **************************************************************************************************
!> \brief ...
!> \param log_unit ...
!> \param dft_control ...
!> \param tddfpt_control ...
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE kernel_info(log_unit, dft_control, tddfpt_control, xc_section)
      INTEGER, INTENT(IN)                                :: log_unit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=4)                                   :: ktype
      LOGICAL                                            :: lsd

      lsd = (dft_control%nspins > 1)
      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         ktype = "FULL"
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(T2,A,T77,A4)") "KERNEL|", TRIM(ktype)
            CALL xc_write(log_unit, xc_section, lsd)
            IF (tddfpt_control%do_hfx) THEN
               IF (tddfpt_control%do_admm) THEN
                  WRITE (log_unit, "(T2,A,T62,A19)") "KERNEL|", "ADMM exact exchange"
                  IF (tddfpt_control%admm_symm) THEN
                     WRITE (log_unit, "(T2,A,T60,A21)") "KERNEL|", "symmetric ADMM kernel"
                  END IF
               ELSE
                  WRITE (log_unit, "(T2,A,T67,A14)") "KERNEL|", "exact exchange"
               END IF
            END IF
            IF (tddfpt_control%do_lrigpw) THEN
               WRITE (log_unit, "(T2,A,T43,A38)") "KERNEL|", "LRI approximation of transition density"
            END IF
         END IF
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         ktype = "sTDA"
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(T2,A,T77,A4)") "KERNEL|", TRIM(ktype)
            IF (tddfpt_control%stda_control%do_ewald) THEN
               WRITE (log_unit, "(T2,A,T78,A3)") "KERNEL| Coulomb term uses Ewald summation"
            ELSE
               WRITE (log_unit, "(T2,A,T78,A3)") "KERNEL| Coulomb term uses direct summation (MIC)"
            END IF
            IF (tddfpt_control%stda_control%do_exchange) THEN
               WRITE (log_unit, "(T2,A,T78,A3)") "KERNEL| Exact exchange term", "YES"
               WRITE (log_unit, "(T2,A,T71,F10.3)") "KERNEL| Short range HFX fraction:", &
                  tddfpt_control%stda_control%hfx_fraction
            ELSE
               WRITE (log_unit, "(T2,A,T79,A2)") "KERNEL| Exact exchange term", "NO"
            END IF
            WRITE (log_unit, "(T2,A,T66,E15.3)") "KERNEL| Transition density filter", &
               tddfpt_control%stda_control%eps_td_filter
         END IF
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_none) THEN
         ktype = "NONE"
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(T2,A,T77,A4)") "KERNEL|", TRIM(ktype)
         END IF
      ELSE
         !CPABORT("Unknown kernel")
      END IF
      !
      IF (log_unit > 0) THEN
         IF (tddfpt_control%rks_triplets) THEN
            WRITE (log_unit, "(T2,A,T74,A7)") "KERNEL| Spin symmetry of excitations", "Triplet"
         ELSE IF (lsd) THEN
            WRITE (log_unit, "(T2,A,T69,A12)") "KERNEL| Spin symmetry of excitations", "Unrestricted"
         ELSE
            WRITE (log_unit, "(T2,A,T74,A7)") "KERNEL| Spin symmetry of excitations", "Singlet"
         END IF
         WRITE (log_unit, "(T2,A,T73,I8)") "TDDFPT| Number of states calculated", tddfpt_control%nstates
         WRITE (log_unit, "(T2,A,T73,I8)") "TDDFPT| Number of Davidson iterations", tddfpt_control%niters
         WRITE (log_unit, "(T2,A,T66,E15.3)") "TDDFPT| Davidson iteration convergence", tddfpt_control%conv
         WRITE (log_unit, "(T2,A,T73,I8)") "TDDFPT| Max. number of Krylov space vectors", tddfpt_control%nkvs
      END IF

   END SUBROUTINE kernel_info

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_s ...
!> \param evects ...
!> \param mos ...
!> \param fingerprint ...
!> \param my_state ...
! **************************************************************************************************
   SUBROUTINE assign_state(qs_env, matrix_s, evects, mos, fingerprint, my_state)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_p_type), DIMENSION(:, :)                :: evects
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fingerprint
      INTEGER, INTENT(INOUT)                             :: my_state

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'assign_state'

      INTEGER                                            :: handle, is, ispin, natom, nspins, nstate
      REAL(KIND=dp)                                      :: dvm, dvp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dv
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chr
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: charges
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: psi0
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p1

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, natom=natom, para_env=para_env)
      nspins = SIZE(mos)
      nstate = SIZE(evects, 2)
      !
      NULLIFY (matrix_p1)
      CALL dbcsr_allocate_matrix_set(matrix_p1, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_p1(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_p1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(matrix_p1(ispin)%matrix, matrix_s(1)%matrix)
      END DO
      ALLOCATE (psi0(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         psi0(ispin)%matrix => mo_coeff
      END DO
      !
      IF (ASSOCIATED(fingerprint)) THEN
         CPASSERT(SIZE(fingerprint, 1) == natom)
         CPASSERT(SIZE(fingerprint, 2) == nspins)
         ALLOCATE (charges(natom, nspins, nstate))
         !
         DO is = 1, nstate
            DO ispin = 1, nspins
               CALL dbcsr_set(matrix_p1(ispin)%matrix, 0.0_dp)
            END DO
            chr => charges(:, :, is)
            CALL build_dm_response(psi0, evects(:, is), matrix_p1)
            CALL compute_charges(matrix_p1, matrix_s(1)%matrix, chr, para_env)
         END DO
         !
         ALLOCATE (dv(nstate))
         DO is = 1, nstate
            dvp = SQRT(SUM((fingerprint(:, :) + charges(:, :, is))**2)/REAL(natom, dp))
            dvm = SQRT(SUM((fingerprint(:, :) - charges(:, :, is))**2)/REAL(natom, dp))
            dv(is) = MIN(dvp, dvm)
         END DO
         !
         my_state = MINVAL(MINLOC(dv))
         fingerprint(:, :) = charges(:, :, my_state)
         !
         DEALLOCATE (charges, dv)
      ELSE
         ALLOCATE (fingerprint(natom, nspins))
         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_p1(ispin)%matrix, 0.0_dp)
         END DO
         CALL build_dm_response(psi0, evects(:, my_state), matrix_p1)
         CALL compute_charges(matrix_p1, matrix_s(1)%matrix, fingerprint, para_env)
      END IF
      !
      CALL dbcsr_deallocate_matrix_set(matrix_p1)
      DEALLOCATE (psi0)

      CALL timestop(handle)

   END SUBROUTINE assign_state

END MODULE qs_tddfpt2_methods

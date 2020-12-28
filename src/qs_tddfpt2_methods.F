!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE header,                          ONLY: tddfpt_header
   USE input_constants,                 ONLY: tddfpt_dipole_velocity,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_eigensolver,          ONLY: tddfpt_davidson_solver,&
                                              tddfpt_orthogonalize_psi1_psi0,&
                                              tddfpt_orthonormalize_psi1_psi1
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_excitation_analysis,&
                                              tddfpt_print_nto_analysis,&
                                              tddfpt_print_summary
   USE qs_tddfpt2_restart,              ONLY: tddfpt_read_restart,&
                                              tddfpt_write_restart
   USE qs_tddfpt2_stda_types,           ONLY: allocate_stda_env,&
                                              deallocate_stda_env,&
                                              stda_env_type,&
                                              stda_init_param
   USE qs_tddfpt2_stda_utils,           ONLY: stda_init_matrices
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_sub_env_init,&
                                              tddfpt_sub_env_release,&
                                              tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: &
        full_kernel_env_type, stda_create_work_matrices, tddfpt_create_kernel_env, &
        tddfpt_create_work_matrices, tddfpt_ground_state_mos, tddfpt_kernel_env_type, &
        tddfpt_release_kernel_env, tddfpt_release_work_matrices, tddfpt_work_matrices
   USE qs_tddfpt2_utils,                ONLY: tddfpt_guess_vectors,&
                                              tddfpt_init_mos,&
                                              tddfpt_oecorr,&
                                              tddfpt_release_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
#include "base/base_uses.f90"

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
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 refactored to be used with Davidson eigensolver [Sergey Chulkov]
!>    * 03.2017 cleaned and refactored [Sergey Chulkov]
!> \note Based on the subroutines tddfpt_env_init(), and tddfpt_env_deallocate().
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt'

      CHARACTER(len=20)                                  :: nstates_str
      INTEGER                                            :: energy_unit, handle, ispin, istate, &
                                                            iter, log_unit, mult, nao, niters, &
                                                            nspins, nstates, nstates_read
      INTEGER, DIMENSION(maxspins)                       :: nactive
      LOGICAL                                            :: do_admm, do_hfx, is_restarted
      REAL(kind=dp)                                      :: conv
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dipole_op_mos_occ, evects, S_evects
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_oep, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env, kernel_env_admm_aux
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section, tddfpt_section, &
                                                            xc_section
      TYPE(stda_env_type), TARGET                        :: stda_kernel
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_kernel_env_type)                       :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! input section print/xc
      CALL tddfpt_input(qs_env, do_hfx, do_admm, xc_section, tddfpt_print_section)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, dft_control=dft_control, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, scf_env=scf_env)
      tddfpt_control => dft_control%tddfpt2_control

      CALL cite_reference(Iannuzzi2005)
      IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         CALL cite_reference(Grimme2013)
         CALL cite_reference(Grimme2016)
      END IF

      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "PROGRAM_BANNER", extension=".tddfptLog")
      CALL tddfpt_header(log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "PROGRAM_BANNER")

      ! obtain occupied and virtual (unoccupied) ground-state Kohn-Sham orbitals
      NULLIFY (gs_mos)
      CALL tddfpt_init_mos(qs_env, gs_mos)

      ! obtain corrected KS-matrix
      CALL tddfpt_oecorr(qs_env, gs_mos, matrix_ks_oep, do_hfx)

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
         CALL tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, qs_env, sub_env)

         CALL tddfpt_construct_ground_state_orb_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                        is_rks_triplets=tddfpt_control%rks_triplets, &
                                                        qs_env=qs_env, sub_env=sub_env, &
                                                        wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub)

         IF (do_admm) THEN
            ! Full kernel with ADMM
            CALL tddfpt_create_kernel_env(kernel_env=full_kernel_env, &
                                          rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                          xc_section=admm_env%xc_section_primary, &
                                          is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)

            CALL tddfpt_construct_aux_fit_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                  rho_aux_fit_struct=work_matrices%rho_aux_fit_struct_sub, &
                                                  qs_env=qs_env, sub_env=sub_env, &
                                                  wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub, &
                                                  wfm_rho_aux_fit=work_matrices%rho_ao_aux_fit_fm_sub, &
                                                  wfm_aux_orb=work_matrices%wfm_aux_orb_sub)

            CALL tddfpt_create_kernel_env(kernel_env=kernel_env_admm_aux, &
                                          rho_struct_sub=work_matrices%rho_aux_fit_struct_sub, &
                                          xc_section=admm_env%xc_section_aux, &
                                          is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)
            kernel_env%full_kernel => full_kernel_env
            kernel_env%admm_kernel => kernel_env_admm_aux
         ELSE
            ! Full kernel
            CALL tddfpt_create_kernel_env(kernel_env=full_kernel_env, &
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
      ELSE
         CPABORT('Unknown kernel type')
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
                                          do_hfx=do_hfx, tddfpt_control=tddfpt_control, &
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

      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit()

      IF (ASSOCIATED(matrix_ks_oep) .AND. tddfpt_control%dipole_form == tddfpt_dipole_velocity) THEN
         CALL cp_warn(__LOCATION__, &
                      "Transition dipole moments and oscillator strengths are likely to be incorrect "// &
                      "when computed using an orbital energy correction XC-potential together with "// &
                      "the velocity form of dipole transition integrals")
      END IF

      CALL tddfpt_print_summary(log_unit, evects, evals, mult, dipole_op_mos_occ)
      CALL tddfpt_print_excitation_analysis(log_unit, evects, evals, gs_mos, matrix_s(1)%matrix, &
                                            min_amplitude=tddfpt_control%min_excitation_amplitude)
      CALL tddfpt_print_nto_analysis(qs_env, evects, evals, gs_mos, matrix_s(1)%matrix, &
                                     tddfpt_print_section)

      ! -- clean up all useless stuff
      DO istate = SIZE(evects, 2), 1, -1
         DO ispin = nspins, 1, -1
            CALL cp_fm_release(evects(ispin, istate)%matrix)
            CALL cp_fm_release(S_evects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (evects, S_evects, evals)

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         IF (do_admm) CALL tddfpt_release_kernel_env(kernel_env%admm_kernel)
         CALL tddfpt_release_kernel_env(kernel_env%full_kernel)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         CALL deallocate_stda_env(stda_kernel)
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
!> \param xc_section ...
!> \param tddfpt_print_section ...
! **************************************************************************************************
   SUBROUTINE tddfpt_input(qs_env, do_hfx, do_admm, xc_section, tddfpt_print_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(INOUT)                             :: do_hfx, do_admm
      TYPE(section_vals_type), POINTER                   :: xc_section, tddfpt_print_section

      CHARACTER(len=20)                                  :: nstates_str
      LOGICAL                                            :: explicit_xc
      REAL(kind=dp)                                      :: C_hf
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, tddfpt_section
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

      NULLIFY (xc_section)
      xc_section => section_vals_get_subs_vals(tddfpt_section, "XC%XC_FUNCTIONAL")
      CALL section_vals_get(xc_section, explicit=explicit_xc)
      IF (explicit_xc) THEN
         xc_section => section_vals_get_subs_vals(tddfpt_section, "XC")
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")

      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
         do_hfx = (C_hf /= 0.0_dp)
      END IF

      do_admm = do_hfx .AND. dft_control%do_admm
      IF (do_admm) THEN
         IF (explicit_xc) THEN
            ! 'admm_env%xc_section_primary' and 'admm_env%xc_section_aux' need to be redefined
            CALL cp_abort(__LOCATION__, &
                          "ADMM is not implemented for a TDDFT kernel XC-functional which is different from "// &
                          "the one used for the ground-state calculation. A ground-state 'admm_env' cannot be reused.")
         END IF
      END IF

      ! reset rks_triplets if UKS is in use
      IF (tddfpt_control%rks_triplets .AND. dft_control%nspins > 1) THEN
         tddfpt_control%rks_triplets = .FALSE.
         CALL cp_warn(__LOCATION__, "Keyword RKS_TRIPLETS has been ignored for spin-polarised calculations")
      END IF

   END SUBROUTINE tddfpt_input

END MODULE qs_tddfpt2_methods

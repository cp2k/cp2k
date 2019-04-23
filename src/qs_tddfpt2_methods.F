!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_methods
   USE admm_types,                      ONLY: admm_type
   USE bibliography,                    ONLY: Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: &
        cholesky_restore, oe_gllb, oe_lb, oe_none, oe_saop, oe_shift, tddfpt_dipole_velocity, &
        tddfpt_kernel_full, tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_vals_create,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_set_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_post_gpw,                 ONLY: make_lumo_gpw
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_eigensolver,          ONLY: tddfpt_davidson_solver,&
                                              tddfpt_orthogonalize_psi1_psi0,&
                                              tddfpt_orthonormalize_psi1_psi1
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_excitation_analysis,&
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
   USE qs_tddfpt2_utils,                ONLY: ev_shift_operator,&
                                              tddfpt_guess_vectors,&
                                              tddfpt_init_ground_state_mos,&
                                              tddfpt_release_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
   USE xc_pot_saop,                     ONLY: add_saop_pot
   USE xtb_matrices,                    ONLY: build_xtb_ks_matrix
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
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 refactored to be used with Davidson eigensolver [Sergey Chulkov]
!>    * 03.2017 cleaned and refactored [Sergey Chulkov]
!> \note Based on the subroutines tddfpt_env_init(), and tddfpt_env_deallocate().
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt', routineP = moduleN//':'//routineN

      CHARACTER(len=20)                                  :: nstates_str
      INTEGER :: energy_unit, handle, ispin, istate, iter, log_unit, mult, niters, nmo_avail, &
         nmo_occ, nmo_virt, nspins, nstates, nstates_read
      LOGICAL                                            :: do_admm, do_hfx, explicit_xc, &
                                                            is_restarted
      REAL(kind=dp)                                      :: C_hf, conv
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      REAL(kind=dp), DIMENSION(:), POINTER               :: evals_virt_spin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: evals_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dipole_op_mos_occ, evects, S_evects
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_virt
      TYPE(cp_fm_type), POINTER                          :: mos_virt_spin
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_oep, matrix_s
      TYPE(dbcsr_type), POINTER                          :: matrix_ks_oep_spin
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env, kernel_env_admm_aux
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, &
                                                            tddfpt_print_section, tddfpt_section, &
                                                            xc_fun_empty, xc_fun_original, &
                                                            xc_section
      TYPE(stda_env_type), TARGET                        :: stda_kernel
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), ALLOCATABLE, &
         DIMENSION(:)                                    :: gs_mos
      TYPE(tddfpt_kernel_env_type)                       :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(Iannuzzi2005)

      NULLIFY (blacs_env, cell, dft_control, input, matrix_ks, matrix_s, mos)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, dft_control=dft_control, &
                      input=input, matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, scf_env=scf_env)
      tddfpt_control => dft_control%tddfpt2_control

      nspins = dft_control%nspins
      CPASSERT(dft_control%nimages <= 1)

      IF (tddfpt_control%nstates <= 0) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL cp_warn(__LOCATION__, "TDDFPT calculation was requested for "// &
                      TRIM(nstates_str)//" excited states: nothing to do.")
         CALL timestop(handle)
         RETURN
      END IF

      tddfpt_section => section_vals_get_subs_vals(input, "PROPERTIES%TDDFPT")
      tddfpt_print_section => section_vals_get_subs_vals(tddfpt_section, "PRINT")

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

         CALL get_qs_env(qs_env, admm_env=admm_env)
      END IF

      ! reset rks_triplets if UKS is in use
      IF (tddfpt_control%rks_triplets .AND. nspins > 1) THEN
         tddfpt_control%rks_triplets = .FALSE.
         CALL cp_warn(__LOCATION__, "Keyword RKS_TRIPLETS has been ignored for spin-polarised calculations")
      END IF

      ! obtain occupied and virtual (unoccupied) ground-state Kohn-Sham orbitals
      ALLOCATE (gs_mos(nspins))

      ! when the number of unoccupied orbitals is limited and OT has been used
      ! for the ground-state DFT calculation,
      ! compute the missing unoccupied orbitals using OT as well.
      NULLIFY (evals_virt, evals_virt_spin, mos_virt, mos_virt_spin)
      IF (ASSOCIATED(scf_env)) THEN
         IF (scf_env%method == ot_method_nr .AND. tddfpt_control%nlumo > 0) THEN
            ! As OT with ADDED_MOS/=0 is currently not implemented, the following block is equivalent to:
            ! nmo_virt = tddfpt_control%nlumo

            ! number of already computed unoccupied orbitals (added_mos) .
            nmo_virt = HUGE(nmo_virt)
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo_avail, homo=nmo_occ)
               nmo_virt = MIN(nmo_virt, nmo_avail - nmo_occ)
            END DO
            ! number of unoccupied orbitals to compute
            nmo_virt = tddfpt_control%nlumo - nmo_virt

            IF (nmo_virt > 0) THEN
               ALLOCATE (evals_virt(nspins), mos_virt(nspins))
               ! the number of actually computed unoccupied orbitals will be stored as 'nmo_avail'
               CALL make_lumo_gpw(qs_env, scf_env, mos_virt, evals_virt, nmo_virt, nmo_avail)
            END IF
         END IF
      END IF

      ! obtain corrected KS-matrix
      ! We should 'save' the energy values?
      NULLIFY (matrix_ks_oep)
      IF (tddfpt_control%oe_corr /= oe_none) THEN
         CALL cp_warn(__LOCATION__, &
                      "Orbital energy correction potential is an experimental feature. "// &
                      "Use it with extreme care")

         IF (do_hfx) THEN
            CALL cp_abort(__LOCATION__, &
                          "Implementation of orbital energy correction XC-potentials is "// &
                          "currently incompartible with exact-exchange functionals")
         END IF

         CALL dbcsr_allocate_matrix_set(matrix_ks_oep, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(matrix_ks_oep(ispin)%matrix)
            CALL dbcsr_copy(matrix_ks_oep(ispin)%matrix, matrix_ks(ispin)%matrix)
         END DO

         ! KS-matrix without XC-terms
         xc_fun_original => section_vals_get_subs_vals(input, "DFT%XC%XC_FUNCTIONAL")
         CALL section_vals_retain(xc_fun_original)
         NULLIFY (xc_fun_empty)
         CALL section_vals_create(xc_fun_empty, xc_fun_original%section)
         CALL section_vals_set_subs_vals(input, "DFT%XC%XC_FUNCTIONAL", xc_fun_empty)
         CALL section_vals_release(xc_fun_empty)

!deb     CALL get_qs_env(qs_env, energy=energy_init)
         NULLIFY (energy)
!deb     CALL allocate_qs_energy(energy)
!deb     CALL set_qs_env(qs_env, qs_energy=energy)
         IF (dft_control%qs_control%semi_empirical) THEN
            CPABORT("TDDFPT with SE not possible")
         ELSEIF (dft_control%qs_control%dftb) THEN
            CPABORT("TDDFPT with DFTB not possible")
         ELSEIF (dft_control%qs_control%xtb) THEN
            CALL build_xtb_ks_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                     ext_ks_matrix=matrix_ks_oep)
         ELSE
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                           ext_ks_matrix=matrix_ks_oep)
         END IF
!deb     CALL deallocate_qs_energy(energy)
!deb     CALL set_qs_env(qs_env, energy=energy_init)

         IF (tddfpt_control%oe_corr == oe_saop .OR. &
             tddfpt_control%oe_corr == oe_lb .OR. &
             tddfpt_control%oe_corr == oe_gllb) THEN
            CALL add_saop_pot(matrix_ks_oep, qs_env, tddfpt_control%oe_corr)
         ELSE IF (tddfpt_control%oe_corr == oe_shift) THEN
            CALL ev_shift_operator(qs_env, matrix_ks_oep, &
                                   tddfpt_control%ev_shift, tddfpt_control%eos_shift)
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "Unimplemented orbital energy correction potential")
         END IF
         CALL section_vals_set_subs_vals(input, "DFT%XC%XC_FUNCTIONAL", xc_fun_original)
         CALL section_vals_release(xc_fun_original)
      END IF

      DO ispin = 1, nspins
         IF (ASSOCIATED(evals_virt)) THEN
            evals_virt_spin => evals_virt(ispin)%array
         ELSE
            NULLIFY (evals_virt_spin)
         END IF
         IF (ASSOCIATED(mos_virt)) THEN
            mos_virt_spin => mos_virt(ispin)%matrix
         ELSE
            NULLIFY (mos_virt_spin)
         END IF
         IF (ASSOCIATED(matrix_ks_oep)) THEN
            matrix_ks_oep_spin => matrix_ks_oep(ispin)%matrix
         ELSE
            NULLIFY (matrix_ks_oep_spin)
         END IF
         CALL tddfpt_init_ground_state_mos(gs_mos=gs_mos(ispin), mo_set=mos(ispin)%mo_set, &
                                           nlumo=tddfpt_control%nlumo, &
                                           blacs_env=blacs_env, cholesky_method=cholesky_restore, &
                                           matrix_ks=matrix_ks(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                                           mos_virt=mos_virt_spin, evals_virt=evals_virt_spin, &
                                           matrix_ks_oep=matrix_ks_oep_spin)
      END DO

      IF (ASSOCIATED(evals_virt)) THEN
         DO ispin = 1, SIZE(evals_virt)
            IF (ASSOCIATED(evals_virt(ispin)%array)) DEALLOCATE (evals_virt(ispin)%array)
         END DO
         DEALLOCATE (evals_virt)
      END IF

      IF (ASSOCIATED(mos_virt)) THEN
         DO ispin = 1, SIZE(mos_virt)
            IF (ASSOCIATED(mos_virt(ispin)%matrix)) CALL cp_fm_release(mos_virt(ispin)%matrix)
         END DO
         DEALLOCATE (mos_virt)
      END IF

      IF (ASSOCIATED(matrix_ks_oep)) matrix_ks => matrix_ks_oep

      ! components of the dipole operator
      CALL tddfpt_dipole_operator(dipole_op_mos_occ, tddfpt_control, gs_mos, qs_env)

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
         CALL allocate_stda_env(qs_env, stda_kernel)
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

      CALL tddfpt_orthogonalize_psi1_psi0(evects, work_matrices%S_C0_C0T, work_matrices%wfm_ao_mo_occ)
      CALL tddfpt_orthonormalize_psi1_psi1(evects, nstates, S_evects, matrix_s(1)%matrix)

      niters = tddfpt_control%niters
      IF (niters > 0) THEN
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "ITERATION_INFO", extension=".tddfptLog")
         energy_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "DETAILED_ENERGY", extension=".tddfptLog")

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1X,A,/)') 'TDDFPT WAVEFUNCTION OPTIMIZATION'
            WRITE (log_unit, '(5X,A,T15,A,T24,A,T40,A)') "Step", "Time", "Convergence", "Conv. states"
            WRITE (log_unit, '(1X,50("-"))')
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
               WRITE (log_unit, '(1X,10("-"),1X,A,1X,11("-"))') "Restart Davidson iterations"
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
               WRITE (log_unit, '(/,1X,A)') "*** TDDFPT run converged in "//TRIM(nstates_str)//" iteration(s) ***"
            ELSE
               WRITE (log_unit, '(/,1X,A)') "*** TDDFPT run did NOT converge after "// &
                  TRIM(nstates_str)//" iteration(s) ***"
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
      CALL tddfpt_print_excitation_analysis(log_unit, evects, gs_mos, matrix_s(1)%matrix, &
                                            min_amplitude=tddfpt_control%min_excitation_amplitude)

      ! -- clean up all useless stuff
      DO istate = SIZE(evects, 2), 1, -1
         DO ispin = nspins, 1, -1
            CALL cp_fm_release(evects(ispin, istate)%matrix)
            CALL cp_fm_release(S_evects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (evects, S_evects, evals)

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
      IF (do_admm) THEN
            CALL tddfpt_release_kernel_env(kernel_env%admm_kernel)
         END IF
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

END MODULE qs_tddfpt2_methods

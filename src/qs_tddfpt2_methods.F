!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
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
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE exstates_types,                  ONLY: excited_energy_type
   USE header,                          ONLY: tddfpt_header,&
                                              tddfpt_soc_header
   USE hfx_types,                       ONLY: compare_hfx_sections
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
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kernel_methods,               ONLY: create_fxc_kernel,&
                                              create_kernel_env
   USE qs_kernel_types,                 ONLY: full_kernel_env_type,&
                                              kernel_env_type,&
                                              release_kernel_env
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tddfpt2_assign,               ONLY: assign_state
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_eigensolver,          ONLY: tddfpt_davidson_solver,&
                                              tddfpt_orthogonalize_psi1_psi0,&
                                              tddfpt_orthonormalize_psi1_psi1
   USE qs_tddfpt2_forces,               ONLY: tddfpt_forces_main
   USE qs_tddfpt2_fprint,               ONLY: tddfpt_print_forces
   USE qs_tddfpt2_lri_utils,            ONLY: tddfpt2_lri_init
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_excitation_analysis,&
                                              tddfpt_print_nto_analysis,&
                                              tddfpt_print_summary
   USE qs_tddfpt2_restart,              ONLY: tddfpt_read_restart,&
                                              tddfpt_write_newtonx_output,&
                                              tddfpt_write_restart
   USE qs_tddfpt2_soc,                  ONLY: tddfpt_soc
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

   PUBLIC :: tddfpt, tddfpt_energies, tddfpt_input, kernel_info

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

      INTEGER                                            :: handle, ispin, istate, log_unit, mult, &
                                                            my_state, nao, nocc, nspins, &
                                                            nstate_max, nstates, nvirt, old_state
      INTEGER, DIMENSION(maxspins)                       :: nactive
      LOGICAL                                            :: do_admm, do_exck, do_hfx, do_hfxlr, &
                                                            do_soc, lmult_tmp, state_change
      REAL(kind=dp)                                      :: gsmin, gsval, xsval
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals, ostrength
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: my_mos
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ, evects, S_evects
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_oep, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env, kernel_env_admm_aux
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: lri_section, soc_section, &
                                                            tddfpt_print_section, tddfpt_section, &
                                                            xc_section
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

      CALL tddfpt_input(qs_env, do_hfx, do_admm, do_exck, &
                        do_hfxlr, xc_section, tddfpt_print_section, lri_section)

      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      mos=mos, &
                      scf_env=scf_env)
      tddfpt_control => dft_control%tddfpt2_control
      tddfpt_control%do_hfx = do_hfx
      tddfpt_control%do_admm = do_admm
      tddfpt_control%do_hfxlr = do_hfxlr
      tddfpt_control%do_exck = do_exck

      soc_section => section_vals_get_subs_vals(tddfpt_section, "SOC")
      CALL section_vals_get(soc_section, explicit=do_soc)

      IF (do_soc) THEN
         ! start with multiplicity that is not specified in input
         ! so that excited-state gradient is for multiplicity given in input
         lmult_tmp = tddfpt_control%rks_triplets
         tddfpt_control%rks_triplets = .NOT. (tddfpt_control%rks_triplets)
      END IF

      CALL cite_reference(Iannuzzi2005)
      IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         CALL cite_reference(Grimme2013)
         CALL cite_reference(Grimme2016)
      END IF

      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "PROGRAM_BANNER", &
                                      extension=".tddfptLog")
      CALL tddfpt_header(log_unit)
      CALL kernel_info(log_unit, dft_control, tddfpt_control, xc_section)
      ! obtain occupied and virtual (unoccupied) ground-state Kohn-Sham orbitals
      NULLIFY (gs_mos)

      CALL tddfpt_init_mos(qs_env, gs_mos, log_unit)
      CALL cp_print_key_finished_output(log_unit, &
                                        logger, &
                                        tddfpt_print_section, &
                                        "PROGRAM_BANNER")

      ! obtain corrected KS-matrix
      CALL tddfpt_oecorr(qs_env, gs_mos, matrix_ks_oep)

      IF ((tddfpt_control%do_lrigpw) .AND. &
          (tddfpt_control%kernel /= tddfpt_kernel_full)) THEN
         CALL cp_abort(__LOCATION__, "LRI only implemented for full kernel")
      END IF

      IF (ASSOCIATED(matrix_ks_oep)) matrix_ks => matrix_ks_oep

      ! components of the dipole operator
      CALL tddfpt_dipole_operator(dipole_op_mos_occ, &
                                  tddfpt_control, &
                                  gs_mos, &
                                  qs_env)

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
      ALLOCATE (my_mos(nspins))
      DO ispin = 1, nspins
         my_mos(ispin) = gs_mos(ispin)%mos_occ
      END DO
      CALL tddfpt_sub_env_init(sub_env, qs_env, mos_occ=my_mos(:), &
                               kernel=tddfpt_control%kernel)
      DEALLOCATE (my_mos)

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         ! create environment for Full Kernel
         IF (dft_control%qs_control%xtb) THEN
            CPABORT("TDDFPT: xTB only works with sTDA Kernel")
         END IF

         ! allocate pools and work matrices
         nstates = tddfpt_control%nstates
        !! Too many states cal lead to Problems
        !! You should be warned if there are mor states
        !! then occ-virt Combinations!!
         CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=nocc)
         CALL cp_fm_get_info(gs_mos(1)%mos_virt, ncol_global=nvirt)
         nstate_max = nocc*nvirt
         IF (nstates > nstate_max) THEN
            CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
            CPWARN("Experimental: CHANGED NSTATES TO ITS MAXIMUM VALUE!")
            nstates = nstate_max
            tddfpt_control%nstates = nstate_max
         END IF
         CALL tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, do_admm, do_exck, &
                                          qs_env, sub_env)

      END IF

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         kernel_env%full_kernel => full_kernel_env
         kernel_env%admm_kernel => kernel_env_admm_aux
         NULLIFY (kernel_env%stda_kernel)
         ! create full_kernel and admm_kernel within tddfpt_energies
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         ! setup for kernel_stda outside tddfpt_energies
         nactive = 0
         CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)
         DO ispin = 1, SIZE(gs_mos)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, &
                                ncol_global=nactive(ispin))
         END DO
         CALL allocate_stda_env(qs_env, stda_kernel, nao, nactive)
         ! sTDA parameters
         CALL stda_init_param(qs_env, stda_kernel, tddfpt_control%stda_control)
         ! allocate pools and work matrices
         nstates = tddfpt_control%nstates
         CALL stda_create_work_matrices(work_matrices, gs_mos, nstates, &
                                        qs_env, sub_env)
         !
         CALL stda_init_matrices(qs_env, stda_kernel, sub_env, &
                                 work_matrices, tddfpt_control)
         !
         kernel_env%stda_kernel => stda_kernel
         NULLIFY (kernel_env%full_kernel)
         NULLIFY (kernel_env%admm_kernel)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_none) THEN
         ! allocate pools and work matrices
         nstates = tddfpt_control%nstates
         CALL stda_create_work_matrices(work_matrices, gs_mos, nstates, &
                                        qs_env, sub_env)
         NULLIFY (kernel_env%full_kernel)
         NULLIFY (kernel_env%admm_kernel)
         NULLIFY (kernel_env%stda_kernel)
      END IF

      ALLOCATE (evects(nspins, nstates))
      ALLOCATE (evals(nstates))

      ALLOCATE (S_evects(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL fm_pool_create_fm( &
               work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
               S_evects(ispin, istate))
         END DO
      END DO

      IF (.NOT. do_soc) THEN
         ! compute tddfpt excitation energies of multiplicity mult
         CALL tddfpt_energies(qs_env, nstates, work_matrices, &
                              tddfpt_control, logger, tddfpt_print_section, evects, evals, &
                              gs_mos, tddfpt_section, S_evects, matrix_s, kernel_env, matrix_ks, &
                              sub_env, ostrength, dipole_op_mos_occ, mult, xc_section, full_kernel_env, &
                              kernel_env_admm_aux)
      ELSE
         CALL tddfpt_soc_energies(qs_env, nstates, work_matrices, &
                                  tddfpt_control, logger, tddfpt_print_section, &
                                  evects, evals, ostrength, &
                                  gs_mos, tddfpt_section, S_evects, matrix_s, kernel_env, matrix_ks, &
                                  sub_env, dipole_op_mos_occ, lmult_tmp, xc_section, full_kernel_env, &
                                  kernel_env_admm_aux)
      END IF

      !print forces for selected states
      IF (calc_forces) THEN
         CALL tddfpt_print_forces(qs_env, evects, evals, ostrength, &
                                  tddfpt_print_section, gs_mos, &
                                  kernel_env, sub_env, work_matrices)
      END IF

      ! excited state potential energy surface
      IF (qs_env%excited_state) THEN
         IF (sub_env%is_split) THEN
            CALL cp_abort(__LOCATION__, &
                          "Excited state forces not possible when states"// &
                          " are distributed to different CPU pools.")
         END IF
         ! for gradients unshifted KS matrix
         IF (ASSOCIATED(matrix_ks_oep)) CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
         CALL get_qs_env(qs_env, exstate_env=ex_env)
         state_change = .FALSE.
         IF (ex_env%state > 0) THEN
            my_state = ex_env%state
         ELSEIF (ex_env%state < 0) THEN
            ! state following
            my_state = ABS(ex_env%state)
            CALL assign_state(qs_env, matrix_s, evects, mos, ex_env%wfn_history, my_state)
            IF (my_state /= ABS(ex_env%state)) THEN
               state_change = .TRUE.
               old_state = ABS(ex_env%state)
            END IF
            ex_env%state = -my_state
         ELSE
            CALL cp_warn(__LOCATION__, &
                         "Active excited state not assigned. Use the first state.")
            my_state = 1
         END IF
         CPASSERT(my_state > 0)
         IF (my_state > nstates) THEN
            CALL cp_warn(__LOCATION__, &
                         "There were not enough excited states calculated.")
            CPABORT("excited state potential energy surface")
         END IF
         !
         ! energy
         ex_env%evalue = evals(my_state)
         ! excitation vector
         CALL cp_fm_release(ex_env%evect)
         ALLOCATE (ex_env%evect(nspins))
         DO ispin = 1, nspins
            CALL cp_fm_get_info(matrix=evects(ispin, 1), &
                                matrix_struct=matrix_struct)
            CALL cp_fm_create(ex_env%evect(ispin), matrix_struct)
            CALL cp_fm_to_fm(evects(ispin, my_state), ex_env%evect(ispin))
         END DO

         IF (log_unit > 0) THEN
            gsval = ex_env%wfn_history%gsval
            gsmin = ex_env%wfn_history%gsmin
            xsval = ex_env%wfn_history%xsval
            WRITE (log_unit, "(1X,A,T40,F10.6,A,T62,F10.6,A)") "Ground state orbital alignment:", &
               gsmin, "[MinVal]", gsval, "[Average]"
            WRITE (log_unit, "(1X,A,T71,F10.6)") "Excitation vector alignment:", xsval
            IF (state_change) THEN
               WRITE (log_unit, "(1X,A,I5,T60,A14,T76,I5)") &
                  "Target state has been changed from state ", &
                  old_state, " to new state ", my_state
            END IF
            WRITE (log_unit, "(1X,A,I4,A,F12.5,A)") "Calculate properties for state:", &
               my_state, "      with excitation energy ", ex_env%evalue*evolt, " eV"
         END IF

         IF (calc_forces) THEN
            CALL tddfpt_forces_main(qs_env, gs_mos, ex_env, kernel_env, &
                                    sub_env, work_matrices)
         END IF
      END IF

      ! clean up
      CALL cp_fm_release(evects)
      CALL cp_fm_release(S_evects)

      DEALLOCATE (evals, ostrength)

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         IF (do_admm) CALL release_kernel_env(kernel_env%admm_kernel)
         IF (tddfpt_control%do_lrigpw) THEN
            CALL lri_env_release(kernel_env%full_kernel%lri_env)
            DEALLOCATE (kernel_env%full_kernel%lri_env)
            CALL lri_density_release(kernel_env%full_kernel%lri_density)
            DEALLOCATE (kernel_env%full_kernel%lri_density)
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

      CALL cp_fm_release(dipole_op_mos_occ)

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
                                                            explicit_root, expot, exvdw, exwfn, &
                                                            same_hfx
      REAL(kind=dp)                                      :: C_hf
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: hfx_section, hfx_section_gs, input, &
                                                            tddfpt_section, xc_root, xc_sub
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
            !TDDFPT only works if the kernel has the same HF section as the DFT%XC one
            IF (do_hfx) THEN
               hfx_section_gs => section_vals_get_subs_vals(input, "DFT%XC%HF")
               CALL compare_hfx_sections(hfx_section, hfx_section_gs, same_hfx)
               IF (.NOT. same_hfx) THEN
                  CPABORT("TDDFPT Kernel must use the same HF section as DFT%XC or no HF at all.")
               END IF
            END IF

            do_admm = do_hfx .AND. dft_control%do_admm
            IF (do_admm) THEN
               ! 'admm_env%xc_section_primary' and 'admm_env%xc_section_aux' need to be redefined
               CALL cp_abort(__LOCATION__, &
                             "ADMM is not implemented for a TDDFT kernel XC-functional which is different from "// &
                             "the one used for the ground-state calculation. A ground-state 'admm_env' cannot be reused.")
            END IF
            ! SET HFX_KERNEL and/or XC_KERNEL
            IF (exk) THEN
               do_exck = .TRUE.
            ELSE
               do_exck = .FALSE.
            END IF
            IF (exhfxk) THEN
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
                  WRITE (log_unit, "(T2,A,T62,A19)") "KERNEL|", "ADMM Exact Exchange"
                  IF (tddfpt_control%admm_xc_correction) THEN
                     WRITE (log_unit, "(T2,A,T60,A21)") "KERNEL|", "Apply ADMM Kernel XC Correction"
                  END IF
                  IF (tddfpt_control%admm_symm) THEN
                     WRITE (log_unit, "(T2,A,T60,A21)") "KERNEL|", "Symmetric ADMM Kernel"
                  END IF
               ELSE
                  WRITE (log_unit, "(T2,A,T67,A14)") "KERNEL|", "Exact Exchange"
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
!> \brief The energy cylculation has been moved to its own subroutine
!> \param qs_env ...
!> \param nstates ...
!> \param work_matrices ...
!> \param tddfpt_control ...
!> \param logger ...
!> \param tddfpt_print_section ...
!> \param evects ...
!> \param evals ...
!> \param gs_mos ...
!> \param tddfpt_section ...
!> \param S_evects ...
!> \param matrix_s ...
!> \param kernel_env ...
!> \param matrix_ks ...
!> \param sub_env ...
!> \param ostrength ...
!> \param dipole_op_mos_occ ...
!> \param mult ...
!> \param xc_section ...
!> \param full_kernel_env ...
!> \param kernel_env_admm_aux ...
! **************************************************************************************************
   SUBROUTINE tddfpt_energies(qs_env, nstates, work_matrices, &
                              tddfpt_control, logger, tddfpt_print_section, evects, evals, &
                              gs_mos, tddfpt_section, S_evects, matrix_s, kernel_env, matrix_ks, &
                              sub_env, ostrength, dipole_op_mos_occ, mult, xc_section, full_kernel_env, &
                              kernel_env_admm_aux)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: nstates
      TYPE(tddfpt_work_matrices)                         :: work_matrices
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: evects
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: S_evects
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: ostrength
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ
      INTEGER                                            :: mult
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env, kernel_env_admm_aux

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_energies'

      CHARACTER(len=20)                                  :: nstates_str
      INTEGER                                            :: energy_unit, handle, iter, log_unit, &
                                                            niters, nocc, nstate_max, &
                                                            nstates_read, nvirt
      LOGICAL                                            :: do_admm, do_exck, do_soc, explicit, &
                                                            is_restarted
      REAL(kind=dp)                                      :: conv
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_oep
      TYPE(section_vals_type), POINTER                   :: lri_section, namd_print_section, &
                                                            soc_section

      CALL timeset(routineN, handle)

      NULLIFY (admm_env, matrix_ks_oep)
      do_admm = tddfpt_control%do_admm
      IF (do_admm) CALL get_qs_env(qs_env, admm_env=admm_env)

      ! setup for full_kernel and admm_kernel within tddfpt_energies due to dependence on multiplicity
      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN

         CALL tddfpt_construct_ground_state_orb_density( &
            rho_orb_struct=work_matrices%rho_orb_struct_sub, &
            rho_xc_struct=work_matrices%rho_xc_struct_sub, &
            is_rks_triplets=tddfpt_control%rks_triplets, &
            qs_env=qs_env, sub_env=sub_env, &
            wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub)

         IF (do_admm) THEN
            ! Full kernel with ADMM
            IF (tddfpt_control%admm_xc_correction) THEN
               CALL create_kernel_env(kernel_env=full_kernel_env, &
                                      rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                      xc_section=admm_env%xc_section_primary, &
                                      is_rks_triplets=tddfpt_control%rks_triplets, &
                                      sub_env=sub_env)
            ELSE
               CALL create_kernel_env(kernel_env=full_kernel_env, &
                                      rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                      xc_section=xc_section, &
                                      is_rks_triplets=tddfpt_control%rks_triplets, &
                                      sub_env=sub_env)
            END IF

            CALL tddfpt_construct_aux_fit_density( &
               rho_orb_struct=work_matrices%rho_orb_struct_sub, &
               rho_aux_fit_struct=work_matrices%rho_aux_fit_struct_sub, &
               local_rho_set=sub_env%local_rho_set_admm, &
               qs_env=qs_env, sub_env=sub_env, &
               wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub, &
               wfm_rho_aux_fit=work_matrices%rho_ao_aux_fit_fm_sub, &
               wfm_aux_orb=work_matrices%wfm_aux_orb_sub)

            CALL create_kernel_env(kernel_env=kernel_env_admm_aux, &
                                   rho_struct_sub=work_matrices%rho_aux_fit_struct_sub, &
                                   xc_section=admm_env%xc_section_aux, &
                                   is_rks_triplets=tddfpt_control%rks_triplets, &
                                   sub_env=sub_env)
            kernel_env%full_kernel => full_kernel_env
            kernel_env%admm_kernel => kernel_env_admm_aux
         ELSE
            ! Full kernel
            CALL create_kernel_env(kernel_env=full_kernel_env, &
                                   rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                   xc_section=xc_section, &
                                   is_rks_triplets=tddfpt_control%rks_triplets, &
                                   sub_env=sub_env)
            kernel_env%full_kernel => full_kernel_env
            NULLIFY (kernel_env%admm_kernel)
         END IF
         ! Fxc from kernel definition
         do_exck = tddfpt_control%do_exck
         kernel_env%full_kernel%do_exck = do_exck
         ! initilize xc kernel
         IF (do_exck) THEN
            CALL create_fxc_kernel(work_matrices%rho_orb_struct_sub, work_matrices%fxc_rspace_sub, &
                                   xc_section, tddfpt_control%rks_triplets, sub_env, qs_env)
         END IF
      END IF

      ! lri input
      IF (tddfpt_control%do_lrigpw) THEN
         lri_section => section_vals_get_subs_vals(tddfpt_section, "LRIGPW")
         CALL tddfpt2_lri_init(qs_env, kernel_env, lri_section, &
                               tddfpt_print_section)
      END IF

      !! Too many states cal lead to Problems
      !! You should be warned if there are mor states
      !! then occ-virt Combinations!!
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=nocc)
      CALL cp_fm_get_info(gs_mos(1)%mos_virt, ncol_global=nvirt)
      nstate_max = nocc*nvirt
      IF (nstates > nstate_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED NSTATES TO ITS MAXIMUM VALUE!")
         nstates = nstate_max
      END IF

      soc_section => section_vals_get_subs_vals(tddfpt_section, "SOC")
      CALL section_vals_get(soc_section, explicit=do_soc)

      ! reuse Ritz vectors from the previous calculation if available
      IF (tddfpt_control%is_restart .AND. .NOT. do_soc) THEN
         nstates_read = tddfpt_read_restart( &
                        evects=evects, &
                        evals=evals, &
                        gs_mos=gs_mos, &
                        logger=logger, &
                        tddfpt_section=tddfpt_section, &
                        tddfpt_print_section=tddfpt_print_section, &
                        fm_pool_ao_mo_occ=work_matrices%fm_pool_ao_mo_occ, &
                        blacs_env_global=blacs_env)
      ELSE
         nstates_read = 0
      END IF

      is_restarted = nstates_read >= nstates

      ! build the list of missed singly excited states and sort them in ascending order
      ! according to their excitation energies
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                      "GUESS_VECTORS", extension=".tddfptLog")
      CALL tddfpt_guess_vectors(evects=evects, evals=evals, &
                                gs_mos=gs_mos, log_unit=log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, &
                                        tddfpt_print_section, "GUESS_VECTORS")

      CALL tddfpt_orthogonalize_psi1_psi0(evects, work_matrices%S_C0_C0T)
      CALL tddfpt_orthonormalize_psi1_psi1(evects, &
                                           nstates, &
                                           S_evects, &
                                           matrix_s(1)%matrix)

      niters = tddfpt_control%niters
      IF (niters > 0) THEN
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                         "ITERATION_INFO", extension=".tddfptLog")
         energy_unit = cp_print_key_unit_nr(logger, &
                                            tddfpt_print_section, &
                                            "DETAILED_ENERGY", &
                                            extension=".tddfptLog")

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
            conv = tddfpt_davidson_solver( &
                   evects=evects, &
                   evals=evals, &
                   S_evects=S_evects, &
                   gs_mos=gs_mos, &
                   tddfpt_control=tddfpt_control, &
                   matrix_ks=matrix_ks, &
                   qs_env=qs_env, &
                   kernel_env=kernel_env, &
                   sub_env=sub_env, &
                   logger=logger, &
                   iter_unit=log_unit, &
                   energy_unit=energy_unit, &
                   tddfpt_print_section=tddfpt_print_section, &
                   work_matrices=work_matrices)

            ! at this point at least one of the following conditions are met:
            ! a) convergence criteria has been achieved;
            ! b) maximum number of iterations has been reached;
            ! c) Davidson iterations must be restarted due to lack of Krylov vectors or numerical instability

            CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)
            ! terminate the loop if either (a) or (b) is true ...
            IF ((conv <= tddfpt_control%conv &
                 .AND. is_restarted) .OR. iter >= niters) EXIT

            ! ... otherwise restart Davidson iterations
            is_restarted = .TRUE.
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(1X,25("-"),1X,A,1X,25("-"))') "Restart Davidson iterations"
               CALL m_flush(log_unit)
            END IF
         END DO

         ! write TDDFPT restart file at the last iteration if requested to do so
         CALL cp_iterate(logger%iter_info, increment=0, last=.TRUE.)
         CALL tddfpt_write_restart(evects=evects, &
                                   evals=evals, &
                                   gs_mos=gs_mos, &
                                   logger=logger, &
                                   tddfpt_print_section=tddfpt_print_section)

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

         CALL cp_print_key_finished_output(energy_unit, logger, &
                                           tddfpt_print_section, "DETAILED_ENERGY")
         CALL cp_print_key_finished_output(log_unit, logger, &
                                           tddfpt_print_section, "ITERATION_INFO")
      ELSE
         CALL cp_warn(__LOCATION__, &
                      "Skipping TDDFPT wavefunction optimization")
      END IF

      IF (ASSOCIATED(matrix_ks_oep)) THEN
         IF (tddfpt_control%dipole_form == tddfpt_dipole_velocity) THEN
            CALL cp_warn(__LOCATION__, &
                         "Transition dipole moments and oscillator strengths are likely to be incorrect "// &
                         "when computed using an orbital energy correction XC-potential together with "// &
                         "the velocity form of dipole transition integrals")
         END IF
      END IF

      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit(logger)

      namd_print_section => section_vals_get_subs_vals( &
                            tddfpt_print_section, &
                            "NAMD_PRINT")
      CALL section_vals_get(namd_print_section, explicit=explicit)
      IF (explicit) THEN
         CALL tddfpt_write_newtonx_output(evects, &
                                          evals, &
                                          gs_mos, &
                                          logger, &
                                          tddfpt_print_section, &
                                          matrix_s(1)%matrix, &
                                          S_evects, &
                                          sub_env)
      END IF
      ALLOCATE (ostrength(nstates))
      ostrength = 0.0_dp
      CALL tddfpt_print_summary(log_unit, &
                                evects, &
                                evals, &
                                ostrength, &
                                mult, &
                                dipole_op_mos_occ, &
                                tddfpt_control%dipole_form)
      CALL tddfpt_print_excitation_analysis( &
         log_unit, &
         evects, &
         evals, &
         gs_mos, &
         matrix_s(1)%matrix, &
         min_amplitude=tddfpt_control%min_excitation_amplitude)
      CALL tddfpt_print_nto_analysis(qs_env, &
                                     evects, evals, &
                                     ostrength, &
                                     gs_mos, &
                                     matrix_s(1)%matrix, &
                                     tddfpt_print_section)

      IF (tddfpt_control%do_lrigpw) THEN
         CALL lri_print_stat(qs_env, &
                             ltddfpt=.TRUE., &
                             tddfpt_lri_env=kernel_env%full_kernel%lri_env)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_energies

! **************************************************************************************************
!> \brief Perform singlet and triplet computations for subsequent TDDFPT-SOC calculation.
!> \param qs_env  Quickstep environment
!> \param nstates number of requested  exited states
!> \param work_matrices ...
!> \param tddfpt_control ...
!> \param logger ...
!> \param tddfpt_print_section ...
!> \param evects Eigenvector of the requested multiplicity
!> \param evals Eigenvalue of the requested multiplicity
!> \param ostrength Oscillatorstrength
!> \param gs_mos ...
!> \param tddfpt_section ...
!> \param S_evects ...
!> \param matrix_s ...
!> \param kernel_env ...
!> \param matrix_ks ...
!> \param sub_env ...
!> \param dipole_op_mos_occ ...
!> \param lmult_tmp ...
!> \param xc_section ...
!> \param full_kernel_env ...
!> \param kernel_env_admm_aux ...
!> \par History
!>    * 02.2023 created [Jan-Robert Vogt]
!> \note Based on tddfpt2_methods and xas_tdp_utils.
!> \note only the values of one multiplicity will be passed back for force calculations!
! **************************************************************************************************

   SUBROUTINE tddfpt_soc_energies(qs_env, nstates, work_matrices, &
                                  tddfpt_control, logger, tddfpt_print_section, &
                                  evects, evals, ostrength, &
                                  gs_mos, tddfpt_section, S_evects, matrix_s, kernel_env, matrix_ks, &
                                  sub_env, dipole_op_mos_occ, lmult_tmp, xc_section, full_kernel_env, &
                                  kernel_env_admm_aux)

      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      INTEGER, INTENT(in)                                :: nstates
      TYPE(tddfpt_work_matrices)                         :: work_matrices
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: evects
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals, ostrength
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: S_evects
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ
      LOGICAL, INTENT(in)                                :: lmult_tmp
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env, kernel_env_admm_aux

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_soc_energies'

      INTEGER                                            :: handle, ispin, istate, log_unit, mult, &
                                                            nspins
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_mult, ostrength_mult
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: evects_mult

      CALL timeset(routineN, handle)

      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                      "PROGRAM_BANNER", &
                                      extension=".tddfptLog")
      CALL tddfpt_soc_header(log_unit)

      nspins = SIZE(gs_mos)
      ALLOCATE (evects_mult(nspins, nstates))
      ALLOCATE (evals_mult(nstates))

      ! First multiplicity
      IF (lmult_tmp) THEN
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                      TDDFPT SINGLET ENERGIES                                 -", &
               "-------------------------------------------------------------------------------"
         END IF
         mult = 1
      ELSE
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                      TDDFPT TRIPLET ENERGIES                                 -", &
               "-------------------------------------------------------------------------------"
         END IF
         mult = 3
      END IF

      CALL tddfpt_energies(qs_env, nstates, work_matrices, tddfpt_control, logger, &
                           tddfpt_print_section, evects_mult, evals_mult, &
                           gs_mos, tddfpt_section, S_evects, matrix_s, &
                           kernel_env, matrix_ks, sub_env, ostrength_mult, &
                           dipole_op_mos_occ, mult, xc_section, full_kernel_env, &
                           kernel_env_admm_aux)

      ! Clean up in between for full kernel
      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         IF (tddfpt_control%do_admm) CALL release_kernel_env(kernel_env%admm_kernel)
         CALL release_kernel_env(kernel_env%full_kernel)
         CALL tddfpt_release_work_matrices(work_matrices, sub_env)
         CALL tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, tddfpt_control%do_hfx, tddfpt_control%do_admm, &
                                          tddfpt_control%do_exck, qs_env, sub_env)
      END IF

      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL cp_fm_release(S_evects(ispin, istate))
         END DO
      END DO

      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL fm_pool_create_fm( &
               work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
               S_evects(ispin, istate))
         END DO
      END DO

      tddfpt_control%rks_triplets = lmult_tmp

      ! Second multiplicity
      IF (lmult_tmp) THEN
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "                      singlet excitations finished                             ", &
               "                                                                               ", &
               "-------------------------------------------------------------------------------", &
               "-                      TDDFPT TRIPLET ENERGIES                                -", &
               "-------------------------------------------------------------------------------"
         END IF !log_unit
         mult = 3
      ELSE
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "                      triplet excitations finished                             ", &
               "                                                                               ", &
               "-------------------------------------------------------------------------------", &
               "-                      TDDFPT SINGLET ENERGIES                                -", &
               "-------------------------------------------------------------------------------"
         END IF !log_unit
         mult = 1
      END IF

      CALL tddfpt_energies(qs_env, nstates, work_matrices, tddfpt_control, logger, &
                           tddfpt_print_section, evects, evals, &
                           gs_mos, tddfpt_section, S_evects, matrix_s, &
                           kernel_env, matrix_ks, sub_env, ostrength, &
                           dipole_op_mos_occ, mult, xc_section, full_kernel_env, &
                           kernel_env_admm_aux)

      ! Compute perturbative SOC correction
      ! Order should always be singlet triplet in tddfpt_soc
      IF (lmult_tmp) THEN
         CALL tddfpt_soc(qs_env, evals_mult, evals, evects_mult, evects, gs_mos) !mult=singlet
      ELSE
         CALL tddfpt_soc(qs_env, evals, evals_mult, evects, evects_mult, gs_mos) !mult=triplet
      END IF

      ! deallocate the additional multiplicity
      DO ispin = 1, SIZE(evects_mult, 1)
         DO istate = 1, SIZE(evects_mult, 2)
            CALL cp_fm_release(evects_mult(ispin, istate))
         END DO
      END DO
      DEALLOCATE (evects_mult, evals_mult, ostrength_mult)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_soc_energies

END MODULE qs_tddfpt2_methods

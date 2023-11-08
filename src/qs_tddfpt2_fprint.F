!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_fprint
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE efield_utils,                    ONLY: calculate_ecore_efield
   USE exstates_types,                  ONLY: excited_energy_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_environment_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: evolt
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_external_potential,           ONLY: external_c_potential,&
                                              external_e_potential
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              deallocate_qs_force,&
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              sum_qs_force,&
                                              total_qs_force,&
                                              zero_qs_force
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_matrix_w,                     ONLY: qs_energies_compute_matrix_w
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_scf,                          ONLY: scf
   USE qs_tddfpt2_forces,               ONLY: tddfpt_forces_main
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos,&
                                              tddfpt_work_matrices
   USE response_solver,                 ONLY: response_equation,&
                                              response_force,&
                                              response_force_xtb
   USE ri_environment_methods,          ONLY: build_ri_matrices
   USE xtb_matrices,                    ONLY: build_xtb_ks_matrix,&
                                              build_xtb_matrices
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_fprint'

   PUBLIC :: tddfpt_print_forces

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calculate and print forces of selected excited states.
!> \param qs_env             Information on Kinds and Particles
!> \param evects             TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                           SIZE(evects,2) -- number of excited states to print)
!> \param evals              TDDFPT eigenvalues
!> \param ostrength ...
!> \param print_section      ...
!> \param gs_mos             molecular orbitals optimised for the ground state
!> \param kernel_env ...
!> \param sub_env ...
!> \param work_matrices ...
!> \par History
!>    * 10.2022 created [JGH]
! **************************************************************************************************
   SUBROUTINE tddfpt_print_forces(qs_env, evects, evals, ostrength, print_section, &
                                  gs_mos, kernel_env, sub_env, work_matrices)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals, ostrength
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_forces'
      LOGICAL, PARAMETER                                 :: debug_forces = .FALSE.

      INTEGER                                            :: handle, iounit, is, ispin, istate, iw, &
                                                            iwunit, kstate, natom, nspins, nstates
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: alist, natom_of_kind, state_list
      REAL(KIND=dp)                                      :: eener, threshold
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: gs_force, ks_force, td_force
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(section_vals_type), POINTER                   :: force_section, tdlr_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_section, &
                                           "FORCES"), cp_p_file)) THEN

         IF (sub_env%is_split) THEN
            CALL cp_abort(__LOCATION__, "Excited state forces not possible when states"// &
                          " are distributed to different CPU pools.")
         END IF

         nspins = SIZE(evects, 1)
         nstates = SIZE(evects, 2)
         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                     TDDFPT PROPERTIES: Nuclear Forces                       -", &
               "-------------------------------------------------------------------------------"
         END IF
         force_section => section_vals_get_subs_vals(print_section, "FORCES")
         CALL section_vals_val_get(force_section, "THRESHOLD", r_val=threshold)
         tdlr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT%LINRES")
         ALLOCATE (state_list(nstates))
         CALL build_state_list(force_section, state_list)
         ! screen with oscillator strength
         ! Warning: if oscillator strength are not calculated they are set to zero and forces are
         !          only calculated if threshold is also 0
         DO istate = 1, nstates
            IF (ostrength(istate) < threshold) state_list(istate) = 0
         END DO
         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A,T61,E20.8)") " Screening threshold for oscillator strength ", threshold
            ALLOCATE (alist(nstates))
            is = 0
            DO istate = 1, nstates
               IF (state_list(istate) == 1) THEN
                  is = is + 1
                  alist(is) = istate
               END IF
            END DO
            WRITE (iounit, "(1X,A,T71,I10)") " List of states requested for force calculation ", is
            WRITE (iounit, "(16I5)") alist(1:is)
         END IF

         iwunit = cp_print_key_unit_nr(logger, force_section, "", &
                                       extension=".tdfrc", file_status='REPLACE')

         ! prepare force array
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, natom=natom)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
         NULLIFY (td_force, gs_force)
         CALL allocate_qs_force(gs_force, natom_of_kind)
         CALL allocate_qs_force(td_force, natom_of_kind)
         ! ground state forces
         CALL gs_forces(qs_env, gs_force)
         ! Swap force arrays
         CALL get_qs_env(qs_env, force=ks_force)
         CALL set_qs_env(qs_env, force=td_force)

         CALL get_qs_env(qs_env, exstate_env=ex_env, para_env=para_env)
         kstate = ex_env%state
         eener = ex_env%evalue
         ! Start force loop over states
         DO istate = 1, nstates
            IF (state_list(istate) == 1) THEN
               IF (iounit > 0) THEN
                  WRITE (iounit, "(1X,A,I3,T30,F10.5,A,T50,A,T71,F10.7)") "  STATE NR. ", istate, &
                     evals(istate)*evolt, " eV", "Oscillator strength:", ostrength(istate)
               END IF
               IF (iwunit > 0) THEN
                  WRITE (iwunit, "(1X,A,I3,T30,F10.5,A,T50,A,T71,F10.7)") " # STATE NR. ", istate, &
                     evals(istate)*evolt, " eV", "Oscillator strength:", ostrength(istate)
               END IF
               ex_env%state = istate
               ex_env%evalue = evals(istate)
               CALL cp_fm_release(ex_env%evect)
               ALLOCATE (ex_env%evect(nspins))
               DO ispin = 1, nspins
                  CALL cp_fm_get_info(matrix=evects(ispin, 1), matrix_struct=matrix_struct)
                  CALL cp_fm_create(ex_env%evect(ispin), matrix_struct)
                  CALL cp_fm_to_fm(evects(ispin, istate), ex_env%evect(ispin))
               END DO
               ! force array
               CALL zero_qs_force(td_force)
               !
               CALL tddfpt_forces_main(qs_env, gs_mos, ex_env, kernel_env, sub_env, work_matrices)
               !
               IF (debug_forces) THEN
                  iw = iounit
               ELSE
                  iw = -1
               END IF
               CALL response_equation(qs_env, p_env, ex_env%cpmos, iw, tdlr_section)
               !
               CALL get_qs_env(qs_env, dft_control=dft_control)
               IF (dft_control%qs_control%semi_empirical) THEN
                  CPABORT("Not available")
               ELSEIF (dft_control%qs_control%dftb) THEN
                  CPABORT("Not available")
               ELSEIF (dft_control%qs_control%xtb) THEN
                  CALL response_force_xtb(qs_env, p_env, ex_env%matrix_hz, ex_env, debug=debug_forces)
               ELSE
                  CALL response_force(qs_env=qs_env, vh_rspace=ex_env%vh_rspace, &
                                      vxc_rspace=ex_env%vxc_rspace, vtau_rspace=ex_env%vtau_rspace, &
                                      vadmm_rspace=ex_env%vadmm_rspace, matrix_hz=ex_env%matrix_hz, &
                                      matrix_pz=ex_env%matrix_px1, matrix_pz_admm=p_env%p1_admm, &
                                      matrix_wz=p_env%w1, p_env=p_env, ex_env=ex_env, debug=debug_forces)
               END IF
               CALL p_env_release(p_env)
               !
               CALL replicate_qs_force(td_force, para_env)
               CALL sum_qs_force(td_force, gs_force)
               CALL pforce(iwunit, td_force, atomic_kind_set, natom)
               !
            ELSE
               IF (iwunit > 0) THEN
                  WRITE (iwunit, "(1X,A,I3,T30,F10.5,A,T50,A,T71,F10.7)") " # STATE NR. ", istate, &
                     evals(istate)*evolt, " eV", "Oscillator strength:", ostrength(istate)
               END IF
            END IF
         END DO
         CALL set_qs_env(qs_env, force=ks_force)
         CALL deallocate_qs_force(gs_force)
         CALL deallocate_qs_force(td_force)
         DEALLOCATE (state_list)

         ex_env%state = kstate
         ex_env%evalue = eener

         CALL cp_print_key_finished_output(iwunit, logger, force_section, "")

         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A)") &
               "-------------------------------------------------------------------------------"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_print_forces

! **************************************************************************************************
!> \brief   Calculate the Quickstep forces. Adapted from qs_forces
!> \param   qs_env ...
!> \param   gs_force ...
!> \date    11.2022
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE gs_forces(qs_env, gs_force)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: gs_force

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gs_forces'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_w_b, matrix_w_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      ! Swap froce arrays
      CALL get_qs_env(qs_env, force=force)
      CALL zero_qs_force(gs_force)
      CALL set_qs_env(qs_env, force=gs_force)

      ! Check for incompatible options
      CALL get_qs_env(qs_env, dft_control=dft_control)
      CPASSERT(.NOT. dft_control%qs_control%cdft)
      CPASSERT(.NOT. qs_env%run_rtp)
      CPASSERT(.NOT. dft_control%qs_control%mulliken_restraint)
      CPASSERT(.NOT. dft_control%dft_plus_u)
      CPASSERT(.NOT. qs_env%energy_correction)
      IF (ASSOCIATED(qs_env%mp2_env)) THEN
         CPABORT("TDDFPT| MP2 not available")
      END IF

      ! Save current W matrix
      CALL get_qs_env(qs_env=qs_env, matrix_w_kp=matrix_w_b)
      NULLIFY (matrix_w_kp)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      CALL set_ks_env(ks_env, matrix_w_kp=matrix_w_kp)

      ! recalculate energy with forces
      CPASSERT(.NOT. dft_control%qs_control%do_ls_scf)
      CPASSERT(.NOT. dft_control%qs_control%do_almo_scf)
      CALL qs_env_rebuild_pw_env(qs_env)
      CALL qs_energies_init(qs_env, .TRUE.)
      CALL scf(qs_env)
      CALL qs_energies_compute_matrix_w(qs_env, .TRUE.)

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      IF (dft_control%qs_control%semi_empirical) THEN
         CPABORT("TDDFPT| SE not available")
      ELSEIF (dft_control%qs_control%dftb) THEN
         CPABORT("TDDFPT| DFTB not available")
      ELSEIF (dft_control%qs_control%xtb) THEN
         CALL build_xtb_matrices(qs_env=qs_env, para_env=para_env, &
                                 calculate_forces=.TRUE.)
         CALL build_xtb_ks_matrix(qs_env, calculate_forces=.TRUE., just_energy=.FALSE.)
      ELSE
         CALL build_core_hamiltonian_matrix(qs_env=qs_env, calculate_forces=.TRUE.)
         CALL calculate_ecore_self(qs_env)
         CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=.TRUE.)
         CALL calculate_ecore_efield(qs_env, calculate_forces=.TRUE.)
         CALL external_e_potential(qs_env)
         IF (.NOT. dft_control%qs_control%gapw) THEN
            CALL external_c_potential(qs_env, calculate_forces=.TRUE.)
         END IF
         IF (dft_control%qs_control%rigpw) THEN
            CALL get_qs_env(qs_env=qs_env, lri_env=lri_env)
            CALL build_ri_matrices(lri_env, qs_env, calculate_forces=.TRUE.)
         END IF
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.TRUE., just_energy=.FALSE.)
      END IF

      ! replicate forces (get current pointer)
      NULLIFY (gs_force)
      CALL get_qs_env(qs_env=qs_env, force=gs_force)
      CALL replicate_qs_force(gs_force, para_env)
      ! Swap back force array
      CALL set_qs_env(qs_env=qs_env, force=force)

      ! deallocate W Matrix and bring back saved one
      CALL get_qs_env(qs_env=qs_env, matrix_w_kp=matrix_w_kp)
      CALL dbcsr_deallocate_matrix_set(matrix_w_kp)
      NULLIFY (matrix_w_kp)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      CALL set_ks_env(ks_env, matrix_w_kp=matrix_w_b)

      CALL timestop(handle)

   END SUBROUTINE gs_forces

! **************************************************************************************************
!> \brief building a state list
!> \param section input section
!> \param state_list ...
! **************************************************************************************************
   SUBROUTINE build_state_list(section, state_list)

      TYPE(section_vals_type), POINTER                   :: section
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: state_list

      INTEGER                                            :: i, is, k, n_rep, nstate
      INTEGER, DIMENSION(:), POINTER                     :: indexes
      LOGICAL                                            :: explicit

      nstate = SIZE(state_list)

      CALL section_vals_val_get(section, "LIST", explicit=explicit, n_rep_val=n_rep)
      IF (explicit) THEN
         state_list = 0
         DO i = 1, n_rep
            CALL section_vals_val_get(section, "LIST", i_rep_val=i, i_vals=indexes)
            DO is = 1, SIZE(indexes)
               k = indexes(is)
               IF (k <= 0 .OR. k > nstate) THEN
                  CALL cp_warn(__LOCATION__, "State List contains invalid state.")
                  CPABORT("TDDFPT Print Forces: Invalid State")
               END IF
               state_list(k) = 1
            END DO
         END DO
      ELSE
         state_list = 1
      END IF

   END SUBROUTINE build_state_list

! **************************************************************************************************
!> \brief ...
!> \param iwunit ...
!> \param td_force ...
!> \param atomic_kind_set ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE pforce(iwunit, td_force, atomic_kind_set, natom)
      INTEGER, INTENT(IN)                                :: iwunit
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: td_force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER, INTENT(IN)                                :: natom

      INTEGER                                            :: iatom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: force

      ALLOCATE (force(3, natom))
      CALL total_qs_force(force, td_force, atomic_kind_set)

      IF (iwunit > 0) THEN
         WRITE (iwunit, *) natom
         WRITE (iwunit, *)
         DO iatom = 1, natom
            WRITE (iwunit, "(3F24.14)") - force(1:3, iatom)
         END DO
      END IF
      DEALLOCATE (force)

   END SUBROUTINE pforce

END MODULE qs_tddfpt2_fprint

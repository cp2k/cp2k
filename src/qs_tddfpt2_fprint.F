!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_fprint
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE exstates_types,                  ONLY: excited_energy_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              deallocate_qs_force,&
                                              qs_force_type,&
                                              sum_qs_force,&
                                              total_qs_force,&
                                              zero_qs_force
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_tddfpt2_forces,               ONLY: tddfpt_forces_main
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos,&
                                              tddfpt_work_matrices
   USE response_solver,                 ONLY: response_equation,&
                                              response_force,&
                                              response_force_xtb
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
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
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
                                                            iwunit, natom, nkind, nspins, nstates, kstate
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: alist, natom_of_kind, state_list
      REAL(KIND=dp)                                      :: threshold, eener
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: ks_force, td_force
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(section_vals_type), POINTER                   :: force_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_section, &
                                           "FORCES"), cp_p_file)) THEN

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
         CALL get_qs_env(qs_env, force=ks_force, atomic_kind_set=atomic_kind_set, natom=natom)
         nkind = SIZE(atomic_kind_set)
         ALLOCATE (natom_of_kind(nkind))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
         NULLIFY (td_force)
         CALL allocate_qs_force(td_force, natom_of_kind)
         CALL set_qs_env(qs_env, force=td_force)

         CALL get_qs_env(qs_env, exstate_env=ex_env)
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
               IF (ASSOCIATED(ex_env%evect)) THEN
                  DO ispin = 1, SIZE(ex_env%evect)
                     CALL cp_fm_release(ex_env%evect(ispin)%matrix)
                     DEALLOCATE (ex_env%evect(ispin)%matrix)
                  END DO
                  DEALLOCATE (ex_env%evect)
               END IF
               ALLOCATE (ex_env%evect(nspins))
               DO ispin = 1, nspins
                  CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, matrix_struct=matrix_struct)
                  ALLOCATE (ex_env%evect(ispin)%matrix)
                  CALL cp_fm_create(ex_env%evect(ispin)%matrix, matrix_struct)
                  CALL cp_fm_to_fm(evects(ispin, istate)%matrix, ex_env%evect(ispin)%matrix)
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
               CALL response_equation(qs_env, p_env, ex_env%cpmos, iw)
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
         CALL deallocate_qs_force(td_force)
         DEALLOCATE (natom_of_kind, state_list)

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

      WRITE (iwunit, *) natom
      WRITE (iwunit, *)
      DO iatom = 1, natom
         WRITE (iwunit, "(3F24.14)") -force(1:3, iatom)
      END DO
      DEALLOCATE (force)

   END SUBROUTINE pforce

END MODULE qs_tddfpt2_fprint

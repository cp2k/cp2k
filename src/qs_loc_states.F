!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
! **************************************************************************************************
MODULE qs_loc_states
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_p_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE molecular_dipoles,               ONLY: calculate_molecular_dipole
   USE molecular_moments,               ONLY: calculate_molecular_moments
   USE molecular_states,                ONLY: construct_molecular_states
   USE molecule_types,                  ONLY: molecule_type
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_loc_methods,                  ONLY: qs_loc_driver
   USE qs_loc_molecules,                ONLY: wfc_to_molecule
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE wannier_states,                  ONLY: construct_wannier_states
   USE wannier_states_types,            ONLY: wannier_centres_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_states'
   PUBLIC :: get_localization_info

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Performs localization of the orbitals
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param loc_section ...
!> \param mo_local ...
!> \param wf_r ...
!> \param wf_g ...
!> \param particles ...
!> \param coeff ...
!> \param evals ...
!> \param marked_states ...
! **************************************************************************************************
   SUBROUTINE get_localization_info(qs_env, qs_loc_env, loc_section, mo_local, &
                                    wf_r, wf_g, particles, coeff, evals, marked_states)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: loc_section
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_local
      TYPE(pw_p_type)                                    :: wf_r, wf_g
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: coeff
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: evals
      INTEGER, DIMENSION(:, :, :), POINTER               :: marked_states

      CHARACTER(len=*), PARAMETER :: routineN = 'get_localization_info', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, mystate, ns, output_unit
      INTEGER, DIMENSION(:), POINTER                     :: lstates, marked_states_spin
      LOGICAL                                            :: do_homo
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: scenter
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(section_vals_type), POINTER                   :: loc_print_section
      TYPE(wannier_centres_type), DIMENSION(:), POINTER  :: wc

      CALL timeset(routineN, handle)
      NULLIFY (mos, ks_rmpv, dft_control, loc_print_section, marked_states_spin, &
               matrix_s, scenter, wc)
      CALL get_qs_env(qs_env, mos=mos, matrix_ks=ks_rmpv, molecule_set=molecule_set, &
                      dft_control=dft_control, matrix_s=matrix_s)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      loc_print_section => section_vals_get_subs_vals(loc_section, "PRINT")
      do_homo = qs_loc_env%localized_wfn_control%do_homo
      IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                           "WANNIER_STATES"), cp_p_file)) THEN
         CALL get_qs_env(qs_env=qs_env, WannierCentres=wc)
         IF (.NOT. ASSOCIATED(wc)) THEN
            ALLOCATE (wc(dft_control%nspins))
         ENDIF
      ENDIF
      !
      DO ispin = 1, dft_control%nspins
         !
         IF (do_homo) THEN
            qs_loc_env%tag_mo = "HOMO"
         ELSE
            qs_loc_env%tag_mo = "LUMO"
         END IF

         IF (qs_loc_env%do_localize) THEN
            ! Do the Real localization..
            IF (output_unit > 0 .AND. do_homo) WRITE (output_unit, "(/,T2,A,I3)") &
               "LOCALIZATION| Computing localization properties "// &
               "for OCCUPIED ORBITALS. Spin:", ispin
            IF (output_unit > 0 .AND. (.NOT. do_homo)) WRITE (output_unit, "(/,T2,A,I3)") &
               "LOCALIZATION| Computing localization properties "// &
               "for UNOCCUPIED ORBITALS. Spin:", ispin

            scenter => qs_loc_env%localized_wfn_control%centers_set(ispin)%array

            CALL qs_loc_driver(qs_env, qs_loc_env, loc_print_section, &
                               myspin=ispin, ext_mo_coeff=mo_local(ispin)%matrix)

            ! maps wfc to molecules, and compute the molecular dipoles if required
            IF ((BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                                  "MOLECULAR_DIPOLES"), cp_p_file) .OR. &
                 BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                                  "MOLECULAR_MOMENTS"), cp_p_file) .OR. &
                 BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                                  "MOLECULAR_STATES"), cp_p_file))) THEN
               CALL wfc_to_molecule(qs_loc_env, scenter, molecule_set, ispin, dft_control%nspins)
            END IF

            ! Compute the wannier states
            IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                                 "WANNIER_STATES"), cp_p_file)) THEN
               ns = SIZE(qs_loc_env%localized_wfn_control%loc_states, 1)
               IF (.NOT. ASSOCIATED(wc(ispin)%centres)) THEN
                  ALLOCATE (wc(ispin)%WannierHamDiag(ns))
                  ALLOCATE (wc(ispin)%centres(3, ns))
               ENDIF

               wc(ispin)%centres(:, :) = scenter(1 + (ispin - 1)*3:ispin*3, :)
               lstates => qs_loc_env%localized_wfn_control%loc_states(:, ispin)
               CALL construct_wannier_states(mo_local(ispin)%matrix, &
                                             ks_rmpv(ispin)%matrix, qs_env, loc_print_section=loc_print_section, &
                                             WannierCentres=wc(ispin), ns=ns, states=lstates)
            ENDIF
            ! Compute the molecular states
            IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                                 "MOLECULAR_STATES"), cp_p_file)) THEN
               CALL construct_molecular_states( &
                  molecule_set, mo_local(ispin)%matrix, coeff(ispin)%matrix, &
                  evals(ispin)%array, ks_rmpv(ispin)%matrix, matrix_s(1)%matrix, qs_env, wf_r, wf_g, &
                  loc_print_section=loc_print_section, particles=particles, tag=TRIM(qs_loc_env%tag_mo), &
                  marked_states=marked_states_spin, ispin=ispin)
               IF (ASSOCIATED(marked_states_spin)) THEN
                  IF (.NOT. ASSOCIATED(marked_states)) THEN
                     ALLOCATE (marked_states(SIZE(marked_states_spin), dft_control%nspins, 2))
                  END IF
                  mystate = 1
                  IF (qs_loc_env%tag_mo == "LUMO") mystate = 2
                  marked_states(:, ispin, mystate) = marked_states_spin(:)
                  DEALLOCATE (marked_states_spin)
               END IF
            ENDIF
         END IF
      ENDDO

      ! Compute molecular dipoles
      IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                           "MOLECULAR_DIPOLES"), cp_p_file)) THEN
         CALL calculate_molecular_dipole(qs_env, qs_loc_env, loc_print_section, molecule_set)
      END IF

      ! Compute molecular multipole moments
      IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                           "MOLECULAR_MOMENTS"), cp_p_file)) THEN
         CALL calculate_molecular_moments(qs_env, qs_loc_env, mo_local, loc_print_section, molecule_set)
      END IF
      !
      IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                           "WANNIER_STATES"), cp_p_file)) THEN
         CALL set_qs_env(qs_env=qs_env, WannierCentres=wc)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE get_localization_info

END MODULE qs_loc_states

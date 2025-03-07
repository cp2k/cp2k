!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Perform a QUICKSTEP wavefunction optimization (single point)
!> \par History
!>      none
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_energy
   USE almo_scf,                        ONLY: almo_entry_scf
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_external_control,             ONLY: external_control
   USE dm_ls_scf,                       ONLY: ls_scf
   USE energy_corrections,              ONLY: energy_correction
   USE excited_states,                  ONLY: excited_state_energy
   USE input_constants,                 ONLY: smeagol_runtype_emtransport
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE lri_environment_methods,         ONLY: lri_print_stat
   USE mp2,                             ONLY: mp2_main
   USE qs_active_space_methods,         ONLY: active_space_main
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_energy_utils,                 ONLY: qs_energies_properties
   USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_harris_utils,                 ONLY: harris_energy_correction
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_matrix_w,                     ONLY: compute_matrix_w
   USE qs_nonscf,                       ONLY: nonscf
   USE qs_scf,                          ONLY: scf
   USE qs_tddfpt2_smearing_methods,     ONLY: deallocate_fermi_params
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy'

   PUBLIC :: qs_energies

CONTAINS

! **************************************************************************************************
!> \brief   Driver routine for QUICKSTEP single point wavefunction optimization.
!> \param qs_env ...
!> \param consistent_energies ...
!> \param calc_forces ...
!> \date    29.10.2002
!> \par History
!>          - consistent_energies option added (25.08.2005, TdK)
!>          - introduced driver for energy in order to properly decide between
!>            SCF or RTP (fschiff 02.09)
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE qs_energies(qs_env, consistent_energies, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: consistent_energies, calc_forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_energies'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_consistent_energies, &
                                                            do_excited_state, loverlap_deltat, &
                                                            my_calc_forces, run_rtp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: excited_state_section

      CALL timeset(routineN, handle)

      my_calc_forces = .FALSE.
      IF (PRESENT(calc_forces)) my_calc_forces = calc_forces

      do_consistent_energies = .FALSE.
      IF (PRESENT(consistent_energies)) do_consistent_energies = consistent_energies

      CALL qs_env_rebuild_pw_env(qs_env)

      CALL get_qs_env(qs_env=qs_env, run_rtp=run_rtp)
      IF (.NOT. run_rtp) THEN

         NULLIFY (dft_control, energy)
         CALL qs_energies_init(qs_env, my_calc_forces)
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, scf_control=scf_control, energy=energy)

         ! *** check if only overlap matrix is needed for couplings
         loverlap_deltat = .FALSE.
         NULLIFY (excited_state_section)
         excited_state_section => section_vals_get_subs_vals(qs_env%input, "DFT%EXCITED_STATES")
         CALL section_vals_get(excited_state_section, explicit=do_excited_state)
         IF (do_excited_state) THEN
            CALL section_vals_val_get(excited_state_section, "OVERLAP_DELTAT", &
                                      l_val=loverlap_deltat)
         END IF

         ! *** Perform a SCF run ***
         IF (.NOT. loverlap_deltat) THEN
            IF (scf_control%non_selfconsistent .AND. .NOT. scf_control%force_scf_calculation) THEN
               CALL nonscf(qs_env)
            ELSE IF (dft_control%qs_control%do_ls_scf) THEN
               CALL ls_scf(qs_env)
            ELSE IF (dft_control%qs_control%do_almo_scf) THEN
               CALL almo_entry_scf(qs_env, calc_forces=my_calc_forces)
            ELSE
               ! current-induced forces
               IF (dft_control%smeagol_control%smeagol_enabled .AND. &
                   dft_control%smeagol_control%run_type == smeagol_runtype_emtransport) THEN
                  dft_control%smeagol_control%emforces = my_calc_forces
               END IF

               CALL scf(qs_env)
            END IF
         END IF

         IF (do_consistent_energies) THEN
            CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)
         END IF

         IF (.NOT. (dft_control%qs_control%do_ls_scf .OR. dft_control%qs_control%do_almo_scf)) THEN
            ! Compute MP2 energy
            CALL qs_energies_mp2(qs_env, my_calc_forces)

            IF (.NOT. ASSOCIATED(qs_env%mp2_env)) THEN
               ! do not overwrite w matrix computed by SMEAGOL (current-induced forces)
               IF (.NOT. (dft_control%smeagol_control%smeagol_enabled .AND. &
                          dft_control%smeagol_control%run_type == smeagol_runtype_emtransport)) THEN
                  ! if calculate forces, time to compute the w matrix
                  CALL compute_matrix_w(qs_env, my_calc_forces)
               END IF
            END IF
         END IF

         ! Check for energy correction
         IF (qs_env%harris_method) THEN
            CALL harris_energy_correction(qs_env, my_calc_forces)
         END IF

         ! Do active space calculation
         CALL active_space_main(qs_env)

         ! Check for energy correction
         IF (qs_env%energy_correction) THEN
            CALL energy_correction(qs_env, ec_init=.TRUE., calculate_forces=.FALSE.)
         END IF

         IF (.NOT. loverlap_deltat) THEN
            CALL qs_energies_properties(qs_env, calc_forces)

            CALL excited_state_energy(qs_env, calculate_forces=.FALSE.)
         END IF

         IF (dft_control%tddfpt2_control%do_smearing) THEN
            IF (.NOT. ASSOCIATED(dft_control%tddfpt2_control%smeared_occup)) &
               CPABORT("Smearing occupation not associated.")
            CALL deallocate_fermi_params(dft_control%tddfpt2_control%smeared_occup)
         END IF
         IF (dft_control%qs_control%lrigpw) THEN
            CALL lri_print_stat(qs_env)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies

! **************************************************************************************************
!> \brief Enters the mp2 part of cp2k
!> \param qs_env ...
!> \param calc_forces ...
! **************************************************************************************************

   SUBROUTINE qs_energies_mp2(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      LOGICAL                                            :: should_stop

      ! Compute MP2 energy

      IF (ASSOCIATED(qs_env%mp2_env)) THEN

         CALL external_control(should_stop, "MP2", target_time=qs_env%target_time, &
                               start_time=qs_env%start_time)

         CALL mp2_main(qs_env=qs_env, calc_forces=calc_forces)
      END IF

   END SUBROUTINE qs_energies_mp2

END MODULE qs_energy

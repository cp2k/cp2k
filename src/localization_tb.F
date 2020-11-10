!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Localization/Wannier functions for TB
!> \par History
!> \author JHU (03.2019)
! **************************************************************************************************
MODULE localization_tb
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_loc_dipole,                   ONLY: loc_dipole
   USE qs_loc_states,                   ONLY: get_localization_info
   USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                              qs_loc_env_destroy,&
                                              qs_loc_env_new_type
   USE qs_loc_utils,                    ONLY: loc_write_restart,&
                                              qs_loc_control_init,&
                                              qs_loc_init,&
                                              retain_history
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'localization_tb'
   PUBLIC :: wfn_localization_tb

CONTAINS

! **************************************************************************************************
!> \brief wfn localization
!> \param qs_env ...
!> \param tb_type ...
!> \par History
!>      03.2019 initial version
!> \author JHU
!> \note
! **************************************************************************************************
   SUBROUTINE wfn_localization_tb(qs_env, tb_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*)                                   :: tb_type

      CHARACTER(len=*), PARAMETER :: routineN = 'wfn_localization_tb'

      INTEGER                                            :: handle, iounit, ispin, nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: marked_states
      LOGICAL                                            :: do_homo, do_kpoints, explicit, &
                                                            loc_explicit
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: occupied_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: homo_localized, mo_loc_history, &
                                                            occupied_orbs
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env_homo
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, loc_print_section, &
                                                            localize_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
      localize_section => section_vals_get_subs_vals(dft_section, "LOCALIZE")
      loc_print_section => section_vals_get_subs_vals(localize_section, "PRINT")
      CALL section_vals_get(localize_section, explicit=loc_explicit)
      !
      IF (TRIM(tb_type) == "xTB") THEN
         ! okay
      ELSE
         IF (loc_explicit) THEN
            CPWARN("Wfn localization for this TB type not implemented")
            loc_explicit = .FALSE.
         END IF
      END IF

      IF (loc_explicit) THEN
         CALL section_vals_val_get(localize_section, "LIST", explicit=explicit)
         IF (explicit) THEN
            CPWARN("Localization using LIST of states not implemented for TB methods")
         END IF
         CALL section_vals_val_get(localize_section, "ENERGY_RANGE", explicit=explicit)
         IF (explicit) THEN
            CPWARN("Localization using energy range not implemented for TB methods")
         END IF
         CALL section_vals_val_get(localize_section, "LIST_UNOCCUPIED", explicit=explicit)
         IF (explicit) THEN
            CPWARN("Localization of unoccupied states not implemented for TB methods")
         END IF
         ! localize all occupied states
         IF (iounit > 0) THEN
            WRITE (iounit, "(/,T11,A)") " +++++++++++++ Start Localization of Orbitals +++++++++++++"
         END IF
         !
         CALL get_qs_env(qs_env, &
                         dft_control=dft_control, &
                         do_kpoints=do_kpoints, &
                         subsys=subsys, &
                         particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set)
         CALL qs_subsys_get(subsys, particles=particles)

         IF (do_kpoints) THEN
            CPWARN("Localization not implemented for k-point calculations!!")
         ELSEIF (dft_control%restricted) THEN
            IF (iounit > 0) WRITE (iounit, *) &
               " Unclear how we define MOs / localization in the restricted case ... skipping"
         ELSE
            CALL get_qs_env(qs_env, mos=mos)
            nspins = dft_control%nspins
            ALLOCATE (occupied_orbs(nspins))
            ALLOCATE (occupied_evals(nspins))
            ALLOCATE (homo_localized(nspins))
            DO ispin = 1, nspins
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                               eigenvalues=mo_eigenvalues)
               occupied_orbs(ispin)%matrix => mo_coeff
               occupied_evals(ispin)%array => mo_eigenvalues
               CALL cp_fm_create(homo_localized(ispin)%matrix, occupied_orbs(ispin)%matrix%matrix_struct)
               CALL cp_fm_to_fm(occupied_orbs(ispin)%matrix, homo_localized(ispin)%matrix)
            END DO

            CALL get_qs_env(qs_env, mo_loc_history=mo_loc_history)
            do_homo = .TRUE.

            CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
            CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                                   use_data=COMPLEXDATA1D, &
                                   in_space=RECIPROCALSPACE)

            NULLIFY (marked_states, qs_loc_env_homo)
            CALL qs_loc_env_create(qs_loc_env_homo)
            CALL qs_loc_control_init(qs_loc_env_homo, localize_section, do_homo=do_homo)
            CALL qs_loc_init(qs_env, qs_loc_env_homo, localize_section, homo_localized, do_homo, &
                             .FALSE., mo_loc_history=mo_loc_history)
            CALL get_localization_info(qs_env, qs_loc_env_homo, localize_section, homo_localized, &
                                       wf_r, wf_g, particles, occupied_orbs, occupied_evals, marked_states)

            !retain the homo_localized for future use
            IF (qs_loc_env_homo%localized_wfn_control%use_history) THEN
               CALL retain_history(mo_loc_history, homo_localized)
               CALL set_qs_env(qs_env, mo_loc_history=mo_loc_history)
            END IF

            !write restart for localization of occupied orbitals
            CALL loc_write_restart(qs_loc_env_homo, loc_print_section, mos, &
                                   homo_localized, do_homo)
            CALL cp_fm_vect_dealloc(homo_localized)
            DEALLOCATE (occupied_orbs)
            DEALLOCATE (occupied_evals)
            ! Print Total Dipole if the localization has been performed
            IF (qs_loc_env_homo%do_localize) THEN
               CALL loc_dipole(qs_env%input, dft_control, qs_loc_env_homo, logger, qs_env)
            END IF
            CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
            CALL qs_loc_env_destroy(qs_loc_env_homo)
            IF (ASSOCIATED(marked_states)) THEN
               DEALLOCATE (marked_states)
            END IF
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE wfn_localization_tb

! **************************************************************************************************

END MODULE localization_tb

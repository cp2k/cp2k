!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for Resonant Inelastic XRAY Scattering (RIXS) calculations
!> \author BSG (02.2025)
! **************************************************************************************************
MODULE rixs_methods
   USE bibliography,                    ONLY: VazdaCruz2021,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rixs_control_create,&
                                              rixs_control_release,&
                                              rixs_control_type
   USE cp_control_utils,                ONLY: read_rixs_control
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE header,                          ONLY: rixs_header
   USE input_constants,                 ONLY: xas_1s_type,&
                                              xas_2p_type,&
                                              xas_2s_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_tddfpt2_methods,              ONLY: tddfpt
   USE rixs_types,                      ONLY: rixs_env_create,&
                                              rixs_env_release,&
                                              rixs_env_type,&
                                              tddfpt2_valence_type
   USE xas_tdp_methods,                 ONLY: xas_tdp
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rixs_methods'

   PUBLIC :: rixs

CONTAINS

! **************************************************************************************************
!> \brief Driver for RIXS calculations.
!> \param qs_env the inherited qs_environment
!> \author BSG
! **************************************************************************************************

   SUBROUTINE rixs(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rixs'

      INTEGER                                            :: handle, output_unit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: rixs_section, tddfp2_section, &
                                                            xas_tdp_section

      CALL timeset(routineN, handle)

      NULLIFY (rixs_section)
      rixs_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%RIXS")
      output_unit = cp_logger_get_default_io_unit()

      qs_env%do_rixs = .TRUE.

      CALL cite_reference(VazdaCruz2021)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      xas_tdp_section => section_vals_get_subs_vals(rixs_section, "XAS_TDP")
      tddfp2_section => section_vals_get_subs_vals(rixs_section, "TDDFPT")

      CALL rixs_core(rixs_section, qs_env)

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,(T2,A79))") &
            "*******************************************************************************", &
            "!    Normal termination of Resonant Inelastic X-RAY Scattering calculation    !", &
            "*******************************************************************************"
      END IF

      CALL timestop(handle)

   END SUBROUTINE rixs

! **************************************************************************************************
!> \brief Perform RIXS calculation.
!> \param rixs_section ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rixs_core(rixs_section, qs_env)

      TYPE(section_vals_type), POINTER                   :: rixs_section
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rixs_core'

      INTEGER :: ax, c_nstates, current_state_index, handle, iatom, ispin, istate, k, &
         n_states_donor, nactive_max, nao, nex_atoms, nspins, output_unit, td_state, v_nstates
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nactive, nocc
      LOGICAL                                            :: do_sc, do_sg, roks, uks
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: w_i0, w_if
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dip_block, mu_i0
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: mu_if
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER :: core_evect_struct, cRc_struct, dip_0_struct, &
         gs_coeff_struct, i_dip_0_struct, matrix_struct, overlap_struct, Rc_struct
      TYPE(cp_fm_type)                                   :: dip_0
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: core_evects, cRc, i_dip_0, mat_product, &
                                                            overlap_X, Rc, state_gs_coeffs, &
                                                            Y_overlap_X
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: local_gs_coeffs, mo_active
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER         :: valence_evects
      TYPE(cp_fm_type), POINTER                          :: target_ex_coeffs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(donor_state_type), POINTER                    :: current_state
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rixs_control_type), POINTER                   :: rixs_control
      TYPE(rixs_env_type), POINTER                       :: rixs_env
      TYPE(tddfpt2_valence_type), POINTER                :: valence_state
      TYPE(xas_tdp_env_type), POINTER                    :: core_state

      NULLIFY (rixs_control, dft_control, rixs_env)
      NULLIFY (valence_state, core_state)
      NULLIFY (para_env, blacs_env)
      NULLIFY (local_gs_coeffs, mo_active, valence_evects)
      NULLIFY (dipmat, dip_0_struct, i_dip_0_struct, &
               core_evect_struct, gs_coeff_struct)

      output_unit = cp_logger_get_default_io_unit()

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      blacs_env=blacs_env)
      CALL rixs_control_create(rixs_control)
      CALL read_rixs_control(rixs_control, rixs_section, dft_control%qs_control)

      ! create rixs_env
      CALL rixs_env_create(rixs_env)

      ! first, xas_tdp calculation
      CALL xas_tdp(qs_env, rixs_env)

      do_sg = rixs_control%xas_tdp_control%do_singlet
      do_sc = rixs_control%xas_tdp_control%do_spin_cons

      IF (rixs_control%xas_tdp_control%check_only) THEN
         CPWARN("CHECK_ONLY run for XAS_TDP requested, RIXS and TDDFPT will not be performed.")
      ELSE

         ! then, tddfpt calculation
         CALL tddfpt(qs_env, calc_forces=.FALSE., rixs_env=rixs_env)

         IF (output_unit > 0) THEN
            CALL rixs_header(output_unit)
         END IF

         ! timings for rixs only, excluding xas_tdp and tddft calls
         CALL timeset(routineN, handle)

         IF (do_sg) THEN ! singlet
            nspins = 1
         ELSE IF (do_sc) THEN ! spin-conserving
            nspins = 2
         ELSE
            CPABORT("RIXS only implemented for singlet and spin-conserving excitations")
         END IF

         IF (output_unit > 0) THEN
            IF (dft_control%uks) THEN
               uks = .TRUE.
               WRITE (UNIT=output_unit, FMT="(T2,A)") "RIXS| Unrestricted Open-Shell Kohn-Sham"
            ELSE IF (dft_control%roks) THEN
               roks = .TRUE.
               WRITE (UNIT=output_unit, FMT="(T2,A)") "RIXS| Restricted Open-Shell Kohn-Sham"
            ELSE
               WRITE (UNIT=output_unit, FMT="(T2,A)") "RIXS| Restricted Closed-Shell Kohn-Sham"
            END IF
         END IF

         core_state => rixs_env%core_state
         valence_state => rixs_env%valence_state

         ! gs coefficients from tddfpt
         mo_active => valence_state%mos_active
         ! localised gs coefficients from xas_tdp
         local_gs_coeffs => core_state%mo_coeff
         valence_evects => valence_state%evects

         ! find max nactive for open-shell cases
         ALLOCATE (nactive(nspins))
         ALLOCATE (nocc(nspins))
         DO ispin = 1, nspins
            CALL cp_fm_get_info(matrix=local_gs_coeffs(ispin), ncol_global=nocc(ispin))
            CALL cp_fm_get_info(matrix=mo_active(ispin), nrow_global=nao, ncol_global=nactive(ispin))
            ! res tddft
            IF (nocc(ispin) /= nactive(ispin)) THEN
               IF (.NOT. roks) THEN
                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(T2,A)") "RIXS| Reduced Excitation Space TDDFPT requested."
                  END IF
               END IF
            END IF
         END DO
         nactive_max = MAXVAL(nactive)

         nex_atoms = core_state%nex_atoms
         v_nstates = valence_state%nstates
         c_nstates = core_state%nstates

         IF (rixs_control%core_states > 0) THEN
            rixs_control%core_states = MIN(rixs_control%core_states, c_nstates)
         ELSE
            rixs_control%core_states = c_nstates
         END IF
         c_nstates = rixs_control%core_states

         IF (rixs_control%valence_states > 0) THEN
            rixs_control%valence_states = MIN(rixs_control%valence_states, v_nstates)
         ELSE
            rixs_control%valence_states = v_nstates
         END IF
         v_nstates = rixs_control%valence_states

         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,I5,A,I5)") &
               "RIXS| Using ", c_nstates, " core states out of ", core_state%nstates
            WRITE (UNIT=output_unit, FMT="(T2,A,I5,A,I5,/)") &
               "RIXS| Using ", v_nstates, " valence states out of ", valence_state%nstates
         END IF

         dipmat => core_state%dipmat

         ALLOCATE (core_evects(nspins), state_gs_coeffs(nspins))
         ALLOCATE (dip_block(1, nspins), mu_i0(4, c_nstates), mu_if(4, c_nstates, v_nstates), w_i0(c_nstates), w_if(v_nstates))
         w_if(:) = valence_state%evals(1:v_nstates)*evolt
         ALLOCATE (i_dip_0(nspins))
         ALLOCATE (Rc(nspins), cRc(nspins))
         ALLOCATE (overlap_X(nspins), Y_overlap_X(nspins), mat_product(nspins))

         CALL cp_fm_struct_create(core_evect_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=c_nstates)
         CALL cp_fm_struct_create(gs_coeff_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=1)

         ! looping over ex_atoms and ex_kinds is enough as excited atoms have to be unique
         current_state_index = 1
         DO iatom = 1, nex_atoms
            current_state => core_state%donor_states(current_state_index)

            IF (current_state%state_type == xas_1s_type) THEN
               n_states_donor = 1
            ELSE IF (current_state%state_type == xas_2s_type) THEN
               n_states_donor = 1
            ELSE IF (current_state%state_type == xas_2p_type) THEN
               n_states_donor = 3
            END IF

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,A,A,I3,A,A)") &
                  "RIXS| Calculating dipole moment from core-excited state ", &
                  core_state%state_type_char(current_state%state_type), " for atom ", &
                  current_state%at_index, " of kind ", TRIM(current_state%at_symbol)
            END IF

            IF (do_sg) THEN ! singlet
               target_ex_coeffs => current_state%sg_coeffs
               w_i0(:) = current_state%sg_evals(1:c_nstates)*evolt
            ELSE IF (do_sc) THEN ! spin-conserving
               target_ex_coeffs => current_state%sc_coeffs
               w_i0(:) = current_state%sc_evals(1:c_nstates)*evolt
            END IF

            ! reshape sc and sg coeffs (separate spins to columns)
            DO ispin = 1, nspins
               CALL cp_fm_create(core_evects(ispin), core_evect_struct)
               CALL cp_fm_to_fm_submat(msource=target_ex_coeffs, mtarget=core_evects(ispin), s_firstrow=1, &
                                       s_firstcol=(c_nstates*(ispin - 1) + 1), t_firstrow=1, t_firstcol=1, nrow=nao, ncol=c_nstates)
            END DO
            DO ispin = 1, nspins
               CALL cp_fm_create(state_gs_coeffs(ispin), gs_coeff_struct)
               IF (roks) THEN
                  ! store same coeffs for both spins, easier later on
                  CALL cp_fm_to_fm_submat(msource=current_state%gs_coeffs, mtarget=state_gs_coeffs(ispin), s_firstrow=1, &
                                          s_firstcol=1, t_firstrow=1, t_firstcol=1, nrow=nao, ncol=1)
               ELSE
                  CALL cp_fm_to_fm_submat(msource=current_state%gs_coeffs, mtarget=state_gs_coeffs(ispin), s_firstrow=1, &
                                          s_firstcol=ispin, t_firstrow=1, t_firstcol=1, nrow=nao, ncol=1)
               END IF
            END DO

            ! initialise matrices for i->0
            CALL cp_fm_struct_create(dip_0_struct, para_env=para_env, context=blacs_env, &
                                     nrow_global=nao, ncol_global=1)
            CALL cp_fm_create(dip_0, dip_0_struct)
            CALL cp_fm_struct_create(i_dip_0_struct, para_env=para_env, context=blacs_env, &
                                     nrow_global=c_nstates, ncol_global=1)
            DO ispin = 1, nspins
               CALL cp_fm_create(i_dip_0(ispin), i_dip_0_struct)
            END DO

            DO ispin = 1, nspins
               ! initialise cRc (nocc x nactive)
               CALL cp_fm_struct_create(Rc_struct, para_env=para_env, context=blacs_env, &
                                        nrow_global=nao, ncol_global=nactive(ispin))
               CALL cp_fm_struct_create(cRc_struct, para_env=para_env, context=blacs_env, &
                                        nrow_global=n_states_donor, ncol_global=nactive(ispin))
               CALL cp_fm_create(Rc(ispin), Rc_struct)
               CALL cp_fm_create(cRc(ispin), cRc_struct)

               CALL cp_fm_struct_create(overlap_struct, para_env=para_env, context=blacs_env, &
                                        nrow_global=c_nstates, ncol_global=nactive(ispin))

               CALL cp_fm_struct_create(matrix_struct, para_env=para_env, context=blacs_env, &
                                        nrow_global=c_nstates, ncol_global=nocc(ispin))

               ! initialise overlap_X
               CALL cp_fm_create(overlap_X(ispin), Rc_struct) ! nao x nactive

               ! initialise Y_overlap_X
               CALL cp_fm_create(Y_overlap_X(ispin), overlap_struct) ! c_nstates x nactive

               ! initialise mat_product = Y_overlap_X * cRc
               CALL cp_fm_create(mat_product(ispin), matrix_struct) ! c_nstates x nocc

               CALL cp_fm_struct_release(cRc_struct)
               CALL cp_fm_struct_release(Rc_struct)
               CALL cp_fm_struct_release(overlap_struct)
               CALL cp_fm_struct_release(matrix_struct)
            END DO

            mu_i0 = 0.0_dp
            mu_if = 0.0_dp
            ! 0 -> i
            DO ax = 1, 3

               ! i*R*0
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(dipmat(ax)%matrix, state_gs_coeffs(ispin), dip_0, ncol=1)
                  CALL parallel_gemm('T', 'N', c_nstates, 1, nao, 1.0_dp, core_evects(ispin), dip_0, 0.0_dp, i_dip_0(ispin))
               END DO

               DO istate = 1, c_nstates
                  dip_block = 0.0_dp
                  DO ispin = 1, nspins
                     CALL cp_fm_get_submatrix(fm=i_dip_0(ispin), target_m=dip_block, start_row=istate, &
                                              start_col=1, n_rows=1, n_cols=1)
                     mu_i0(ax, istate) = mu_i0(ax, istate) + dip_block(1, 1)
                  END DO ! ispin
                  mu_i0(4, istate) = mu_i0(4, istate) + mu_i0(ax, istate)**2
               END DO ! istate

            END DO ! ax

            ! i -> f
            DO td_state = 1, v_nstates

               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T9,A,I3,A,F10.4)") &
                     "to valence-excited state ", td_state, " with energy ", w_if(td_state)
               END IF

               DO ispin = 1, nspins
                  ! 1) build S*X (nao x nactive)
                  CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, valence_evects(ispin, td_state), overlap_X(ispin), &
                                               ncol=nactive(ispin))
                  ! 2) build Y^T*S*X (c_states x nactive)
                  CALL parallel_gemm('T', 'N', c_nstates, nactive(ispin), nao, 1.0_dp, core_evects(ispin), &
                                     overlap_X(ispin), 0.0_dp, Y_overlap_X(ispin))
               END DO

               DO ax = 1, 3

                  DO ispin = 1, nspins
                     ! 3) build c_i^T*R*c (k x nactive)
                     CALL cp_dbcsr_sm_fm_multiply(dipmat(ax)%matrix, mo_active(ispin), Rc(ispin), ncol=nactive(ispin))
                     CALL parallel_gemm('T', 'N', n_states_donor, nactive(ispin), nao, 1.0_dp, state_gs_coeffs(ispin), &
                                        Rc(ispin), 0.0_dp, cRc(ispin))

                     ! 4) build mat_product (Y^T*S*X*c_i^T*R*c) (c_nstates x k)
                     CALL parallel_gemm('N', 'T', c_nstates, n_states_donor, nactive(ispin), 1.0_dp, Y_overlap_X(ispin), &
                                        cRc(ispin), 0.0_dp, mat_product(ispin))
                  END DO

                  DO istate = 1, c_nstates
                     DO k = 1, n_states_donor
                        dip_block = 0.0_dp
                        DO ispin = 1, nspins
                           CALL cp_fm_get_submatrix(fm=mat_product(ispin), target_m=dip_block, start_row=istate, &
                                                    start_col=k, n_rows=1, n_cols=1)
                           mu_if(ax, istate, td_state) = mu_if(ax, istate, td_state) + dip_block(1, 1)
                        END DO ! ispin
                     END DO ! k
                     mu_if(4, istate, td_state) = mu_if(4, istate, td_state) + mu_if(ax, istate, td_state)**2
                  END DO ! istate (core)

               END DO ! ax

            END DO ! td_state

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,A,/)") "RIXS| Printing spectrum to file"
            END IF
            CALL print_rixs_to_file(current_state, mu_i0, mu_if, w_i0, w_if, rixs_env, rixs_section, rixs_control)

            current_state_index = current_state_index + 1

            ! cleanup
            DO ispin = 1, nspins
               CALL cp_fm_release(core_evects(ispin))
               CALL cp_fm_release(state_gs_coeffs(ispin))
               CALL cp_fm_release(i_dip_0(ispin))
               CALL cp_fm_release(Rc(ispin))
               CALL cp_fm_release(cRc(ispin))
               CALL cp_fm_release(overlap_X(ispin))
               CALL cp_fm_release(Y_overlap_X(ispin))
               CALL cp_fm_release(mat_product(ispin))
            END DO
            CALL cp_fm_struct_release(i_dip_0_struct)
            CALL cp_fm_struct_release(dip_0_struct)
            CALL cp_fm_release(dip_0)

         END DO ! iatom

         NULLIFY (current_state)

         ! cleanup
         CALL cp_fm_struct_release(core_evect_struct)
         CALL cp_fm_struct_release(gs_coeff_struct)

      END IF

      ! more cleanup
      CALL rixs_control_release(rixs_control)
      CALL rixs_env_release(rixs_env)
      NULLIFY (valence_state, core_state)

      CALL timestop(handle)

   END SUBROUTINE rixs_core

!**************************************************************************************************
!> \brief Print RIXS spectrum.
!> \param donor_state ...
!> \param mu_i0 ...
!> \param mu_if ...
!> \param w_i0 ...
!> \param w_if ...
!> \param rixs_env ...
!> \param rixs_section ...
!> \param rixs_control ...
! **************************************************************************************************
   SUBROUTINE print_rixs_to_file(donor_state, mu_i0, mu_if, w_i0, w_if, &
                                 rixs_env, rixs_section, rixs_control)

      TYPE(donor_state_type), POINTER                    :: donor_state
      REAL(dp), DIMENSION(:, :)                          :: mu_i0
      REAL(dp), DIMENSION(:, :, :)                       :: mu_if
      REAL(dp), DIMENSION(:)                             :: w_i0, w_if
      TYPE(rixs_env_type), POINTER                       :: rixs_env
      TYPE(section_vals_type), POINTER                   :: rixs_section
      TYPE(rixs_control_type), POINTER                   :: rixs_control

      INTEGER                                            :: f, i, output_unit, rixs_unit
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      rixs_unit = cp_print_key_unit_nr(logger, rixs_section, "PRINT%SPECTRUM", &
                                       extension=".rixs", file_position="APPEND", &
                                       file_action="WRITE", file_form="FORMATTED")

      output_unit = cp_logger_get_default_io_unit()

      IF (rixs_unit > 0) THEN

         WRITE (rixs_unit, FMT="(A,/,T2,A,A,A,I3,A,A,A/,A)") &
            "=====================================================================================", &
            "Excitation from ground-state (", &
            rixs_env%core_state%state_type_char(donor_state%state_type), " for atom ", &
            donor_state%at_index, " of kind ", TRIM(donor_state%at_symbol), &
            ") to core-excited state i ", &
            "====================================================================================="

         WRITE (rixs_unit, FMT="(T3,A)") &
            "w_0i (eV)            mu^x_0i (a.u.)  mu^y_0i (a.u.)  mu^z_0i (a.u.)  mu^2_0i (a.u.)"
         DO i = 1, rixs_control%core_states
            WRITE (rixs_unit, FMT="(T2,F10.4,T26,E12.5,T42,E12.5,T58,E12.5,T74,E12.5)") &
               w_i0(i), mu_i0(1, i), mu_i0(2, i), mu_i0(3, i), mu_i0(4, i)
         END DO

         WRITE (rixs_unit, FMT="(A,/,T2,A,/,A)") &
            "=====================================================================================", &
            "Emission from core-excited state i to valence-excited state f ", &
            "====================================================================================="

         WRITE (rixs_unit, FMT="(T3,A)") &
            "w_0i (eV) w_if (eV)  mu^x_if (a.u.)  mu^y_if (a.u.)  mu^z_if (a.u.)  mu^2_if (a.u.)"

         DO i = 1, rixs_control%core_states
            DO f = 1, rixs_control%valence_states
               WRITE (rixs_unit, FMT="(T2,F10.4,T14,F8.4,T26,E12.5,T42,E12.5,T58,E12.5,T74,E12.5)") &
                  w_i0(i), w_if(f), mu_if(1, i, f), mu_if(2, i, f), mu_if(3, i, f), mu_if(4, i, f)
            END DO
         END DO

      END IF

      CALL cp_print_key_finished_output(rixs_unit, logger, rixs_section, "PRINT%SPECTRUM")

   END SUBROUTINE print_rixs_to_file

END MODULE rixs_methods

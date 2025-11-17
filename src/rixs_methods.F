!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE header,                          ONLY: rixs_header
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

   PUBLIC :: rixs, rixs_core

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

      INTEGER :: ax, current_state_index, fstate, handle, iatom, ispin, istate, nao, nex_atoms, &
         nocc_max, nspins, nstates, nvirt, output_unit, td_state
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nocc
      LOGICAL                                            :: do_sc, do_sg, roks, uks
      REAL(dp)                                           :: mu_xyz
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: w_i0, w_if
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dip_block, mu_i0
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: mu_if
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: core_evect_struct, dip_0_struct, &
                                                            dip_f_struct, gs_coeff_struct, &
                                                            i_dip_0_struct, i_dip_f_struct
      TYPE(cp_fm_type)                                   :: dip_0
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: core_evects, dip_f, i_dip_0, i_dip_f, &
                                                            state_gs_coeffs
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: local_gs_coeffs, mo_coeffs
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
      NULLIFY (local_gs_coeffs, mo_coeffs, valence_evects)
      NULLIFY (dipmat, dip_0_struct, i_dip_0_struct, dip_f_struct, i_dip_f_struct, &
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
            END IF
         END IF

         core_state => rixs_env%core_state
         valence_state => rixs_env%valence_state

         ! gs coefficients from tddfpt
         mo_coeffs => valence_state%mos_occ
         ! localised gs coefficients from xas_tdp
         local_gs_coeffs => core_state%mo_coeff
         valence_evects => valence_state%evects

         IF (rixs_control%xas_tdp_control%do_loc) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "RIXS| Found localised XAS_TDP orbitals"
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "RIXS| Rotating TDDFPT vectors..."
            END IF
            CALL rotate_vectors(valence_state%evects, local_gs_coeffs, mo_coeffs, matrix_s(1)%matrix, output_unit)
         END IF

         ! find max nocc for open-shell cases
         ALLOCATE (nocc(nspins))
         DO ispin = 1, nspins
            CALL cp_fm_get_info(matrix=valence_state%mos_occ(ispin), nrow_global=nao, ncol_global=nocc(ispin))
         END DO
         nocc_max = MAXVAL(nocc)

         nex_atoms = core_state%nex_atoms
         nstates = valence_state%nstates
         dipmat => core_state%dipmat

         ALLOCATE (core_evects(nspins), state_gs_coeffs(nspins))
         nvirt = core_state%nvirt
         ALLOCATE (dip_block(1, nspins), mu_i0(4, nvirt), mu_if(4, nvirt, nstates), w_i0(nvirt), w_if(nstates))
         mu_i0 = 0.0_dp
         mu_if = 0.0_dp
         w_if(:) = valence_state%evals(:)*evolt
         ALLOCATE (i_dip_0(nspins))
         ALLOCATE (dip_f(nspins), i_dip_f(nspins))

         CALL cp_fm_struct_create(core_evect_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=nvirt)
         CALL cp_fm_struct_create(gs_coeff_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=1)

         ! looping over ex_atoms and ex_kinds is enough as excited atoms have to be unique
         current_state_index = 1
         DO iatom = 1, nex_atoms
            current_state => core_state%donor_states(current_state_index)
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,A,A,A,A,I5)") &
                  "RIXS| Calculating dipole moment from core-excited state ", &
                  core_state%state_type_char(current_state%state_type), " of ", TRIM(current_state%at_symbol), &
                  " with index ", current_state%kind_index
            END IF

            IF (do_sg) THEN ! singlet
               target_ex_coeffs => current_state%sg_coeffs
               w_i0(:) = current_state%sg_evals(:)*evolt
            ELSE IF (do_sc) THEN ! spin-conserving
               target_ex_coeffs => current_state%sc_coeffs
               w_i0(:) = current_state%sc_evals(:)*evolt
            END IF

            ! reshape sc and sg coeffs (separate spins to columns)
            DO ispin = 1, nspins
               CALL cp_fm_create(core_evects(ispin), core_evect_struct)
               CALL cp_fm_to_fm_submat(msource=target_ex_coeffs, mtarget=core_evects(ispin), s_firstrow=1, &
                                       s_firstcol=(nvirt*(ispin - 1) + 1), t_firstrow=1, t_firstcol=1, nrow=nao, ncol=nvirt)
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
                                     nrow_global=nvirt, ncol_global=1)
            DO ispin = 1, nspins
               CALL cp_fm_create(i_dip_0(ispin), i_dip_0_struct)
            END DO

            ! initialise matrices for i->f
            DO ispin = 1, nspins
               CALL cp_fm_struct_create(dip_f_struct, para_env=para_env, context=blacs_env, &
                                        nrow_global=nao, ncol_global=nocc(ispin))
               CALL cp_fm_struct_create(i_dip_f_struct, para_env=para_env, context=blacs_env, &
                                        nrow_global=nvirt, ncol_global=nocc(ispin))
               CALL cp_fm_create(dip_f(ispin), dip_f_struct)
               CALL cp_fm_create(i_dip_f(ispin), i_dip_f_struct)
               CALL cp_fm_struct_release(i_dip_f_struct)
               CALL cp_fm_struct_release(dip_f_struct)
            END DO

            ! 0 -> i
            DO ax = 1, 3

               ! i*R*0
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(dipmat(ax)%matrix, state_gs_coeffs(ispin), dip_0, ncol=1)
                  CALL parallel_gemm('T', 'N', nvirt, 1, nao, 1.0_dp, core_evects(ispin), dip_0, 0.0_dp, i_dip_0(ispin))
               END DO

               DO istate = 1, nvirt
                  dip_block = 0.0_dp
                  DO ispin = 1, nspins
                     CALL cp_fm_get_submatrix(fm=i_dip_0(ispin), target_m=dip_block, start_row=istate, &
                                              start_col=1, n_rows=1, n_cols=1)
                     mu_i0(ax, istate) = dip_block(1, 1)
                  END DO ! ispin
                  mu_xyz = mu_i0(ax, istate)
                  mu_i0(4, istate) = mu_i0(4, istate) + mu_xyz
               END DO ! istate

            END DO ! ax

            ! i -> f
            DO td_state = 1, nstates

               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T9,A,I3,A,F10.4)") &
                     "to valence-excited state ", td_state, " with energy ", w_if(td_state)
               END IF

               DO ax = 1, 3

                  ! core_evects x dipmat x valence_evects (per spin)
                  DO ispin = 1, nspins
                     CALL cp_dbcsr_sm_fm_multiply(dipmat(ax)%matrix, valence_evects(ispin, td_state), dip_f(ispin), &
                                                  ncol=nocc(ispin))
                     CALL parallel_gemm('T', 'N', nvirt, nocc(ispin), nao, 1.0_dp, core_evects(ispin), &
                                        dip_f(ispin), 0.0_dp, i_dip_f(ispin))
                  END DO

                  DO istate = 1, nvirt
                     DO fstate = 1, nocc_max
                        dip_block = 0.0_dp
                        DO ispin = 1, nspins
                           IF (fstate <= nocc(ispin)) THEN
                              CALL cp_fm_get_submatrix(fm=i_dip_f(ispin), target_m=dip_block, start_row=istate, &
                                                       start_col=fstate, n_rows=1, n_cols=1)
                              mu_if(ax, istate, td_state) = mu_if(ax, istate, td_state) + dip_block(1, 1)
                           END IF
                        END DO ! ispin
                     END DO ! fstate (tddft)
                     mu_xyz = mu_if(ax, istate, td_state)
                     mu_if(4, istate, td_state) = mu_if(4, istate, td_state) + mu_xyz
                  END DO ! istate (core)

               END DO ! ax

            END DO ! td_state

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,A,/)") "RIXS| Printing spectrum to file"
            END IF
            CALL print_rixs_to_file(current_state, mu_i0, mu_if, w_i0, w_if, rixs_env, rixs_section)

            current_state_index = current_state_index + 1

            ! cleanup
            DO ispin = 1, nspins
               CALL cp_fm_release(core_evects(ispin))
               CALL cp_fm_release(state_gs_coeffs(ispin))
               CALL cp_fm_release(i_dip_0(ispin))
               CALL cp_fm_release(i_dip_f(ispin))
               CALL cp_fm_release(dip_f(ispin))
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

! **************************************************************************************************
!> \brief Rotate vectors. Returns rotated mo_occ and evects.
!> \param evects ...
!> \param mo_ref ...
!> \param mo_occ ...
!> \param overlap_matrix ...
!> \param unit_nr ...
! **************************************************************************************************

   SUBROUTINE rotate_vectors(evects, mo_ref, mo_occ, overlap_matrix, unit_nr)
      TYPE(cp_fm_type), DIMENSION(:, :)                  :: evects
      TYPE(cp_fm_type), DIMENSION(:)                     :: mo_ref, mo_occ
      TYPE(dbcsr_type), POINTER                          :: overlap_matrix
      INTEGER                                            :: unit_nr

      INTEGER                                            :: ispin, istate, nactive, ncol, nrow, &
                                                            nspins, nstates
      REAL(kind=dp)                                      :: diff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: emat_struct
      TYPE(cp_fm_type)                                   :: emat, rotated_mo_coeffs, smo
      TYPE(cp_fm_type), POINTER                          :: current_evect
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (emat_struct, para_env, blacs_env, current_evect)

      nspins = SIZE(evects, DIM=1)
      DO ispin = 1, nspins

         CALL cp_fm_get_info(matrix=mo_occ(ispin), nrow_global=nrow, ncol_global=ncol, &
                             para_env=para_env, context=blacs_env)
         CALL cp_fm_create(smo, mo_occ(ispin)%matrix_struct)

         ! rotate mo_occ
         ! smo = matrix_s x mo_occ
         CALL cp_dbcsr_sm_fm_multiply(overlap_matrix, mo_occ(ispin), smo, ncol, alpha=1.0_dp, beta=0.0_dp)
         CALL cp_fm_struct_create(emat_struct, nrow_global=ncol, ncol_global=ncol, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(emat, emat_struct)
         ! emat = mo_ref^T x smo
         CALL parallel_gemm('T', 'N', ncol, ncol, nrow, 1.0_dp, mo_ref(ispin), smo, 0.0_dp, emat)
         CALL cp_fm_create(rotated_mo_coeffs, mo_occ(ispin)%matrix_struct)
         ! rotated_mo_coeffs = cpmos x emat
         CALL parallel_gemm('N', 'N', nrow, ncol, ncol, 1.0_dp, mo_occ(ispin), emat, 0.0_dp, rotated_mo_coeffs)

         diff = MAXVAL(ABS(rotated_mo_coeffs%local_data - mo_occ(ispin)%local_data))
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(T9,A,I2,A,F10.6,/)") "For spin ", ispin, ": Max difference between orbitals = ", diff
         END IF

         CALL cp_fm_to_fm(rotated_mo_coeffs, mo_occ(ispin))

         ! rotation of transition vectors
         CALL cp_fm_get_info(matrix=evects(ispin, 1), ncol_global=nactive)
         IF (nactive /= ncol) THEN
            CALL cp_warn(__LOCATION__, "RIXS with reduced occupied state TDDFPT not possible")
            CPABORT("rotate_vectors in rixs_method")
         END IF

         nstates = SIZE(evects, DIM=2)
         DO istate = 1, nstates
            ASSOCIATE (current_evect => evects(ispin, istate))
               CALL parallel_gemm('N', 'N', nrow, ncol, ncol, 1.0_dp, current_evect, emat, 0.0_dp, smo)
               CALL cp_fm_to_fm(smo, current_evect)
            END ASSOCIATE
         END DO

         CALL cp_fm_struct_release(emat_struct)
         CALL cp_fm_release(smo)
         CALL cp_fm_release(emat)
         CALL cp_fm_release(rotated_mo_coeffs)

      END DO ! ispin

   END SUBROUTINE rotate_vectors

!**************************************************************************************************
!> \brief Print RIXS spectrum.
!> \param donor_state ...
!> \param mu_i0 ...
!> \param mu_if ...
!> \param w_i0 ...
!> \param w_if ...
!> \param rixs_env ...
!> \param rixs_section ...
! **************************************************************************************************
   SUBROUTINE print_rixs_to_file(donor_state, mu_i0, mu_if, w_i0, w_if, &
                                 rixs_env, rixs_section)

      TYPE(donor_state_type), POINTER                    :: donor_state
      REAL(dp), DIMENSION(:, :)                          :: mu_i0
      REAL(dp), DIMENSION(:, :, :)                       :: mu_if
      REAL(dp), DIMENSION(:)                             :: w_i0, w_if
      TYPE(rixs_env_type), POINTER                       :: rixs_env
      TYPE(section_vals_type), POINTER                   :: rixs_section

      INTEGER                                            :: f, i, output_unit, rixs_unit
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      rixs_unit = cp_print_key_unit_nr(logger, rixs_section, "PRINT%SPECTRUM", &
                                       extension=".rixs", file_position="APPEND", &
                                       file_action="WRITE", file_form="FORMATTED")

      output_unit = cp_logger_get_default_io_unit()

      IF (rixs_unit > 0) THEN

         WRITE (rixs_unit, FMT="(A,/,T2,A,A,A,A,A,I5,A/,A)") &
            "====================================================================================", &
            "Excitation from ground-state (", &
            rixs_env%core_state%state_type_char(donor_state%state_type), " of kind ", &
            TRIM(donor_state%at_symbol), " with index ", donor_state%kind_index, &
            ") to core-excited state i ", &
            "===================================================================================="

         WRITE (rixs_unit, FMT="(T3,A)") &
            "w_0i (eV)            mu^x_0i (a.u.)  mu^y_0i (a.u.)  mu^z_0i (a.u.)  mu^2_0i (a.u.)"
         DO i = 1, SIZE(mu_i0, DIM=2)
            WRITE (rixs_unit, FMT="(T2,F10.4,T26,E12.5,T42,E12.5,T58,E12.5,T74,E12.5)") &
               w_i0(i), mu_i0(1, i), mu_i0(2, i), mu_i0(3, i), mu_i0(4, i)
         END DO

         WRITE (rixs_unit, FMT="(A,/,T2,A,/,A)") &
            "====================================================================================", &
            "Emission from core-excited state i to valence-excited state f ", &
            "===================================================================================="

         WRITE (rixs_unit, FMT="(T3,A)") &
            "w_0i (eV) w_if (eV)  mu^x_if (a.u.)  mu^y_if (a.u.)  mu^z_if (a.u.)  mu^2_if (a.u.)"

         DO i = 1, SIZE(mu_if, DIM=2)
            DO f = 1, SIZE(mu_if, DIM=3)
               WRITE (rixs_unit, FMT="(T2,F10.4,T14,F8.4,T26,E12.5,T42,E12.5,T58,E12.5,T74,E12.5)") &
                  w_i0(i), w_if(f), mu_if(1, i, f), mu_if(2, i, f), mu_if(3, i, f), mu_if(4, i, f)
            END DO
         END DO

      END IF

      CALL cp_print_key_finished_output(rixs_unit, logger, rixs_section, "PRINT%SPECTRUM")

   END SUBROUTINE print_rixs_to_file

END MODULE rixs_methods

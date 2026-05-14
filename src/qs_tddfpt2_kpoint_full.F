!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_kpoint_full
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_constants,                 ONLY: tddfpt_kernel_full
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_p_type,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathlib,                         ONLY: diag
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kernel_methods,               ONLY: create_kernel_env
   USE qs_kernel_types,                 ONLY: full_kernel_env_type,&
                                              release_kernel_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_tddfpt2_kpoint_utils,         ONLY: &
        tddfpt_kpoint_apply_energy_diff, tddfpt_kpoint_create_matrix_set, &
        tddfpt_kpoint_hmat_from_rs, tddfpt_kpoint_hmat_to_response, &
        tddfpt_kpoint_init_transitions, tddfpt_kpoint_release_matrix_set, &
        tddfpt_kpoint_response_add_scaled_sigma, tddfpt_kpoint_response_copy_scaled_coeff, &
        tddfpt_kpoint_response_create, tddfpt_kpoint_response_release, &
        tddfpt_kpoint_response_set_unit, tddfpt_kpoint_response_type, &
        tddfpt_kpoint_response_update_rho, tddfpt_kpoint_response_zero, &
        tddfpt_kpoint_transition_offset, tddfpt_kpoint_transition_type
   USE qs_tddfpt2_operators,            ONLY: tddfpt_apply_xc
   USE string_utilities,                ONLY: uppercase
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_kpoint_full'
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3

   PUBLIC :: tddfpt_kpoint_full_reference

CONTAINS

! **************************************************************************************************
!> \brief Diagonalize the local/semi-local TDDFPT kernel in the full k-point transition basis.
!> \param qs_env ...
!> \param logger ...
!> \param tddfpt_control ...
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_full_reference(qs_env, logger, tddfpt_control, xc_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_kpoint_full_reference'

      INTEGER                                            :: handle, ideriv, ikernel, img, ispin, &
                                                            istate, itrans, log_unit, nimages, &
                                                            nresponse, nstates, ntrans_kpoint, &
                                                            ntrans_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: imag_index, inds
      INTEGER, DIMENSION(2)                              :: homo_spin, nao_spin, nmo_spin, nvirt_spin
      LOGICAL                                            :: debug_do_hartree, debug_do_xc
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: active_imag
      REAL(kind=dp) :: checksum, coeff_abs, debug_density_scale, debug_density_weight_power, &
         debug_sigma_scale, debug_sigma_weight_power, max_coeff_abs, symmetry_error, ximag, xreal
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals, oscillator_strength, sorted_evals
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvectors, kernel_matrix, &
                                                            state_dipole_im, state_dipole_re, &
                                                            transition_dipole_im, &
                                                            transition_dipole_re
      REAL(kind=dp), DIMENSION(nderivs)                  :: transition_dipole_abs
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: hmat_ext, pmat_ext
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: hmat_rs, matrix_s_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(full_kernel_env_type), POINTER                :: kernel_env_ptr
      TYPE(full_kernel_env_type), TARGET                 :: kernel_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp
      TYPE(pw_c1d_gs_type)                               :: work_v_gspace
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: kernel_rspace, work_v_xc, work_v_xc_tau
      TYPE(pw_r3d_rs_type), POINTER                      :: hartree_rspace
      TYPE(qs_rho_type), POINTER                         :: rho0, rho1
      TYPE(tddfpt_kpoint_response_type)                  :: column_response, kernel_response
      TYPE(tddfpt_kpoint_transition_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: transitions

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, hartree_rspace, hmat_rs, kernel_env_ptr, kpoints, matrix_s_kp, &
               poisson_env, pw_env, rho0, rho1, sab_kp)
      NULLIFY (kernel_rspace, work_v_xc, work_v_xc_tau)

      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints, matrix_s_kp=matrix_s_kp, &
                      pw_env=pw_env, rho=rho0)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(matrix_s_kp))
      CPASSERT(ASSOCIATED(pw_env))
      CPASSERT(ASSOCIATED(rho0))
      CPASSERT(tddfpt_control%kernel == tddfpt_kernel_full)
      CALL read_debug_options(debug_do_xc, debug_do_hartree, debug_density_weight_power, &
                              debug_sigma_weight_power, debug_density_scale, debug_sigma_scale)

      CALL get_kpoint_info(kpoints, sab_nl=sab_kp)
      CPASSERT(ASSOCIATED(sab_kp))
      nimages = dft_control%nimages

      CALL tddfpt_kpoint_init_transitions(qs_env, transitions, ntrans_kpoint, &
                                          nmo_spin, homo_spin, nvirt_spin, nao_spin)
      MARK_USED(ntrans_kpoint)
      MARK_USED(nmo_spin)
      MARK_USED(homo_spin)
      MARK_USED(nvirt_spin)
      MARK_USED(nao_spin)
      ntrans_total = SIZE(transitions)
      ALLOCATE (active_imag(ntrans_total), imag_index(ntrans_total))
      CALL init_response_layout(kpoints, transitions, active_imag, imag_index, nresponse)
      nstates = MIN(tddfpt_control%nstates, nresponse)
      IF (tddfpt_control%nstates > nresponse) THEN
         CPWARN("Requested more TDDFPT states than k-point occupied-virtual transitions")
      END IF
      ALLOCATE (transition_dipole_re(ntrans_total, nderivs), &
                transition_dipole_im(ntrans_total, nderivs))
      CALL build_transition_dipoles(qs_env, tddfpt_control, transitions, ntrans_kpoint, &
                                    nmo_spin, homo_spin, nvirt_spin, nao_spin, &
                                    transition_dipole_re, transition_dipole_im)

      ALLOCATE (rho1)
      CALL qs_rho_create(rho1)
      CALL create_kernel_env(kernel_env=kernel_env, xc_section=xc_section, &
                             is_rks_triplets=tddfpt_control%rks_triplets, &
                             rho_struct_sub=rho0, qs_env=qs_env)

      CALL tddfpt_kpoint_response_create(column_response, ntrans_total)
      CALL tddfpt_kpoint_response_create(kernel_response, ntrans_total)
      CALL tddfpt_kpoint_create_matrix_set(qs_env, pmat_ext)
      CALL tddfpt_kpoint_create_matrix_set(qs_env, hmat_ext)
      CALL dbcsr_allocate_matrix_set(hmat_rs, dft_control%nspins, nimages)
      DO ispin = 1, dft_control%nspins
         DO img = 1, nimages
            ALLOCATE (hmat_rs(ispin, img)%matrix)
            CALL dbcsr_create(hmat_rs(ispin, img)%matrix, name="TDDFPT K-POINT RESPONSE MATRIX", &
                              template=matrix_s_kp(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(hmat_rs(ispin, img)%matrix, sab_kp)
            CALL dbcsr_set(hmat_rs(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
      CPASSERT(ASSOCIATED(auxbas_pw_pool))
      CPASSERT(ASSOCIATED(poisson_env))
      ALLOCATE (kernel_rspace(dft_control%nspins), work_v_xc(dft_control%nspins), &
                work_v_xc_tau(dft_control%nspins))
      DO ispin = 1, dft_control%nspins
         CALL auxbas_pw_pool%create_pw(kernel_rspace(ispin))
         CALL auxbas_pw_pool%create_pw(work_v_xc(ispin))
         CALL auxbas_pw_pool%create_pw(work_v_xc_tau(ispin))
      END DO
      ALLOCATE (hartree_rspace)
      CALL auxbas_pw_pool%create_pw(hartree_rspace)
      CALL auxbas_pw_pool%create_pw(work_v_gspace)

      ALLOCATE (kernel_matrix(nresponse, nresponse), eigenvectors(nresponse, nresponse), &
                evals(nresponse), sorted_evals(nresponse), inds(nresponse))
      kernel_matrix = 0.0_dp

      DO itrans = 1, ntrans_total
         CALL tddfpt_kpoint_response_set_unit(column_response, itrans)
         CALL tddfpt_kpoint_apply_energy_diff(column_response, transitions)
         CALL apply_local_kernel(qs_env, column_response, kernel_response, transitions, kernel_env, &
                                 rho1, pmat_ext, hmat_rs, hmat_ext, kernel_rspace, work_v_xc, &
                                 work_v_xc_tau, hartree_rspace, work_v_gspace, poisson_env, &
                                 tddfpt_control%rks_triplets, debug_do_xc, debug_do_hartree, &
                                 debug_density_weight_power, debug_sigma_weight_power, &
                                 debug_density_scale, debug_sigma_scale)
         kernel_matrix(1:ntrans_total, itrans) = column_response%sigma_re(:)
         CALL set_imag_rows(kernel_matrix, itrans, column_response%sigma_im, active_imag, &
                            imag_index)

         IF (active_imag(itrans)) THEN
            CALL tddfpt_kpoint_response_set_unit(column_response, itrans, imaginary=.TRUE.)
            CALL tddfpt_kpoint_apply_energy_diff(column_response, transitions)
            CALL apply_local_kernel(qs_env, column_response, kernel_response, transitions, &
                                    kernel_env, rho1, pmat_ext, hmat_rs, hmat_ext, kernel_rspace, &
                                    work_v_xc, work_v_xc_tau, hartree_rspace, work_v_gspace, &
                                    poisson_env, tddfpt_control%rks_triplets, debug_do_xc, &
                                    debug_do_hartree, debug_density_weight_power, &
                                    debug_sigma_weight_power, debug_density_scale, debug_sigma_scale)
            kernel_matrix(1:ntrans_total, imag_index(itrans)) = column_response%sigma_re(:)
            CALL set_imag_rows(kernel_matrix, imag_index(itrans), column_response%sigma_im, &
                               active_imag, imag_index)
         END IF
      END DO

      symmetry_error = MAXVAL(ABS(kernel_matrix(:, :) - TRANSPOSE(kernel_matrix(:, :))))
      IF (symmetry_error > 1.0E-6_dp) THEN
         CPWARN("K-point TDDFPT full-kernel matrix is not symmetric within the S(k) metric")
      END IF
      kernel_matrix(:, :) = 0.5_dp*(kernel_matrix(:, :) + TRANSPOSE(kernel_matrix(:, :)))
      CALL diag(nresponse, kernel_matrix, evals, eigenvectors)
      sorted_evals(:) = evals(:)
      CALL sort(sorted_evals, nresponse, inds)

      ALLOCATE (state_dipole_re(nstates, nderivs), state_dipole_im(nstates, nderivs), &
                oscillator_strength(nstates))
      state_dipole_re = 0.0_dp
      state_dipole_im = 0.0_dp
      oscillator_strength = 0.0_dp
      DO istate = 1, nstates
         DO itrans = 1, ntrans_total
            xreal = eigenvectors(itrans, inds(istate))
            ximag = 0.0_dp
            IF (active_imag(itrans)) ximag = eigenvectors(imag_index(itrans), inds(istate))
            DO ideriv = 1, nderivs
               state_dipole_re(istate, ideriv) = state_dipole_re(istate, ideriv) + &
                                                 xreal*transition_dipole_re(itrans, ideriv) - &
                                                 ximag*transition_dipole_im(itrans, ideriv)
               state_dipole_im(istate, ideriv) = state_dipole_im(istate, ideriv) + &
                                                 xreal*transition_dipole_im(itrans, ideriv) + &
                                                 ximag*transition_dipole_re(itrans, ideriv)
            END DO
         END DO
         oscillator_strength(istate) = 2.0_dp/3.0_dp*MAX(0.0_dp, sorted_evals(istate))* &
                                       SUM(state_dipole_re(istate, :)**2 + &
                                           state_dipole_im(istate, :)**2)
      END DO

      checksum = SQRT(SUM(sorted_evals(1:nstates)**2))
      log_unit = cp_logger_get_default_io_unit(logger)
      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                 TDDFPT K-point Local/Semilocal Full Kernel                  -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, "(1X,A)") &
            "Explicit full transition-space reference solver; HFX, SOC, smearing and metals are guarded."
         WRITE (log_unit, "(1X,A,ES12.4E2)") "Metric symmetry error:", symmetry_error
         WRITE (log_unit, '(/,T10,A,T19,A,T37,A,T69,A,T86,A)') "State", "Excitation", &
            "Transition dipole (a.u.)", "Oscillator", "Dominant transition"
         WRITE (log_unit, '(T10,A,T19,A,T37,A,T49,A,T61,A,T67,A,T86,A,T95,A,T104,A,T113,A)') &
            "number", "energy (eV)", "x", "y", "z", "strength (a.u.)", "k", "spin", "occ", &
            "virt"
         WRITE (log_unit, '(T10,112("-"))')
      END IF

      DO istate = 1, nstates
         ikernel = 1
         max_coeff_abs = -1.0_dp
         DO itrans = 1, ntrans_total
            coeff_abs = ABS(eigenvectors(itrans, inds(istate)))
            IF (active_imag(itrans)) coeff_abs = coeff_abs + &
                                                 ABS(eigenvectors(imag_index(itrans), inds(istate)))
            IF (coeff_abs > max_coeff_abs) THEN
               max_coeff_abs = coeff_abs
               ikernel = itrans
            END IF
         END DO
         transition_dipole_abs(:) = SQRT(state_dipole_re(istate, :)**2 + &
                                         state_dipole_im(istate, :)**2)
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(1X,A,T9,I7,T19,F11.5,T31,3(1X,ES11.4E2),T69,ES12.5E2,'// &
                   'T86,I7,T95,I7,T104,I7,T113,I7)') "TDDFPT|", istate, &
               sorted_evals(istate)*evolt, transition_dipole_abs, oscillator_strength(istate), &
               transitions(ikernel)%ikp, transitions(ikernel)%ispin, transitions(ikernel)%iocc, &
               transitions(ikernel)%ivirt
         END IF
      END DO

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,E14.6)') 'TDDFPT : CheckSum  =', checksum
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      CALL auxbas_pw_pool%give_back_pw(work_v_gspace)
      CALL auxbas_pw_pool%give_back_pw(hartree_rspace)
      DEALLOCATE (hartree_rspace)
      DO ispin = 1, dft_control%nspins
         CALL auxbas_pw_pool%give_back_pw(work_v_xc_tau(ispin))
         CALL auxbas_pw_pool%give_back_pw(work_v_xc(ispin))
         CALL auxbas_pw_pool%give_back_pw(kernel_rspace(ispin))
      END DO
      DEALLOCATE (kernel_rspace, work_v_xc, work_v_xc_tau)
      CALL dbcsr_deallocate_matrix_set(hmat_rs)
      CALL tddfpt_kpoint_release_matrix_set(hmat_ext)
      CALL tddfpt_kpoint_release_matrix_set(pmat_ext)
      CALL tddfpt_kpoint_response_release(kernel_response)
      CALL tddfpt_kpoint_response_release(column_response)
      kernel_env_ptr => kernel_env
      CALL release_kernel_env(kernel_env_ptr)
      CALL qs_rho_release(rho1)
      DEALLOCATE (rho1)
      DEALLOCATE (eigenvectors, evals, inds, kernel_matrix, oscillator_strength, sorted_evals, &
                  state_dipole_im, state_dipole_re, transition_dipole_im, transition_dipole_re, &
                  transitions)
      DEALLOCATE (active_imag, imag_index)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_full_reference

! **************************************************************************************************
!> \brief Build the real-linear response layout.
!> \param kpoints ...
!> \param transitions ...
!> \param active_imag ...
!> \param imag_index ...
!> \param nresponse ...
! **************************************************************************************************
   SUBROUTINE init_response_layout(kpoints, transitions, active_imag, imag_index, nresponse)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      LOGICAL, DIMENSION(:), INTENT(out)                 :: active_imag
      INTEGER, DIMENSION(:), INTENT(out)                 :: imag_index
      INTEGER, INTENT(out)                               :: nresponse

      INTEGER                                            :: ikp, itrans
      REAL(kind=dp), DIMENSION(3)                        :: two_k

      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(SIZE(active_imag) == SIZE(transitions))
      CPASSERT(SIZE(imag_index) == SIZE(transitions))

      nresponse = SIZE(transitions)
      active_imag(:) = .FALSE.
      imag_index(:) = 0
      DO itrans = 1, SIZE(transitions)
         ikp = transitions(itrans)%ikp
         two_k(:) = 2.0_dp*kpoints%xkp(:, ikp)
         active_imag(itrans) = ANY(ABS(two_k(:) - ANINT(two_k(:))) > 1.0E-10_dp)
         IF (active_imag(itrans)) THEN
            nresponse = nresponse + 1
            imag_index(itrans) = nresponse
         END IF
      END DO

   END SUBROUTINE init_response_layout

! **************************************************************************************************
!> \brief Copy imaginary response components into active rows of the real-linear matrix.
!> \param response_matrix ...
!> \param column ...
!> \param sigma_im ...
!> \param active_imag ...
!> \param imag_index ...
! **************************************************************************************************
   SUBROUTINE set_imag_rows(response_matrix, column, sigma_im, active_imag, imag_index)
      REAL(kind=dp), DIMENSION(:, :), INTENT(inout)      :: response_matrix
      INTEGER, INTENT(in)                                :: column
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: sigma_im
      LOGICAL, DIMENSION(:), INTENT(in)                  :: active_imag
      INTEGER, DIMENSION(:), INTENT(in)                  :: imag_index

      INTEGER                                            :: itrans

      DO itrans = 1, SIZE(sigma_im)
         IF (active_imag(itrans)) response_matrix(imag_index(itrans), column) = sigma_im(itrans)
      END DO

   END SUBROUTINE set_imag_rows

! **************************************************************************************************
!> \brief Read internal TDDFPT k-point kernel diagnostics from the environment.
!> \param do_xc ...
!> \param do_hartree ...
!> \param density_weight_power ...
!> \param sigma_weight_power ...
!> \param density_scale ...
!> \param sigma_scale ...
! **************************************************************************************************
   SUBROUTINE read_debug_options(do_xc, do_hartree, density_weight_power, sigma_weight_power, &
                                 density_scale, sigma_scale)
      LOGICAL, INTENT(out)                               :: do_xc, do_hartree
      REAL(kind=dp), INTENT(out)                         :: density_weight_power, &
                                                            sigma_weight_power, density_scale, &
                                                            sigma_scale

      CHARACTER(len=64)                                  :: env_value
      INTEGER                                            :: env_status

      do_xc = .TRUE.
      do_hartree = .TRUE.
      density_weight_power = -0.5_dp
      sigma_weight_power = 0.5_dp
      density_scale = 1.0_dp
      sigma_scale = 1.0_dp

      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TDDFPT_KP_KERNEL_PART", env_value, &
                                    STATUS=env_status)
      IF (env_status == 0) THEN
         CALL uppercase(env_value)
         SELECT CASE (TRIM(ADJUSTL(env_value)))
         CASE ("GAP", "NONE", "ENERGY", "DIAGONAL")
            do_xc = .FALSE.
            do_hartree = .FALSE.
         CASE ("XC")
            do_xc = .TRUE.
            do_hartree = .FALSE.
         CASE ("HARTREE", "COULOMB")
            do_xc = .FALSE.
            do_hartree = .TRUE.
         CASE ("FULL")
            do_xc = .TRUE.
            do_hartree = .TRUE.
         CASE DEFAULT
            CPWARN("Ignoring unknown CP2K_TDDFPT_KP_KERNEL_PART diagnostic value")
         END SELECT
      END IF

      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TDDFPT_KP_DENSITY_WEIGHT_POWER", env_value, &
                                    STATUS=env_status)
      IF (env_status == 0) READ (env_value, *, IOSTAT=env_status) density_weight_power
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TDDFPT_KP_SIGMA_WEIGHT_POWER", env_value, &
                                    STATUS=env_status)
      IF (env_status == 0) READ (env_value, *, IOSTAT=env_status) sigma_weight_power
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TDDFPT_KP_DENSITY_SCALE", env_value, &
                                    STATUS=env_status)
      IF (env_status == 0) READ (env_value, *, IOSTAT=env_status) density_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TDDFPT_KP_SIGMA_SCALE", env_value, &
                                    STATUS=env_status)
      IF (env_status == 0) READ (env_value, *, IOSTAT=env_status) sigma_scale

   END SUBROUTINE read_debug_options

! **************************************************************************************************
!> \brief Build transition dipoles for the explicit k-point transition basis.
!> \param qs_env ...
!> \param tddfpt_control ...
!> \param transitions ...
!> \param ntrans_kpoint ...
!> \param nmo_spin ...
!> \param homo_spin ...
!> \param nvirt_spin ...
!> \param nao_spin ...
!> \param transition_dipole_re ...
!> \param transition_dipole_im ...
! **************************************************************************************************
   SUBROUTINE build_transition_dipoles(qs_env, tddfpt_control, transitions, ntrans_kpoint, &
                                       nmo_spin, homo_spin, nvirt_spin, nao_spin, &
                                       transition_dipole_re, transition_dipole_im)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      INTEGER, INTENT(in)                                :: ntrans_kpoint
      INTEGER, DIMENSION(2), INTENT(in)                  :: nmo_spin, homo_spin, nvirt_spin, nao_spin
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: transition_dipole_re, &
                                                            transition_dipole_im

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_transition_dipoles'

      INTEGER                                            :: gap_index, handle, ideriv, ikp, &
                                                            ikp_local, iocc, ispin, ivirt, nao, &
                                                            nkp, nspins, spin_offset
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: my_kpgrp
      REAL(kind=dp)                                      :: dipole_im, dipole_re, gap, &
                                                            oscillator_factor, spin_factor
      REAL(kind=dp), DIMENSION(:), POINTER               :: wkp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_all
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, moment_struct
      TYPE(cp_fm_type)                                   :: fm_dummy, fm_tmp, mo_coeff_im_global, &
                                                            mo_coeff_re_global, moment_im, &
                                                            moment_re
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: overlap_deriv
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_inter_kp, para_env_kp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp, sab_orb
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, blacs_env_all, cell_to_index, cmatrix, dft_control, fm_struct, kp, &
               kp_env, kpoints, ks_env, mo_coeff_im, mo_coeff_re, moment_struct, mos_kp, &
               overlap_deriv, para_env, para_env_inter_kp, para_env_kp, rmatrix, sab_kp, sab_orb, &
               wkp)
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints, ks_env=ks_env, &
                      sab_orb=sab_orb)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))

      CALL get_kpoint_info(kpoints, nkp=nkp, kp_range=kp_range, kp_env=kp_env, &
                           para_env=para_env, blacs_env_all=blacs_env_all, &
                           para_env_inter_kp=para_env_inter_kp, para_env_kp=para_env_kp, &
                           blacs_env=blacs_env, wkp=wkp, cell_to_index=cell_to_index, &
                           sab_nl=sab_kp)
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(para_env_inter_kp))
      CPASSERT(ASSOCIATED(para_env_kp))
      CPASSERT(ASSOCIATED(blacs_env_all))
      CPASSERT(ASSOCIATED(blacs_env))
      CPASSERT(ASSOCIATED(kp_env))
      CPASSERT(ASSOCIATED(ks_env))
      CPASSERT(ASSOCIATED(sab_orb))
      CPASSERT(ASSOCIATED(sab_kp))
      CPASSERT(ASSOCIATED(cell_to_index))

      nspins = dft_control%nspins
      transition_dipole_re = 0.0_dp
      transition_dipole_im = 0.0_dp

      CALL build_overlap_matrix(ks_env, matrixkp_s=overlap_deriv, nderivative=1, &
                                basis_type_a="ORB", basis_type_b="ORB", sab_nl=sab_orb, &
                                ext_kpoints=kpoints)

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=overlap_deriv(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=overlap_deriv(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_kp)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_kp)

      DO ikp = 1, nkp
         my_kpgrp = (ikp >= kp_range(1) .AND. ikp <= kp_range(2))
         IF (my_kpgrp) THEN
            ikp_local = ikp - kp_range(1) + 1
            kp => kp_env(ikp_local)%kpoint_env
            mos_kp => kp%mos
         ELSE
            NULLIFY (kp, mos_kp)
         END IF
         DO ispin = 1, nspins
            spin_offset = tddfpt_kpoint_transition_offset(ispin, homo_spin, nvirt_spin)
            nao = nao_spin(ispin)

            CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_spin(ispin), &
                                     para_env=para_env, context=blacs_env_all)
            CALL cp_fm_create(mo_coeff_re_global, fm_struct)
            CALL cp_fm_create(mo_coeff_im_global, fm_struct)
            CALL cp_fm_create(fm_tmp, fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            CALL cp_fm_struct_create(moment_struct, nrow_global=nmo_spin(ispin), &
                                     ncol_global=nmo_spin(ispin), para_env=para_env, &
                                     context=blacs_env_all)
            CALL cp_fm_create(moment_re, moment_struct)
            CALL cp_fm_create(moment_im, moment_struct)
            CALL cp_fm_struct_release(moment_struct)

            IF (my_kpgrp) THEN
               CALL get_mo_set(mos_kp(1, ispin), mo_coeff=mo_coeff_re)
               CALL get_mo_set(mos_kp(2, ispin), mo_coeff=mo_coeff_im)
               CPASSERT(ASSOCIATED(mo_coeff_re))
               CPASSERT(ASSOCIATED(mo_coeff_im))
               CALL cp_fm_copy_general(mo_coeff_re, mo_coeff_re_global, para_env)
               CALL cp_fm_copy_general(mo_coeff_im, mo_coeff_im_global, para_env)
            ELSE
               CALL cp_fm_copy_general(fm_dummy, mo_coeff_re_global, para_env)
               CALL cp_fm_copy_general(fm_dummy, mo_coeff_im_global, para_env)
            END IF

            spin_factor = 1.0_dp
            IF (nspins == 1) THEN
               IF (tddfpt_control%rks_triplets) THEN
                  spin_factor = 0.0_dp
               ELSE
                  spin_factor = 2.0_dp
               END IF
            END IF

            DO ideriv = 1, nderivs
               CALL dbcsr_set(rmatrix, 0.0_dp)
               CALL dbcsr_set(cmatrix, 0.0_dp)
               CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=overlap_deriv, &
                                   ispin=ideriv + 1, xkp=kpoints%xkp(:, ikp), &
                                   cell_to_index=cell_to_index, sab_nl=sab_kp)

               ! Project the complex AO operator as C^H A C; the imaginary MO
               ! coefficients on the bra side carry the conjugation signs.
               CALL cp_dbcsr_sm_fm_multiply(rmatrix, mo_coeff_re_global, fm_tmp, nmo_spin(ispin))
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  1.0_dp, mo_coeff_re_global, fm_tmp, 0.0_dp, moment_re)
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  -1.0_dp, mo_coeff_im_global, fm_tmp, 0.0_dp, moment_im)

               CALL cp_dbcsr_sm_fm_multiply(rmatrix, mo_coeff_im_global, fm_tmp, nmo_spin(ispin))
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  1.0_dp, mo_coeff_re_global, fm_tmp, 1.0_dp, moment_im)
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  1.0_dp, mo_coeff_im_global, fm_tmp, 1.0_dp, moment_re)

               CALL cp_dbcsr_sm_fm_multiply(cmatrix, mo_coeff_re_global, fm_tmp, nmo_spin(ispin))
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  1.0_dp, mo_coeff_re_global, fm_tmp, 1.0_dp, moment_im)
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  1.0_dp, mo_coeff_im_global, fm_tmp, 1.0_dp, moment_re)

               CALL cp_dbcsr_sm_fm_multiply(cmatrix, mo_coeff_im_global, fm_tmp, nmo_spin(ispin))
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  -1.0_dp, mo_coeff_re_global, fm_tmp, 1.0_dp, moment_re)
               CALL parallel_gemm("T", "N", nmo_spin(ispin), nmo_spin(ispin), nao, &
                                  1.0_dp, mo_coeff_im_global, fm_tmp, 1.0_dp, moment_im)

               DO iocc = 1, homo_spin(ispin)
                  DO ivirt = homo_spin(ispin) + 1, nmo_spin(ispin)
                     gap_index = (ikp - 1)*ntrans_kpoint + spin_offset + &
                                 (iocc - 1)*nvirt_spin(ispin) + ivirt - homo_spin(ispin)
                     gap = transitions(gap_index)%energy
                     IF (gap <= 0.0_dp) &
                        CPABORT("K-point TDDFPT requires positive occupied-virtual energy gaps")
                     oscillator_factor = SQRT(spin_factor*wkp(ikp))/gap
                     CALL cp_fm_get_element(moment_re, ivirt, iocc, dipole_re)
                     CALL cp_fm_get_element(moment_im, ivirt, iocc, dipole_im)
                     transition_dipole_re(gap_index, ideriv) = oscillator_factor*dipole_re
                     transition_dipole_im(gap_index, ideriv) = oscillator_factor*dipole_im
                  END DO
               END DO
            END DO

            CALL cp_fm_release(moment_im)
            CALL cp_fm_release(moment_re)
            CALL cp_fm_release(fm_tmp)
            CALL cp_fm_release(mo_coeff_im_global)
            CALL cp_fm_release(mo_coeff_re_global)
         END DO
      END DO

      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      CALL dbcsr_deallocate_matrix_set(overlap_deriv)

      CALL timestop(handle)

   END SUBROUTINE build_transition_dipoles

! **************************************************************************************************
!> \brief Apply local/semi-local Hartree/XC kernel to one k-point response vector.
!> \param qs_env ...
!> \param response ...
!> \param kernel_response ...
!> \param transitions ...
!> \param kernel_env ...
!> \param rho1 ...
!> \param pmat_ext ...
!> \param hmat_rs ...
!> \param hmat_ext ...
!> \param kernel_rspace ...
!> \param work_v_xc ...
!> \param work_v_xc_tau ...
!> \param hartree_rspace ...
!> \param work_v_gspace ...
!> \param poisson_env ...
!> \param is_rks_triplets ...
!> \param do_xc ...
!> \param do_hartree ...
!> \param density_weight_power ...
!> \param sigma_weight_power ...
!> \param density_scale ...
!> \param sigma_scale ...
! **************************************************************************************************
   SUBROUTINE apply_local_kernel(qs_env, response, kernel_response, transitions, kernel_env, rho1, &
                                 pmat_ext, hmat_rs, hmat_ext, kernel_rspace, work_v_xc, &
                                 work_v_xc_tau, hartree_rspace, work_v_gspace, poisson_env, &
                                 is_rks_triplets, do_xc, do_hartree, density_weight_power, &
                                 sigma_weight_power, density_scale, sigma_scale)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response, kernel_response
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      TYPE(full_kernel_env_type), INTENT(in)             :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(cp_fm_type), DIMENSION(:, :, :), &
         INTENT(inout)                                   :: pmat_ext
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: hmat_rs
      TYPE(cp_fm_type), DIMENSION(:, :, :), &
         INTENT(inout)                                   :: hmat_ext
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: kernel_rspace, work_v_xc, work_v_xc_tau
      TYPE(pw_r3d_rs_type), POINTER                      :: hartree_rspace
      TYPE(pw_c1d_gs_type), INTENT(inout)                :: work_v_gspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      LOGICAL, INTENT(in)                                :: is_rks_triplets, do_xc, do_hartree
      REAL(kind=dp), INTENT(in)                          :: density_weight_power, &
                                                            sigma_weight_power, density_scale, &
                                                            sigma_scale

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_local_kernel'

      INTEGER                                            :: handle, img, ispin, nspins
      REAL(kind=dp)                                      :: alpha, ehartree
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hmat_spin
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(tddfpt_kpoint_response_type)                  :: density_response

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, hmat_spin, pw_env, rho_g, weights)
      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, xcint_weights=weights)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(pw_env))
      nspins = dft_control%nspins

      IF ((.NOT. do_xc) .AND. ((.NOT. do_hartree) .OR. is_rks_triplets)) THEN
         CALL tddfpt_kpoint_response_zero(kernel_response)
         CALL timestop(handle)
         RETURN
      END IF

      CALL tddfpt_kpoint_response_create(density_response, response%ntrans)
      CALL tddfpt_kpoint_response_copy_scaled_coeff(density_response, response, transitions, &
                                                    density_weight_power)
      density_response%coeff_re(:) = density_scale*density_response%coeff_re(:)
      density_response%coeff_im(:) = density_scale*density_response%coeff_im(:)
      CALL tddfpt_kpoint_response_update_rho(qs_env, density_response, transitions, rho1, pmat_ext)

      DO ispin = 1, nspins
         CALL pw_zero(kernel_rspace(ispin))
      END DO
      IF (do_xc) THEN
         CALL tddfpt_apply_xc(A_ia_rspace=kernel_rspace, kernel_env=kernel_env, &
                              rho_ia_struct=rho1, is_rks_triplets=is_rks_triplets, &
                              pw_env=pw_env, weights=weights, &
                              work_v_xc=work_v_xc, work_v_xc_tau=work_v_xc_tau)
      END IF

      IF (do_hartree .AND. (.NOT. is_rks_triplets)) THEN
         CALL qs_rho_get(rho1, rho_g=rho_g)
         CPASSERT(ASSOCIATED(rho_g))
         DO ispin = 2, nspins
            CALL pw_axpy(rho_g(ispin), rho_g(1))
         END DO
         CALL pw_poisson_solve(poisson_env, rho_g(1), ehartree, work_v_gspace)
         CALL pw_transfer(work_v_gspace, hartree_rspace)
         IF (nspins > 1) THEN
            alpha = 1.0_dp
         ELSE
            alpha = 2.0_dp
         END IF
         DO ispin = 1, nspins
            CALL pw_axpy(hartree_rspace, kernel_rspace(ispin), alpha)
         END DO
      END IF

      DO ispin = 1, nspins
         DO img = 1, SIZE(hmat_rs, 2)
            CALL dbcsr_set(hmat_rs(ispin, img)%matrix, 0.0_dp)
         END DO
         CALL pw_scale(kernel_rspace(ispin), kernel_rspace(ispin)%pw_grid%dvol)
         hmat_spin => hmat_rs(ispin, :)
         CALL integrate_v_rspace(v_rspace=kernel_rspace(ispin), hmat_kp=hmat_spin, &
                                 qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)
      END DO

      CALL tddfpt_kpoint_hmat_from_rs(qs_env, hmat_rs, hmat_ext)
      CALL tddfpt_kpoint_response_zero(kernel_response)
      CALL tddfpt_kpoint_hmat_to_response(qs_env, transitions, hmat_ext, kernel_response)
      kernel_response%sigma_re(:) = sigma_scale*kernel_response%sigma_re(:)
      kernel_response%sigma_im(:) = sigma_scale*kernel_response%sigma_im(:)
      CALL tddfpt_kpoint_response_add_scaled_sigma(response, kernel_response, transitions, &
                                                   sigma_weight_power)

      CALL tddfpt_kpoint_response_release(density_response)

      CALL timestop(handle)

   END SUBROUTINE apply_local_kernel

END MODULE qs_tddfpt2_kpoint_full

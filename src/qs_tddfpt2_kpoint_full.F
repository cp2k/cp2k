!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_kpoint_full
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_constants,                 ONLY: tddfpt_kernel_full
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathlib,                         ONLY: diag_complex
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
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_tddfpt2_kpoint_utils,         ONLY: &
        tddfpt_kpoint_apply_energy_diff, tddfpt_kpoint_create_matrix_set, &
        tddfpt_kpoint_hmat_from_rs, tddfpt_kpoint_hmat_to_response, &
        tddfpt_kpoint_init_transitions, tddfpt_kpoint_release_matrix_set, &
        tddfpt_kpoint_response_create, tddfpt_kpoint_response_release, &
        tddfpt_kpoint_response_set_unit, tddfpt_kpoint_response_type, &
        tddfpt_kpoint_response_update_rho, tddfpt_kpoint_response_zero, &
        tddfpt_kpoint_transition_type
   USE qs_tddfpt2_operators,            ONLY: tddfpt_apply_xc
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_kpoint_full'

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

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenvectors, kernel_matrix
      INTEGER                                            :: handle, ikernel, img, ispin, istate, &
                                                            itrans, log_unit, nimages, nstates, &
                                                            ntrans_kpoint, ntrans_total
      INTEGER, DIMENSION(2)                              :: homo_spin, nao_spin, nmo_spin, nvirt_spin
      REAL(kind=dp)                                      :: checksum
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
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
      nstates = MIN(tddfpt_control%nstates, ntrans_total)
      IF (tddfpt_control%nstates > ntrans_total) THEN
         CPWARN("Requested more TDDFPT states than k-point occupied-virtual transitions")
      END IF

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

      ALLOCATE (kernel_matrix(ntrans_total, ntrans_total), eigenvectors(ntrans_total, ntrans_total), &
                evals(ntrans_total))
      kernel_matrix = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      DO itrans = 1, ntrans_total
         CALL tddfpt_kpoint_response_set_unit(column_response, itrans)
         CALL tddfpt_kpoint_apply_energy_diff(column_response, transitions)
         CALL apply_local_kernel(qs_env, column_response, kernel_response, transitions, kernel_env, &
                                 rho1, pmat_ext, hmat_rs, hmat_ext, kernel_rspace, work_v_xc, &
                                 work_v_xc_tau, hartree_rspace, work_v_gspace, poisson_env, &
                                 tddfpt_control%rks_triplets)
         kernel_matrix(:, itrans) = CMPLX(column_response%sigma_re(:), column_response%sigma_im(:), &
                                          KIND=dp)
      END DO

      kernel_matrix(:, :) = 0.5_dp*(kernel_matrix(:, :) + CONJG(TRANSPOSE(kernel_matrix(:, :))))
      CALL diag_complex(kernel_matrix, eigenvectors, evals)

      checksum = SQRT(SUM(evals(1:nstates)**2))
      log_unit = cp_logger_get_default_io_unit(logger)
      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                 TDDFPT K-point Local/Semilocal Full Kernel                  -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, "(1X,A)") &
            "Explicit full transition-space reference solver; HFX, SOC, smearing and metals are guarded."
         WRITE (log_unit, '(/,T10,A,T19,A,T39,A)') "State", "Excitation", "Dominant transition"
         WRITE (log_unit, '(T10,A,T19,A,T39,A,T48,A,T57,A,T66,A)') "number", "energy (eV)", &
            "k", "spin", "occ", "virt"
         WRITE (log_unit, '(T10,68("-"))')
      END IF

      DO istate = 1, nstates
         ikernel = MAXLOC(ABS(eigenvectors(:, istate)), DIM=1)
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(1X,A,T9,I7,T19,F11.5,T39,I7,T48,I7,T57,I7,T66,I7)') &
               "TDDFPT|", istate, evals(istate)*evolt, transitions(ikernel)%ikp, &
               transitions(ikernel)%ispin, transitions(ikernel)%iocc, transitions(ikernel)%ivirt
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
      DEALLOCATE (eigenvectors, evals, kernel_matrix, transitions)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_full_reference

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
! **************************************************************************************************
   SUBROUTINE apply_local_kernel(qs_env, response, kernel_response, transitions, kernel_env, rho1, &
                                 pmat_ext, hmat_rs, hmat_ext, kernel_rspace, work_v_xc, &
                                 work_v_xc_tau, hartree_rspace, work_v_gspace, poisson_env, &
                                 is_rks_triplets)
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
      LOGICAL, INTENT(in)                                :: is_rks_triplets

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_local_kernel'

      INTEGER                                            :: handle, img, ispin, itrans, nspins
      REAL(kind=dp)                                      :: alpha, ehartree
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hmat_spin
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(tddfpt_kpoint_response_type)                  :: weighted_response

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, hmat_spin, pw_env, rho_g, weights)
      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, xcint_weights=weights)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(pw_env))
      nspins = dft_control%nspins

      CALL tddfpt_kpoint_response_create(weighted_response, response%ntrans)
      weighted_response%coeff_re(:) = SQRT(transitions(:)%weight)*response%coeff_re(:)
      weighted_response%coeff_im(:) = SQRT(transitions(:)%weight)*response%coeff_im(:)
      CALL tddfpt_kpoint_response_update_rho(qs_env, weighted_response, transitions, rho1, pmat_ext)

      DO ispin = 1, nspins
         CALL pw_zero(kernel_rspace(ispin))
      END DO
      CALL tddfpt_apply_xc(A_ia_rspace=kernel_rspace, kernel_env=kernel_env, &
                           rho_ia_struct=rho1, is_rks_triplets=is_rks_triplets, &
                           pw_env=pw_env, weights=weights, &
                           work_v_xc=work_v_xc, work_v_xc_tau=work_v_xc_tau)

      IF (.NOT. is_rks_triplets) THEN
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
      DO itrans = 1, response%ntrans
         response%sigma_re(itrans) = response%sigma_re(itrans) + &
                                     SQRT(transitions(itrans)%weight)*kernel_response%sigma_re(itrans)
         response%sigma_im(itrans) = response%sigma_im(itrans) + &
                                     SQRT(transitions(itrans)%weight)*kernel_response%sigma_im(itrans)
      END DO

      CALL tddfpt_kpoint_response_release(weighted_response)

      CALL timestop(handle)

   END SUBROUTINE apply_local_kernel

END MODULE qs_tddfpt2_kpoint_full

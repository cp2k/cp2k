!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2024 created [Maximilian Graml]
! **************************************************************************************************

MODULE bse_main

   USE bse_full_diag,                   ONLY: create_A_and_B,&
                                              create_hermitian_form_of_ABBA,&
                                              diagonalize_A,&
                                              diagonalize_C
   USE bse_iterative,                   ONLY: solve_bse_iteratively
   USE bse_print,                       ONLY: print_BSE_start_flag
   USE bse_util,                        ONLY: adapt_BSE_input_params,&
                                              bse_window_from_dft,&
                                              deallocate_matrices_bse,&
                                              estimate_BSE_resources,&
                                              get_bse_spin_block_layout,&
                                              invert_static_screening,&
                                              screen_slabs,&
                                              truncate_BSE_matrices
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level
   USE input_constants,                 ONLY: bse_abba,&
                                              bse_both,&
                                              bse_fulldiag,&
                                              bse_iterdiag,&
                                              bse_tda
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mo_window,                       ONLY: mo_window_type
   USE mp2_types,                       ONLY: mp2_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_main'

   PUBLIC :: start_bse_calculation

CONTAINS

! **************************************************************************************************
!> \brief Main subroutine managing BSE calculations
!> \param fm_mat_S_ia_bse ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param bse_lev_virt ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param mo_coeff ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE start_bse_calculation(fm_mat_S_ia_bse, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_Q_static_bse_gemm, &
                                    Eigenval, Eigenval_scf, &
                                    homo, virtual, dimen_RI, dimen_RI_red, bse_lev_virt, &
                                    mp2_env, qs_env, mo_coeff, unit_nr)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S_ia_bse, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_Q_static_bse_gemm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: Eigenval, Eigenval_scf
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      INTEGER, DIMENSION(:), INTENT(IN)                  :: bse_lev_virt
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'start_bse_calculation'

      INTEGER                                            :: handle, ispin, n_ov_joint, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo_red_arr, n_mo, n_ov_arr, &
                                                            offsets_arr, virt_red_arr
      LOGICAL                                            :: my_do_abba, my_do_fulldiag, &
                                                            my_do_iterat_diag, my_do_tda
      REAL(KIND=dp)                                      :: diag_runtime_est
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Eigenval_reduced, Eigenval_reduced_1, &
         Eigenval_reduced_2, Eigenval_reduced_joint, Eigenval_reduced_qp, Eigenval_reduced_scf
      TYPE(cp_fm_type)                                   :: fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                            fm_eps_inv, fm_inv_sqrt_A_minus_B, &
                                                            fm_sqrt_A_minus_B
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: fm_mat_S_ab_trunc_arr, &
         fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr, fm_mat_S_ia_trunc_arr, &
         fm_mat_S_ij_trunc_arr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_window_type)                               :: window
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      para_env => fm_mat_S_ia_bse(1)%matrix_struct%para_env

      ! what the solver reads from QS and the GW, held by the environment (pointers, none owned);
      ! the para_env of the BSE is the RPA's, the one the slabs live on
      mp2_env%bse%nspins = nspins
      mp2_env%bse%dimen_RI = dimen_RI
      mp2_env%bse%dimen_RI_red = dimen_RI_red
      mp2_env%bse%para_env => fm_mat_S_ia_bse(1)%matrix_struct%para_env
      mp2_env%bse%blacs_env => fm_mat_S_ia_bse(1)%matrix_struct%context
      CALL get_qs_env(qs_env, subsys=mp2_env%bse%subsys, pw_env=mp2_env%bse%pw_env, &
                      dft_control=mp2_env%bse%dft_control, mos=mp2_env%bse%mos, &
                      matrix_s=mp2_env%bse%matrix_s, sab_orb=mp2_env%bse%sab_orb, &
                      exstate_env=mp2_env%bse%exstate_env)
      mp2_env%bse%gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
      IF (ALLOCATED(mp2_env%bse%gw_corr_lev_virt)) DEALLOCATE (mp2_env%bse%gw_corr_lev_virt)
      ALLOCATE (mp2_env%bse%gw_corr_lev_virt(nspins))
      mp2_env%bse%gw_corr_lev_virt(:) = bse_lev_virt(:)

      my_do_fulldiag = .FALSE.
      my_do_iterat_diag = .FALSE.
      my_do_tda = .FALSE.
      my_do_abba = .FALSE.
      !Method: Iterative or full diagonalization
      SELECT CASE (mp2_env%bse%bse_diag_method)
      CASE (bse_iterdiag)
         my_do_iterat_diag = .TRUE.
      CASE (bse_fulldiag)
         my_do_fulldiag = .TRUE.
      END SELECT
      !Approximation: TDA and/or full ABBA matrix
      SELECT CASE (mp2_env%bse%flag_tda)
      CASE (bse_tda)
         my_do_tda = .TRUE.
      CASE (bse_abba)
         my_do_abba = .TRUE.
      CASE (bse_both)
         my_do_tda = .TRUE.
         my_do_abba = .TRUE.
      END SELECT

      CALL print_BSE_start_flag(my_do_tda, my_do_abba, unit_nr)

      ! Link BSE debug flag against debug print level
      logger => cp_get_default_logger()
      IF (logger%iter_info%print_level == debug_print_level) THEN
         mp2_env%bse%bse_debug_print = .TRUE.
      END IF

      CALL fm_mat_S_ia_bse(1)%matrix_struct%para_env%sync()

      ALLOCATE (homo_red_arr(nspins), virt_red_arr(nspins))
      ALLOCATE (fm_mat_S_ia_trunc_arr(nspins), fm_mat_S_ij_trunc_arr(nspins), fm_mat_S_ab_trunc_arr(nspins))
      ALLOCATE (fm_mat_S_bar_ia_bse_arr(nspins), fm_mat_S_bar_ij_bse_arr(nspins))

      ALLOCATE (n_mo(nspins))
      n_mo(:) = homo(:) + virtual(:)
      CALL bse_window_from_dft(Eigenval_scf(:, 1, :), n_mo, homo, mp2_env%bse%bse_cutoff_occ, &
                               mp2_env%bse%bse_cutoff_empty, window)
      DEALLOCATE (n_mo)

      IF (nspins > 1) THEN
         IF (my_do_iterat_diag) THEN
            CPABORT("Iterative BSE is not yet available for open-shell references")
         END IF
         CALL cp_warn(__LOCATION__, &
                      "Open-shell (UKS/LSD) BSE is a recent addition and has not been "// &
                      "extensively validated. Verify results carefully before using them "// &
                      "for production calculations.")
         ! === Open-shell path: combined-window cutoff, per-spin truncation, joint A (+B for ABBA) ===
         DO ispin = 1, nspins
            CALL truncate_BSE_matrices(fm_mat_S_ia_bse(ispin), fm_mat_S_ij_bse(ispin), &
                                       fm_mat_S_ab_bse(ispin), &
                                       fm_mat_S_ia_trunc_arr(ispin), fm_mat_S_ij_trunc_arr(ispin), &
                                       fm_mat_S_ab_trunc_arr(ispin), &
                                       Eigenval_scf(:, 1, ispin), Eigenval(:, 1, ispin), &
                                       Eigenval_reduced_scf, Eigenval_reduced_qp, homo(ispin), virtual(ispin), &
                                       dimen_RI, unit_nr, ispin, homo_red_arr(ispin), &
                                       virt_red_arr(ispin), mp2_env%bse, window, .FALSE.)
            ! the truncation returns both energy sets; USE_KS_ENERGIES picks the one the BSE sees
            IF (mp2_env%bse%use_ks_energies) THEN
               CALL MOVE_ALLOC(Eigenval_reduced_scf, Eigenval_reduced)
               DEALLOCATE (Eigenval_reduced_qp)
            ELSE
               CALL MOVE_ALLOC(Eigenval_reduced_qp, Eigenval_reduced)
               DEALLOCATE (Eigenval_reduced_scf)
            END IF
            IF (ispin == 1) THEN
               ALLOCATE (Eigenval_reduced_1(SIZE(Eigenval_reduced)))
               Eigenval_reduced_1(:) = Eigenval_reduced(:)
            ELSE
               ALLOCATE (Eigenval_reduced_2(SIZE(Eigenval_reduced)))
               Eigenval_reduced_2(:) = Eigenval_reduced(:)
            END IF
            DEALLOCATE (Eigenval_reduced)
         END DO
         ! Flat eigenvalue layout: [sigma=1 levels, sigma=2 levels]
         ALLOCATE (Eigenval_reduced_joint(SIZE(Eigenval_reduced_1) + SIZE(Eigenval_reduced_2)))
         Eigenval_reduced_joint(1:SIZE(Eigenval_reduced_1)) = Eigenval_reduced_1
         Eigenval_reduced_joint(SIZE(Eigenval_reduced_1) + 1:) = Eigenval_reduced_2
         DEALLOCATE (Eigenval_reduced_1, Eigenval_reduced_2)

         ALLOCATE (n_ov_arr(nspins), offsets_arr(nspins))
         CALL get_bse_spin_block_layout(homo_red_arr, virt_red_arr, n_ov_arr, offsets_arr, n_ov_joint)

         CALL adapt_BSE_input_params(n_ov_joint, 1, unit_nr, mp2_env%bse)

         ! one [1+Q(0)]^-1 for both spins, then the screened slabs per spin
         CALL invert_static_screening(fm_mat_Q_static_bse_gemm, dimen_RI_red, fm_eps_inv)
         DO ispin = 1, nspins
            CALL screen_slabs(fm_eps_inv, fm_mat_S_ij_trunc_arr(ispin), fm_mat_S_ia_trunc_arr(ispin), &
                              dimen_RI_red, homo_red_arr(ispin), virt_red_arr(ispin), &
                              fm_mat_S_bar_ij_bse_arr(ispin), fm_mat_S_bar_ia_bse_arr(ispin))
         END DO
         CALL cp_fm_release(fm_eps_inv)
         CALL cp_fm_release(fm_mat_Q_static_bse_gemm)

         IF (my_do_fulldiag) THEN
            CALL estimate_BSE_resources(n_ov_joint, unit_nr, my_do_abba, para_env, diag_runtime_est)
            CALL create_A_and_B(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr, &
                                fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr, fm_A_BSE, fm_B_BSE, &
                                my_do_abba, Eigenval_reduced_joint, unit_nr, homo_red_arr, virt_red_arr, &
                                dimen_RI, mp2_env, para_env, qs_env)
            IF (my_do_abba) THEN
               CALL create_hermitian_form_of_ABBA(fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                                  unit_nr, mp2_env, diag_runtime_est)
               CALL cp_fm_release(fm_B_BSE)
            END IF
            NULLIFY (dft_control, tddfpt_control)
            CALL get_qs_env(qs_env, dft_control=dft_control)
            tddfpt_control => dft_control%tddfpt2_control
            ! 4th arg (homo_irred) = full per-spin occupied counts: per-spin absolute-MO labels for
            ! the amplitude table, and N_e = SUM(homo_irred) for the TRK print.
            IF (my_do_tda .AND. (.NOT. tddfpt_control%do_bse)) THEN
               CALL diagonalize_A(fm_A_BSE, homo_red_arr, virt_red_arr, homo, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
            END IF
            CALL cp_fm_release(fm_A_BSE)
            IF (my_do_abba) THEN
               CALL diagonalize_C(fm_C_BSE, homo_red_arr, virt_red_arr, homo, &
                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
               CALL cp_fm_release(fm_C_BSE)
            END IF
         END IF

         DO ispin = 1, nspins
            CALL cp_fm_release(fm_mat_S_bar_ia_bse_arr(ispin))
            CALL cp_fm_release(fm_mat_S_bar_ij_bse_arr(ispin))
            CALL cp_fm_release(fm_mat_S_ia_trunc_arr(ispin))
            CALL cp_fm_release(fm_mat_S_ij_trunc_arr(ispin))
            CALL cp_fm_release(fm_mat_S_ab_trunc_arr(ispin))
         END DO
         IF (mp2_env%bse%do_nto_analysis) DEALLOCATE (mp2_env%bse%bse_nto_state_list_final)
         DEALLOCATE (Eigenval_reduced_joint, n_ov_arr, offsets_arr)

      ELSE
         ! === Closed-shell n_spin=1 path (bit-identical) ===
         CALL truncate_BSE_matrices(fm_mat_S_ia_bse(1), fm_mat_S_ij_bse(1), fm_mat_S_ab_bse(1), &
                                    fm_mat_S_ia_trunc_arr(1), fm_mat_S_ij_trunc_arr(1), &
                                    fm_mat_S_ab_trunc_arr(1), &
                                    Eigenval_scf(:, 1, 1), Eigenval(:, 1, 1), Eigenval_reduced_scf, &
                                    Eigenval_reduced_qp, homo(1), virtual(1), dimen_RI, unit_nr, &
                                    1, homo_red_arr(1), virt_red_arr(1), mp2_env%bse, &
                                    window, .TRUE.)
         ! the truncation returns both energy sets; USE_KS_ENERGIES picks the one the BSE sees
         IF (mp2_env%bse%use_ks_energies) THEN
            CALL MOVE_ALLOC(Eigenval_reduced_scf, Eigenval_reduced)
            DEALLOCATE (Eigenval_reduced_qp)
         ELSE
            CALL MOVE_ALLOC(Eigenval_reduced_qp, Eigenval_reduced)
            DEALLOCATE (Eigenval_reduced_scf)
         END IF
         CALL invert_static_screening(fm_mat_Q_static_bse_gemm, dimen_RI_red, fm_eps_inv)
         CALL screen_slabs(fm_eps_inv, fm_mat_S_ij_trunc_arr(1), fm_mat_S_ia_trunc_arr(1), &
                           dimen_RI_red, homo_red_arr(1), virt_red_arr(1), &
                           fm_mat_S_bar_ij_bse_arr(1), fm_mat_S_bar_ia_bse_arr(1))
         CALL cp_fm_release(fm_eps_inv)

         CALL adapt_BSE_input_params(homo_red_arr(1), virt_red_arr(1), unit_nr, mp2_env%bse)

         IF (my_do_iterat_diag) THEN
            CALL solve_bse_iteratively(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, &
                                       fm_mat_S_ab_trunc_arr, fm_mat_S_bar_ia_bse_arr, &
                                       fm_mat_S_bar_ij_bse_arr, Eigenval_reduced, &
                                       homo_red_arr, virt_red_arr, homo, dimen_RI, &
                                       my_do_tda, my_do_abba, mp2_env, qs_env, mo_coeff, unit_nr)
         END IF

         IF (my_do_fulldiag) THEN
            n_ov_joint = homo_red_arr(1)*virt_red_arr(1)
            CALL estimate_BSE_resources(n_ov_joint, unit_nr, my_do_abba, para_env, diag_runtime_est)
            CALL create_A_and_B(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr, &
                                fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr, fm_A_BSE, fm_B_BSE, &
                                my_do_abba, Eigenval_reduced, unit_nr, homo_red_arr, virt_red_arr, &
                                dimen_RI, mp2_env, para_env, qs_env)
            IF (my_do_abba) THEN
               CALL create_hermitian_form_of_ABBA(fm_A_BSE, fm_B_BSE, fm_C_BSE, &
                                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                                  unit_nr, mp2_env, diag_runtime_est)
            END IF
            CALL cp_fm_release(fm_B_BSE)

            NULLIFY (dft_control, tddfpt_control)
            CALL get_qs_env(qs_env, dft_control=dft_control)
            tddfpt_control => dft_control%tddfpt2_control
            IF (my_do_tda .AND. (.NOT. tddfpt_control%do_bse)) THEN
               CALL diagonalize_A(fm_A_BSE, homo_red_arr, virt_red_arr, homo, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
            END IF
            CALL cp_fm_release(fm_A_BSE)
            IF (my_do_abba) THEN
               CALL diagonalize_C(fm_C_BSE, homo_red_arr, virt_red_arr, homo, &
                                  fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                  unit_nr, diag_runtime_est, mp2_env, qs_env, mo_coeff)
            END IF
            CALL cp_fm_release(fm_C_BSE)
         END IF

         CALL deallocate_matrices_bse(fm_mat_S_bar_ia_bse_arr(1), fm_mat_S_bar_ij_bse_arr(1), &
                                      fm_mat_S_ia_trunc_arr(1), fm_mat_S_ij_trunc_arr(1), &
                                      fm_mat_S_ab_trunc_arr(1), fm_mat_Q_static_bse_gemm, mp2_env%bse)
         DEALLOCATE (Eigenval_reduced)

      END IF

      DEALLOCATE (homo_red_arr, virt_red_arr)
      DEALLOCATE (fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr)
      DEALLOCATE (fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A53)') 'BSE|', 'The BSE was successfully calculated. Have a nice day!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE start_bse_calculation

END MODULE bse_main

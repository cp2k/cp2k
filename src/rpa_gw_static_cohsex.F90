!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Static Coulomb-hole screened-exchange (COHSEX) quasiparticle energies and MOs in the
!>        O(N^4) (RI_RPA%GW) GW code for molecules
!>
!>        Symbols
!>          μ,ν        atomic orbitals (AOs)
!>          P,Q        RI auxiliary functions (after the RI-basis reduction, dimen_RI)
!>          i / a      occupied / virtual Kohn-Sham MOs,  p  any MO,  m,n  MOs in the window of
!>                     corrected MOs [homo-CORR_MOS_OCC+1, homo+CORR_MOS_VIRT]
!>          ϵ_n        DFT eigenvalue,  ϵ_F  Fermi level (between HOMO and LUMO)
!>          B^P_nm     = Σ_Q (P|Q)^(-1/2) (Q|nm)            RI 3-center integrals (fm_mat_S_gw)
!>
!>        Static screened Coulomb interaction (RPA, ω = 0), closed shell:
!>          Q_PQ(0)    = 4 Σ_ia B^P_ia B^Q_ia / (ϵ_a - ϵ_i)
!>          W^c_PQ(0)  = [1 + Q(0)]^(-1)_PQ - δ_PQ           (-δ_PQ removes the bare Coulomb part)
!>          W^c_mp,pn  = Σ_PQ B^P_mp W^c_PQ(0) B^Q_pn
!>
!>        static COHSEX self-energy = Coulomb hole (COH) + screened exchange (SEX):
!>          Σ^COH_mn   =  1/2 Σ_p     W^c_mp,pn
!>          Σ^SEX_mn   = -    Σ_i^occ W^c_mi,in              (+ bare exchange Σ^x_mn)
!>          Σ^c_mn     =  Σ^COH_mn + Σ^SEX_mn = 1/2 Σ_p sgn(ϵ_p - ϵ_F) W^c_mp,pn
!>          Σ^x_mn     = -Σ_i^occ Σ_P B^P_mi B^P_in          (RI_SIGMA_X; else from HFX)
!>
!>        Quasiparticle Hamiltonian in the window:
!>          H^QP_mn    = ϵ_n δ_mn + Σ^x_mn - v^xc_mn + α_mn Σ^c_mn
!>          α_mn       = α_IP = COHSEX_SCALE_OCC (m,n occupied),  α_EA = COHSEX_SCALE_VIRT (m,n virtual),  H^QP_ia = 0
!>          ϵ_n^DG     = H^QP_nn                              (COHSEX DIAGONAL)
!>          H^QP U     = U ϵ^CN, occ-occ and virt-virt block separately   (COHSEX QP_UPDATE)
!>          C^QP_μn    = Σ_m C_μm U_mn                        (COHSEX QP_UPDATE MOs)
!>        COHSEX DIAGONAL prints ϵ^DG only; COHSEX QP_UPDATE also diagonalizes H^QP, prints ϵ^CN
!>        and the canonical MOs.
!>
!> \par History
!>      09.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE rpa_gw_static_cohsex
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Hedin1965,&
                                              Hybertsen1986,&
                                              Tyagi2024,&
                                              Tyagi2026,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_release_p,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_add_on_diag,&
                                              cp_fm_column_scale
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_redistribution,            ONLY: cp_fm_permute
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_constants,                 ONLY: cohsex_qp_update
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_io,                        ONLY: write_mo_set_to_output_unit
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_static_cohsex'

   PUBLIC :: static_cohsex_qp_energies, rotate_bse_slab, get_static_cohsex_mo_coeff

CONTAINS

! **************************************************************************************************
!> \brief Compute and print the diagonal (ϵ^DG) static COHSEX quasiparticle energies; with
!>        COHSEX QP_UPDATE also the canonical (ϵ^CN) energies and the canonical MOs
!> \param fm_mat_Q on entry the spin-summed 1 + Q_PQ(ω=0); on exit W^c_PQ(0) = [1 + Q(0)]^(-1)_PQ - δ_PQ
!> \param fm_mat_S_gw per-spin B^P_mp
!> \param fm_mat_S_gw_work ...
!> \param dimen_RI ...
!> \param nmo ...
!> \param homo ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param Eigenval ...
!> \param mo_coeff ...
!> \param unit_nr ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param eigenval_dg ϵ_n^DG of all MOs (scissor outside the window), nmo x nspins
!> \param eigenval_cn ϵ_n^CN of all MOs (scissor outside the window), nmo x nspins, QP_UPDATE only
!> \param u_qp_all U_mn of all MOs (δ_mn outside the window), nmo x nmo x nspins, QP_UPDATE only
! **************************************************************************************************
   SUBROUTINE static_cohsex_qp_energies(fm_mat_Q, fm_mat_S_gw, fm_mat_S_gw_work, dimen_RI, nmo, homo, &
                                        gw_corr_lev_occ, gw_corr_lev_virt, Eigenval, mo_coeff, unit_nr, &
                                        qs_env, mp2_env, eigenval_dg, eigenval_cn, u_qp_all)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_Q
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S_gw, fm_mat_S_gw_work
      INTEGER, INTENT(IN)                                :: dimen_RI, nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, gw_corr_lev_occ, gw_corr_lev_virt
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type)                                     :: mp2_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: eigenval_dg, eigenval_cn
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: u_qp_all

      CHARACTER(LEN=*), PARAMETER :: routineN = 'static_cohsex_qp_energies'

      CHARACTER(LEN=10)                                  :: spin_label
      CHARACTER(LEN=3)                                   :: occ_vir
      CHARACTER(LEN=40)                                  :: cohsex_label
      INTEGER                                            :: handle, ispin, m, n, n_mo_win, no, &
                                                            nspins, nv, p, win_start
      LOGICAL                                            :: qp_update
      REAL(KIND=dp)                                      :: shift_occ, shift_virt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_can, e_diag, e_qp, scale_vec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_qp, sigma_c, sigma_x, &
                                                            sigma_x_minus_vxc, u_qp, vxc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_B_re, fm_struct_win
      TYPE(cp_fm_type)                                   :: fm_B_re, fm_sigma, fm_WB_re, fm_work
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_mat_S_gw)
      ! COHSEX QP_UPDATE: diagonalize H^QP (canonical energies and MOs); DIAGONAL: ϵ^DG only
      qp_update = (mp2_env%ri_g0w0%cohsex_flavor == cohsex_qp_update)

      ! static COHSEX; with QP_UPDATE the quasiparticle update; scaled self-energy if requested
      CALL cite_reference(Hedin1965)
      CALL cite_reference(Hybertsen1986)
      IF (qp_update) CALL cite_reference(Tyagi2024)
      IF (mp2_env%ri_g0w0%cohsex_scale_occ /= 1.0_dp .OR. mp2_env%ri_g0w0%cohsex_scale_virt /= 1.0_dp) THEN
         CALL cite_reference(Tyagi2026)
      END IF
      CALL cp_fm_get_info(mo_coeff(1), para_env=para_env)

      ! ==============================================================================
      ! 1. Static correlation part of the screened Coulomb interaction
      !      ε_PQ(0)   = δ_PQ + Q_PQ(0)
      !      W^c_PQ(0) = ε^(-1)_PQ(0) - δ_PQ
      ! ==============================================================================
      CALL cp_fm_cholesky_decompose(fm_mat_Q, n=dimen_RI)
      CALL cp_fm_create(fm_work, fm_mat_Q%matrix_struct)
      CALL cp_fm_cholesky_invert(fm_mat_Q, n=dimen_RI, work=fm_work)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_add_on_diag(fm_mat_Q, -1.0_dp, n_active=dimen_RI)

      DO ispin = 1, nspins

         no = gw_corr_lev_occ(ispin)
         nv = gw_corr_lev_virt(ispin)
         n_mo_win = no + nv
         win_start = homo(ispin) - no + 1

         ! ==============================================================================
         ! 2. Screened 3-center integrals, weighted with the COHSEX occupation factor
         !      (W^c B)^P_mp = Σ_Q W^c_PQ(0) B^Q_mp  s_p,   s_p = 1/2 sgn(ϵ_p - ϵ_F)
         !      s_p = -1/2 for p occupied,  s_p = +1/2 for p virtual (COH)
         ! ==============================================================================
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=n_mo_win*nmo, k=dimen_RI, &
                            alpha=1.0_dp, matrix_a=fm_mat_Q, matrix_b=fm_mat_S_gw(ispin), &
                            beta=0.0_dp, matrix_c=fm_mat_S_gw_work(ispin))
         ALLOCATE (scale_vec(n_mo_win*nmo))
         DO m = 1, n_mo_win
            DO p = 1, nmo
               IF (p <= homo(ispin)) THEN
                  scale_vec((m - 1)*nmo + p) = -0.5_dp
               ELSE
                  scale_vec((m - 1)*nmo + p) = 0.5_dp
               END IF
            END DO
         END DO
         CALL cp_fm_column_scale(fm_mat_S_gw_work(ispin), scale_vec)
         DEALLOCATE (scale_vec)

         ! ==============================================================================
         ! 3. Reorder the compound index for O(N^4) contractions over (p,P)
         !      B'_(p,P),m      = B^P_mp
         !      (W^c B)'_(p,P),m = (W^c B)^P_mp
         !    row index (p-1)*dimen_RI + P, P fast: the occupied p (p <= homo) are the
         !    first homo*dimen_RI rows. cp_fm_permute maps the input indices
         !    (P,1 ; m,p) with reordering [4,1,3,2] to (p,P ; m,1).
         ! ==============================================================================
         NULLIFY (fm_struct_B_re)
         CALL cp_fm_struct_create(fm_struct_B_re, template_fmstruct=fm_mat_S_gw(ispin)%matrix_struct, &
                                  nrow_global=nmo*dimen_RI, ncol_global=n_mo_win)
         CALL cp_fm_create(fm_B_re, fm_struct_B_re, name="static_cohsex_B_reordered", set_zero=.TRUE.)
         CALL cp_fm_create(fm_WB_re, fm_struct_B_re, name="static_cohsex_WB_reordered", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct_B_re)
         CALL cp_fm_permute(fm_B_re, fm_mat_S_gw(ispin), 1.0_dp, 1, nmo, dimen_RI, 1, [4, 1, 3, 2])
         CALL cp_fm_permute(fm_WB_re, fm_mat_S_gw_work(ispin), 1.0_dp, 1, nmo, dimen_RI, 1, [4, 1, 3, 2])

         NULLIFY (fm_struct_win)
         CALL cp_fm_struct_create(fm_struct_win, template_fmstruct=fm_mat_S_gw(ispin)%matrix_struct, &
                                  nrow_global=n_mo_win, ncol_global=n_mo_win)
         CALL cp_fm_create(fm_sigma, fm_struct_win)
         CALL cp_fm_struct_release(fm_struct_win)

         ALLOCATE (sigma_c(n_mo_win, n_mo_win), sigma_x(n_mo_win, n_mo_win))

         ! ==============================================================================
         ! 4. Correlation self-energy
         !      Σ^c_mn = Σ_p s_p Σ_PQ B^P_mp W^c_PQ(0) B^Q_pn
         !             = Σ_(p,P) (W^c B)'_(p,P),m B'_(p,P),n
         ! ==============================================================================
         CALL parallel_gemm(transa="T", transb="N", m=n_mo_win, n=n_mo_win, k=nmo*dimen_RI, &
                            alpha=1.0_dp, matrix_a=fm_WB_re, matrix_b=fm_B_re, &
                            beta=0.0_dp, matrix_c=fm_sigma)
         CALL cp_fm_get_submatrix(fm_sigma, sigma_c)

         ! ==============================================================================
         ! 5. RI exchange self-energy (RI_SIGMA_X)
         !      Σ^x_mn = -Σ_i^occ Σ_P B^P_mi B^P_in = -Σ_(i,P) B'_(i,P),m B'_(i,P),n
         ! ==============================================================================
         sigma_x(:, :) = 0.0_dp
         IF (mp2_env%ri_g0w0%do_ri_Sigma_x) THEN
            CALL parallel_gemm(transa="T", transb="N", m=n_mo_win, n=n_mo_win, k=homo(ispin)*dimen_RI, &
                               alpha=-1.0_dp, matrix_a=fm_B_re, matrix_b=fm_B_re, &
                               beta=0.0_dp, matrix_c=fm_sigma)
            CALL cp_fm_get_submatrix(fm_sigma, sigma_x)
         END IF

         CALL cp_fm_release(fm_sigma)
         CALL cp_fm_release(fm_B_re)
         CALL cp_fm_release(fm_WB_re)

         ! ==============================================================================
         ! 6. v^xc and Σ^x - v^xc in the MO window
         !      v^xc_mn         = Σ_μν C_μm v^xc_μν C_νn
         !      (Σ^x - v^xc)_mn = Σ_μν C_μm (Σ^x - v^xc)_μν C_νn
         !    The AO matrix (Σ^x - v^xc)_μν holds the HFX exchange Σ^x_μν without
         !    RI_SIGMA_X and only -v^xc_μν with RI_SIGMA_X, hence in both cases
         !      Σ^x_mn = Σ^x,RI_mn + (Σ^x - v^xc)_mn + v^xc_mn   (Σ^x,RI = 0 without RI)
         ! ==============================================================================
         ALLOCATE (vxc(n_mo_win, n_mo_win), sigma_x_minus_vxc(n_mo_win, n_mo_win))
         CALL ao_to_mo_window(mp2_env%ri_g0w0%matrix_vxc(ispin)%matrix, mo_coeff(ispin), &
                              win_start, n_mo_win, vxc)
         CALL ao_to_mo_window(mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc(ispin)%matrix, mo_coeff(ispin), &
                              win_start, n_mo_win, sigma_x_minus_vxc)
         sigma_x(:, :) = sigma_x(:, :) + sigma_x_minus_vxc(:, :) + vxc(:, :)

         ! ==============================================================================
         ! 7. static COHSEX quasiparticle Hamiltonian in the window
         !      H^QP_ij = ϵ_i δ_ij + Σ^x_ij - v^xc_ij + α_IP Σ^c_ij   (occ-occ)
         !      H^QP_ab = ϵ_a δ_ab + Σ^x_ab - v^xc_ab + α_EA Σ^c_ab   (virt-virt)
         !      H^QP_ia = 0                                              (occ-virt)
         !    symmetrized, H^QP = (H^QP + H^QP^T)/2
         !      ϵ_n^DG  = H^QP_nn                                        (diagonal static COHSEX)
         ! ==============================================================================
         ALLOCATE (h_qp(n_mo_win, n_mo_win))
         h_qp(:, :) = 0.0_dp
         h_qp(1:no, 1:no) = sigma_x(1:no, 1:no) - vxc(1:no, 1:no) + &
                            mp2_env%ri_g0w0%cohsex_scale_occ*sigma_c(1:no, 1:no)
         h_qp(no + 1:, no + 1:) = sigma_x(no + 1:, no + 1:) - vxc(no + 1:, no + 1:) + &
                                  mp2_env%ri_g0w0%cohsex_scale_virt*sigma_c(no + 1:, no + 1:)
         h_qp(:, :) = 0.5_dp*(h_qp(:, :) + TRANSPOSE(h_qp(:, :)))
         DO n = 1, n_mo_win
            h_qp(n, n) = h_qp(n, n) + Eigenval(win_start + n - 1, ispin)
         END DO

         ALLOCATE (e_diag(n_mo_win))
         DO n = 1, n_mo_win
            e_diag(n) = h_qp(n, n)
         END DO

         ! ϵ^DG of all MOs for a subsequent BSE, scissor shift of the DFT eigenvalue outside the window
         IF (PRESENT(eigenval_dg)) THEN
            DO n = 1, nmo
               IF (n < win_start) THEN
                  eigenval_dg(n, ispin) = Eigenval(n, ispin) + e_diag(1) - Eigenval(win_start, ispin)
               ELSE IF (n >= win_start + n_mo_win) THEN
                  eigenval_dg(n, ispin) = Eigenval(n, ispin) + e_diag(n_mo_win) - &
                                          Eigenval(win_start + n_mo_win - 1, ispin)
               ELSE
                  eigenval_dg(n, ispin) = e_diag(n - win_start + 1)
               END IF
            END DO
         END IF

         IF (qp_update) THEN
            ! ==============================================================================
            ! 8. Canonical static COHSEX (QP_UPDATE): diagonalize the occ-occ and virt-virt blocks
            !      Σ_m H^QP_nm U_mk = U_nk ϵ_k^CN
            ! ==============================================================================
            ALLOCATE (e_can(n_mo_win))
            e_can(:) = 0.0_dp
            IF (para_env%is_source()) THEN
               IF (no > 0) CALL diamat_all(h_qp(1:no, 1:no), e_can(1:no))
               IF (nv > 0) CALL diamat_all(h_qp(no + 1:, no + 1:), e_can(no + 1:))
            END IF
            CALL para_env%bcast(h_qp, para_env%source)
            CALL para_env%bcast(e_can, para_env%source)

            ! ==============================================================================
            ! 9. Energies and MO transformation for all nmo MOs
            !      inside the window:   ϵ_n^QP = ϵ_n^CN,  U_mn from step 8
            !      outside the window:  U_mn = δ_mn and scissor shift of the DFT eigenvalue
            !        ϵ_n^QP = ϵ_n + (ϵ^CN_first - ϵ_first)   below the window (occupied)
            !        ϵ_n^QP = ϵ_n + (ϵ^CN_last  - ϵ_last)    above the window (virtual)
            ! ==============================================================================
            ALLOCATE (e_qp(nmo), u_qp(nmo, nmo))
            shift_occ = e_can(1) - Eigenval(win_start, ispin)
            shift_virt = e_can(n_mo_win) - Eigenval(win_start + n_mo_win - 1, ispin)
            DO n = 1, nmo
               IF (n < win_start) THEN
                  e_qp(n) = Eigenval(n, ispin) + shift_occ
               ELSE IF (n >= win_start + n_mo_win) THEN
                  e_qp(n) = Eigenval(n, ispin) + shift_virt
               ELSE
                  e_qp(n) = e_can(n - win_start + 1)
               END IF
            END DO
            u_qp(:, :) = 0.0_dp
            DO n = 1, nmo
               u_qp(n, n) = 1.0_dp
            END DO
            u_qp(win_start:win_start + n_mo_win - 1, win_start:win_start + n_mo_win - 1) = h_qp(:, :)

            ! for a subsequent BSE: ϵ^CN and U of all MOs
            IF (PRESENT(eigenval_cn)) eigenval_cn(:, ispin) = e_qp(:)
            IF (PRESENT(u_qp_all)) u_qp_all(:, :, ispin) = u_qp(:, :)
         END IF

         ! ==============================================================================
         ! 10. Print per MO: ϵ_n, Σ^c_nn (unscaled), Σ^x_nn, v^xc_nn, ϵ_n^DG and, with QP_UPDATE,
         !     ϵ_n^CN; then HOMO, LUMO and gap. ϵ^CN is sorted in each block, i.e.
         !     ϵ^CN_HOMO = max over the occ block and ϵ^CN_LUMO = min over the virt block
         ! ==============================================================================
         IF (unit_nr > 0) THEN
            spin_label = ""
            IF (nspins == 2 .AND. ispin == 1) spin_label = "Alpha "
            IF (nspins == 2 .AND. ispin == 2) spin_label = "Beta "
            WRITE (unit_nr, "(/,T3,A)") REPEAT("-", 77)
            ! "static COHSEX" for closed shell, "Alpha static COHSEX" / "Beta static COHSEX" otherwise
            cohsex_label = TRIM(ADJUSTL(TRIM(spin_label)//" static COHSEX"))
            WRITE (unit_nr, "(T3,A)") TRIM(cohsex_label)//" quasiparticle energies"
            WRITE (unit_nr, "(T3,A,F6.3,A,F6.3)") "COHSEX_SCALE_OCC =", mp2_env%ri_g0w0%cohsex_scale_occ, &
               ",  COHSEX_SCALE_VIRT =", mp2_env%ri_g0w0%cohsex_scale_virt
            WRITE (unit_nr, "(T3,A)") REPEAT("-", 77)
            IF (qp_update) THEN
               WRITE (unit_nr, "(T3,A)") "DG: diagonal of H^QP,  CN: canonical, eigenvalues of the "// &
                  "occ-occ and virt-virt blocks"
               WRITE (unit_nr, "(/,T3,A,T16,A,T27,A,T38,A,T49,A,T60,A,T71,A)") "MO", "e_n^DFT", "Sigma^c_nn", &
                  "Sigma^x_nn", "v^xc_nn", "e_n^DG", "e_n^CN"
               WRITE (unit_nr, "(T3,A,T16,A,T27,A,T38,A,T49,A,T60,A,T71,A)") "", "(eV)", "(eV)", "(eV)", &
                  "(eV)", "(eV)", "(eV)"
            ELSE
               WRITE (unit_nr, "(T3,A)") "DG: diagonal of H^QP"
               WRITE (unit_nr, "(/,T3,A,T16,A,T27,A,T38,A,T49,A,T60,A)") "MO", "e_n^DFT", "Sigma^c_nn", &
                  "Sigma^x_nn", "v^xc_nn", "e_n^DG"
               WRITE (unit_nr, "(T3,A,T16,A,T27,A,T38,A,T49,A,T60,A)") "", "(eV)", "(eV)", "(eV)", &
                  "(eV)", "(eV)"
            END IF
            DO n = 1, n_mo_win
               IF (n <= no) THEN
                  occ_vir = "occ"
               ELSE
                  occ_vir = "vir"
               END IF
               IF (qp_update) THEN
                  WRITE (unit_nr, "(T3,I5,A,A3,A,6F11.3)") win_start + n - 1, " (", occ_vir, ")", &
                     Eigenval(win_start + n - 1, ispin)*evolt, sigma_c(n, n)*evolt, sigma_x(n, n)*evolt, &
                     vxc(n, n)*evolt, e_diag(n)*evolt, e_can(n)*evolt
               ELSE
                  WRITE (unit_nr, "(T3,I5,A,A3,A,5F11.3)") win_start + n - 1, " (", occ_vir, ")", &
                     Eigenval(win_start + n - 1, ispin)*evolt, sigma_c(n, n)*evolt, sigma_x(n, n)*evolt, &
                     vxc(n, n)*evolt, e_diag(n)*evolt
               END IF
            END DO
            IF (no > 0 .AND. nv > 0) THEN
               WRITE (unit_nr, "(/,T3,A,T66,F15.3)") TRIM(cohsex_label)//" HOMO, diagonal (eV)", &
                  MAXVAL(e_diag(1:no))*evolt
               WRITE (unit_nr, "(T3,A,T66,F15.3)") TRIM(cohsex_label)//" LUMO, diagonal (eV)", &
                  MINVAL(e_diag(no + 1:))*evolt
               WRITE (unit_nr, "(T3,A,T66,F15.3)") TRIM(cohsex_label)//" HOMO-LUMO gap, diagonal (eV)", &
                  (MINVAL(e_diag(no + 1:)) - MAXVAL(e_diag(1:no)))*evolt
               IF (qp_update) THEN
                  WRITE (unit_nr, "(T3,A,T66,F15.3)") TRIM(cohsex_label)//" HOMO, canonical (eV)", &
                     MAXVAL(e_can(1:no))*evolt
                  WRITE (unit_nr, "(T3,A,T66,F15.3)") TRIM(cohsex_label)//" LUMO, canonical (eV)", &
                     MINVAL(e_can(no + 1:))*evolt
                  WRITE (unit_nr, "(T3,A,T66,F15.3)") TRIM(cohsex_label)//" HOMO-LUMO gap, canonical (eV)", &
                     (MINVAL(e_can(no + 1:)) - MAXVAL(e_can(1:no)))*evolt
               END IF
            END IF
            WRITE (unit_nr, "(T3,A)") REPEAT("-", 77)
         END IF

         IF (qp_update) THEN
            ! ==============================================================================
            ! 11. Canonical static COHSEX MOs C^QP = C U via DFT%PRINT%MO and DFT%PRINT%MO_CUBES
            ! ==============================================================================
            CALL print_static_cohsex_mos(qs_env, mo_coeff(ispin), u_qp, e_qp, homo(ispin), nspins, ispin)
            DEALLOCATE (e_can, e_qp, u_qp)
         END IF

         DEALLOCATE (sigma_c, sigma_x, sigma_x_minus_vxc, vxc, h_qp, e_diag)

      END DO

      DO ispin = 1, nspins
         CALL dbcsr_release_p(mp2_env%ri_g0w0%matrix_vxc(ispin)%matrix)
         CALL dbcsr_release_p(mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc(ispin)%matrix)
      END DO
      DEALLOCATE (mp2_env%ri_g0w0%matrix_vxc, mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc)

      CALL timestop(handle)

   END SUBROUTINE static_cohsex_qp_energies

! **************************************************************************************************
!> \brief MO-window block of an AO matrix
!>          A_mn = Σ_μν C_μm A_μν C_νn,   m,n = win_start, ..., win_start + n_mo_win - 1
!> \param matrix_ao ...
!> \param mo_coeff ...
!> \param win_start ...
!> \param n_mo_win ...
!> \param mat_mo ...
! **************************************************************************************************
   SUBROUTINE ao_to_mo_window(matrix_ao, mo_coeff, win_start, n_mo_win, mat_mo)
      TYPE(dbcsr_type), POINTER                          :: matrix_ao
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff
      INTEGER, INTENT(IN)                                :: win_start, n_mo_win
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: mat_mo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ao_to_mo_window'

      INTEGER                                            :: handle, nao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_ao, fm_ao_mo, fm_mo_mo

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(mo_coeff, nrow_global=nao)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=mo_coeff%matrix_struct, &
                               nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(fm_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(matrix_ao, fm_ao)

      ! (A C)_μn = Σ_ν A_μν C_νn,  only the window columns n of C

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=mo_coeff%matrix_struct, &
                               nrow_global=nao, ncol_global=n_mo_win)
      CALL cp_fm_create(fm_ao_mo, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL parallel_gemm(transa="N", transb="N", m=nao, n=n_mo_win, k=nao, alpha=1.0_dp, &
                         matrix_a=fm_ao, matrix_b=mo_coeff, beta=0.0_dp, matrix_c=fm_ao_mo, &
                         b_first_col=win_start)
      CALL cp_fm_release(fm_ao)

      ! A_mn = Σ_μ C_μm (A C)_μn
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=mo_coeff%matrix_struct, &
                               nrow_global=n_mo_win, ncol_global=n_mo_win)
      CALL cp_fm_create(fm_mo_mo, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL parallel_gemm(transa="T", transb="N", m=n_mo_win, n=n_mo_win, k=nao, alpha=1.0_dp, &
                         matrix_a=mo_coeff, matrix_b=fm_ao_mo, beta=0.0_dp, matrix_c=fm_mo_mo, &
                         a_first_col=win_start)
      CALL cp_fm_release(fm_ao_mo)

      CALL cp_fm_get_submatrix(fm_mo_mo, mat_mo)
      CALL cp_fm_release(fm_mo_mo)

      CALL timestop(handle)

   END SUBROUTINE ao_to_mo_window

! **************************************************************************************************
!> \brief Print the canonical static COHSEX MOs with the DFT%PRINT%MO and DFT%PRINT%MO_CUBES
!> \param qs_env ...
!> \param mo_coeff ...
!> \param u_qp ...
!> \param e_qp ...
!> \param homo ...
!> \param nspins ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE print_static_cohsex_mos(qs_env, mo_coeff, u_qp, e_qp, homo, nspins, ispin)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: u_qp
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: e_qp
      INTEGER, INTENT(IN)                                :: homo, nspins, ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_static_cohsex_mos'

      CHARACTER(LEN=default_path_length)                 :: fname, title
      INTEGER                                            :: first_mo, handle, ivector, last_mo, nao, &
                                                            nelectron, nhomo, nlumo, nmo, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL                                            :: print_coeff, print_cubes, &
                                                            print_energies, print_mo, print_occ
      REAL(KIND=dp)                                      :: maxocc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_qp_mo
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type)                                  :: qp_mo_set
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: wf_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      logger => cp_get_default_logger()

      CALL section_vals_val_get(dft_section, "PRINT%MO%ENERGIES", l_val=print_energies)
      CALL section_vals_val_get(dft_section, "PRINT%MO%COEFFICIENTS", l_val=print_coeff)
      CALL section_vals_val_get(dft_section, "PRINT%MO%OCCUPATION_NUMBERS", l_val=print_occ)
      print_mo = (cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO") /= 0) .AND. &
                 (print_energies .OR. print_occ .OR. print_coeff)
      print_cubes = (cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO_CUBES") /= 0)

      IF (print_mo .OR. print_cubes) THEN

         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_global=nmo)

         ! C^QP_μn = Σ_m C_μm U_mn
         CALL get_static_cohsex_mo_coeff(mo_coeff, u_qp, fm_qp_mo)

         ! MO energies ϵ_n^QP, occupations f_n (maxocc for n <= homo) and coefficients C^QP_μn
         IF (print_mo) THEN
            IF (nspins == 1) THEN
               maxocc = 2.0_dp
            ELSE
               maxocc = 1.0_dp
            END IF
            nelectron = NINT(maxocc*REAL(homo, KIND=dp))
            CALL allocate_mo_set(qp_mo_set, nao=nao, nmo=nmo, nelectron=nelectron, &
                                 n_el_f=REAL(nelectron, KIND=dp), maxocc=maxocc, &
                                 flexible_electron_count=0.0_dp)
            CALL init_mo_set(qp_mo_set, fm_struct=fm_qp_mo%matrix_struct, name="static_cohsex")
            qp_mo_set%eigenvalues(:) = e_qp(:)
            CALL cp_fm_to_fm(fm_qp_mo, qp_mo_set%mo_coeff)

            CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
            IF (nspins == 1) THEN
               CALL write_mo_set_to_output_unit(qp_mo_set, qs_kind_set, particle_set, dft_section, &
                                                4, 0, final_mos=.TRUE.)
            ELSE IF (ispin == 1) THEN
               CALL write_mo_set_to_output_unit(qp_mo_set, qs_kind_set, particle_set, dft_section, &
                                                4, 0, final_mos=.TRUE., spin="ALPHA")
            ELSE
               CALL write_mo_set_to_output_unit(qp_mo_set, qs_kind_set, particle_set, dft_section, &
                                                4, 0, final_mos=.TRUE., spin="BETA")
            END IF
            CALL deallocate_mo_set(qp_mo_set)
         END IF

         ! cube files of ψ_n(r) = Σ_μ C^QP_μn ϕ_μ(r) for n = homo-NHOMO+1, ..., homo+NLUMO
         IF (print_cubes) THEN
            CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%NHOMO", i_val=nhomo)
            CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%NLUMO", i_val=nlumo)
            CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%STRIDE", i_vals=stride)
            IF (nhomo == -1) nhomo = homo
            IF (nlumo == -1) nlumo = nmo - homo
            first_mo = MAX(1, homo - nhomo + 1)
            last_mo = MIN(nmo, homo + nlumo)

            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                            cell=cell, dft_control=dft_control, particle_set=particle_set, &
                            pw_env=pw_env, subsys=subsys)
            CALL qs_subsys_get(subsys, particles=particles)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL auxbas_pw_pool%create_pw(wf_r)
            CALL auxbas_pw_pool%create_pw(wf_g)

            DO ivector = first_mo, last_mo
               CALL calculate_wavefunction(fm_qp_mo, ivector, wf_r, wf_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env)
               WRITE (fname, '(A,I5.5,A,I1)') "static_cohsex_WFN_", ivector, "_", ispin
               unit_nr = cp_print_key_unit_nr(logger, dft_section, "PRINT%MO_CUBES", &
                                              extension=".cube", middle_name=TRIM(fname), &
                                              file_position="REWIND", log_filename=.FALSE.)
               WRITE (title, *) "static COHSEX WAVEFUNCTION ", ivector, " spin ", ispin, &
                  " i.e. HOMO ", ivector - homo
               CALL cp_pw_to_cube(wf_r, unit_nr, title, particles, stride=stride)
               CALL cp_print_key_finished_output(unit_nr, logger, dft_section, "PRINT%MO_CUBES")
            END DO

            CALL auxbas_pw_pool%give_back_pw(wf_r)
            CALL auxbas_pw_pool%give_back_pw(wf_g)
         END IF

         CALL cp_fm_release(fm_qp_mo)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_static_cohsex_mos

! **************************************************************************************************
!> \brief Canonical static COHSEX MO coefficients
!>          C^QP_μn = Σ_m C_μm U_mn
!> \param mo_coeff ...
!> \param u_qp ...
!> \param fm_qp_mo ...
! **************************************************************************************************
   SUBROUTINE get_static_cohsex_mo_coeff(mo_coeff, u_qp, fm_qp_mo)
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: u_qp
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_qp_mo

      INTEGER                                            :: nao, nmo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_u

      CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_global=nmo)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=mo_coeff%matrix_struct, &
                               nrow_global=nmo, ncol_global=nmo)
      CALL cp_fm_create(fm_u, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_submatrix(fm_u, u_qp)
      CALL cp_fm_create(fm_qp_mo, mo_coeff%matrix_struct)
      CALL parallel_gemm(transa="N", transb="N", m=nao, n=nmo, k=nmo, alpha=1.0_dp, &
                         matrix_a=mo_coeff, matrix_b=fm_u, beta=0.0_dp, matrix_c=fm_qp_mo)
      CALL cp_fm_release(fm_u)

   END SUBROUTINE get_static_cohsex_mo_coeff

! **************************************************************************************************
!> \brief Rotate an RI 3-center slab to the canonical static COHSEX MOs
!> \param fm_slab ...
!> \param dimen_RI ...
!> \param n_outer ...
!> \param n_inner ...
!> \param u_outer ...
!> \param u_inner ...
! **************************************************************************************************
   SUBROUTINE rotate_bse_slab(fm_slab, dimen_RI, n_outer, n_inner, u_outer, u_inner)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_slab
      INTEGER, INTENT(IN)                                :: dimen_RI, n_outer, n_inner
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: u_outer, u_inner

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rotate_bse_slab'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! (a) inner index: (P ; o,i) -> (P,o ; i), P fast, reordering [3,1,2,4] and back [3,2,1,4]
      CALL rotate_one_index(fm_slab, dimen_RI, n_outer, n_inner, u_inner, [3, 1, 2, 4], [3, 2, 1, 4], &
                            n_inner)
      ! (b) outer index: (P ; o,i) -> (P,i ; o), P fast, reordering [4,1,2,3] and back [3,2,4,1]
      CALL rotate_one_index(fm_slab, dimen_RI, n_inner, n_outer, u_outer, [4, 1, 2, 3], [3, 2, 4, 1], &
                            n_inner)

      CALL timestop(handle)

   END SUBROUTINE rotate_bse_slab

! **************************************************************************************************
!> \brief rotate one MO index of the compound column index
!> \param fm_slab ...
!> \param dimen_RI ...
!> \param n_spect ...
!> \param n_rot ...
!> \param u_rot ...
!> \param reord_to ...
!> \param reord_back ...
!> \param n_inner size ...
! **************************************************************************************************
   SUBROUTINE rotate_one_index(fm_slab, dimen_RI, n_spect, n_rot, u_rot, reord_to, reord_back, &
                               n_inner)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_slab
      INTEGER, INTENT(IN)                                :: dimen_RI, n_spect, n_rot
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: u_rot
      INTEGER, DIMENSION(4)                              :: reord_to, reord_back
      INTEGER, INTENT(IN)                                :: n_inner

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_slab_out, fm_t, fm_t_rot, fm_u

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=fm_slab%matrix_struct, &
                               nrow_global=n_spect*dimen_RI, ncol_global=n_rot)
      CALL cp_fm_create(fm_t, fm_struct, name="static_cohsex_slab_reordered", set_zero=.TRUE.)
      CALL cp_fm_create(fm_t_rot, fm_struct, name="static_cohsex_slab_rotated", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_permute(fm_t, fm_slab, 1.0_dp, 1, n_inner, dimen_RI, n_rot, reord_to)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, template_fmstruct=fm_slab%matrix_struct, &
                               nrow_global=n_rot, ncol_global=n_rot)
      CALL cp_fm_create(fm_u, fm_struct, name="static_cohsex_U")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_submatrix(fm_u, u_rot)
      CALL parallel_gemm(transa="N", transb="N", m=n_spect*dimen_RI, n=n_rot, k=n_rot, alpha=1.0_dp, &
                         matrix_a=fm_t, matrix_b=fm_u, beta=0.0_dp, matrix_c=fm_t_rot)
      CALL cp_fm_release(fm_u)
      CALL cp_fm_release(fm_t)

      CALL cp_fm_create(fm_slab_out, fm_slab%matrix_struct, name="static_cohsex_slab_back", set_zero=.TRUE.)
      CALL cp_fm_permute(fm_slab_out, fm_t_rot, 1.0_dp, dimen_RI, n_rot, 1, n_inner, reord_back)
      CALL cp_fm_to_fm(fm_slab_out, fm_slab)
      CALL cp_fm_release(fm_slab_out)
      CALL cp_fm_release(fm_t_rot)

   END SUBROUTINE rotate_one_index

END MODULE rpa_gw_static_cohsex

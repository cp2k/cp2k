!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Adaptively Compressed Exchange (ACE) operator for HFX.
!>        Reference: Lin, J. Chem. Theory Comput. 2016, 12, 5, 2242-2249
!>
!>  Algorithm (per spin):
!>
!>    BUILD  (first call, or every rebuild_frequency steps):
!>      1. Full HFX: ks_matrix = K_HFX + H_core, energy%ex = E_x(full)
!>      2. K_AO  = ks_matrix - H_core           (nao x nao, negative semidefinite)
!>      3. C_occ = first nocc columns of mo_coeff, redistributed to a layout
!>                 compatible with K_AO so that PDGEMM works correctly
!>      4. xi    = K_AO * C_occ                 (nao x nocc)
!>      5. M     = C_occ^T * xi                 (nocc x nocc, negative definite)
!>      6. -M    = U^T U  via Cholesky          (U upper triangular, stored in M_fm)
!>      7. W     = xi * U^{-1}                  (nao x nocc, the ACE projector)
!>      8. Apply (see below) to update ks_matrix and energy%ex
!>
!>    APPLY  (all other steps):
!>      ks_matrix = H_core - W * W^T
!>      E_x       = -0.5 * Tr[W^T * P * W]
!>
!>  Important implementation notes:
!>
!>  (a) Distribution compatibility.
!>      mo_coeff is created by the SCF pool machinery with its own block-cyclic
!>      layout; K_AO is freshly created with the QS blacs_env default layout.
!>      ScaLAPACK PDGEMM requires all three matrices (A,B,C) to share the same
!>      process grid; mismatched block sizes cause PDGEMM to silently return zero.
!>      Fix: cp_fm_to_fm redistributes mo_coeff into C_occ_fm (same struct as
!>      K_fm) before any GEMM call.  cp_fm_to_fm also trims nmo → nocc columns.
!>
!>  (b) Zero mo_coeff on the first SCF step.
!>      With an ATOMIC initial guess, qs_scf_initialization::calculate_first_
!>      density_matrix constructs rho_ao from atomic superposition without
!>      diagonalising any Hamiltonian.  mo_coeff is therefore zero until after
!>      the first diagonalisation (standard diag) or ot_scf_mini step (OT).
!>      Fix: detect ||mo_coeff||_F ≈ 0, leave ks_matrix as the full HFX result
!>      from step 1, and signal the caller to retry the BUILD on the next step.
!>
!>  (c) OT dual-storage.
!>      When use_mo_coeff_b = .TRUE. (OT), the live orbitals are in the sparse
!>      mo_coeff_b (dbcsr).  mo_coeff (fm) may be stale.
!>      Fix: copy_dbcsr_to_fm(mo_coeff_b, mo_coeff) before reading mo_coeff.
!>
!>  (d) Forces.
!>      Ionic forces require derivatives of HFX integrals w.r.t. nuclear
!>      positions; the APPLY path K_ACE = -W W^T contains no such derivatives.
!>      When calculate_forces = .TRUE., fall back to full HFX.
!>
!> \author Ritama Kar
! **************************************************************************************************

MODULE hfx_ace_methods

    USE cp_blacs_env,           ONLY: cp_blacs_env_type
    USE cp_control_types,       ONLY: dft_control_type
    USE cp_dbcsr_api,           ONLY: dbcsr_add,               &
                                      dbcsr_copy,              &
                                      dbcsr_create,            &
                                      dbcsr_p_type,            &
                                      dbcsr_release,           &
                                      dbcsr_set,               &
                                      dbcsr_type
    USE cp_dbcsr_operations,    ONLY: copy_dbcsr_to_fm,        &
                                      cp_dbcsr_plus_fm_fm_t
    USE cp_fm_basic_linalg,     ONLY: cp_fm_column_scale,      &
                                      cp_fm_gemm,              &
                                      cp_fm_scale,             &
                                      cp_fm_trace,             &
                                      cp_fm_triangular_multiply
    USE cp_fm_cholesky,         ONLY: cp_fm_cholesky_decompose
    USE cp_fm_struct,           ONLY: cp_fm_struct_create,     &
                                      cp_fm_struct_release,    &
                                      cp_fm_struct_type
    USE cp_fm_types,            ONLY: cp_fm_create,            &
                                      cp_fm_get_info,          &
                                      cp_fm_release,           &
                                      cp_fm_to_fm,             &
                                      cp_fm_type
    USE cp_log_handling,        ONLY: cp_get_default_logger,   &
                                      cp_logger_get_default_io_unit, &
                                      cp_logger_type
    USE hfx_admm_utils,         ONLY: hfx_ks_matrix
    USE hfx_types,              ONLY: hfx_type
    USE input_section_types,    ONLY: section_vals_type
    USE kinds,                  ONLY: dp
    USE message_passing,        ONLY: mp_para_env_type
    USE pw_types,               ONLY: pw_r3d_rs_type
    USE qs_energy_types,        ONLY: qs_energy_type
    USE qs_environment_types,   ONLY: get_qs_env,              &
                                      qs_environment_type
    USE qs_mo_types,            ONLY: get_mo_set,              &
                                      mo_set_type
    USE qs_rho_types,           ONLY: qs_rho_get,              &
                                      qs_rho_type
#include "./base/base_uses.f90"
 
    IMPLICIT NONE
    PRIVATE
 
    CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ace_methods'
 
    ! -----------------------------------------------------------------------
    ! Module-level state: persists across SCF steps within one run.
    !
    !   ace_W(1,ispin)    ACE projector W, shape nao x nocc
    !   ace_is_built      .FALSE. until at least one successful BUILD
    !   ace_step_counter  counts calls since last BUILD; triggers rebuild at
    !                     every rebuild_frequency steps
    ! -----------------------------------------------------------------------
    TYPE(cp_fm_type), ALLOCATABLE, SAVE :: ace_W(:,:)
    LOGICAL,                       SAVE :: ace_is_built     = .FALSE.
    INTEGER,                       SAVE :: ace_step_counter = 0
 
    ! -----------------------------------------------------------------------
    ! Debug flags  (set all .FALSE. for production)
    !
    !   DBG_ROUTING  BUILD / APPLY / DEFER decision, step counters, exit energy
    !   DBG_BUILD    E_x(full), orbital info, matrix norms during BUILD
    !   DBG_ENERGY   E_x(ACE); on BUILD steps also |E_x(ACE) - E_x(full)|
    ! -----------------------------------------------------------------------
    LOGICAL, PARAMETER, PRIVATE :: DBG_ROUTING = .FALSE.
    LOGICAL, PARAMETER, PRIVATE :: DBG_BUILD   = .FALSE.
    LOGICAL, PARAMETER, PRIVATE :: DBG_ENERGY  = .TRUE.
 
    PUBLIC :: hfx_ace_ks_matrix, hfx_ace_release
 
 CONTAINS
 
 ! **************************************************************************************************
 !> \brief Main ACE entry point, replacing hfx_ks_matrix in qs_ks_methods.
 ! **************************************************************************************************
    SUBROUTINE hfx_ace_ks_matrix(qs_env, ks_matrix, rho, energy,          &
                                  calculate_forces, just_energy,            &
                                  v_rspace_new, v_tau_rspace,               &
                                  ace_rebuild_frequency, ext_xc_section)
 
       TYPE(qs_environment_type),  POINTER                :: qs_env
       TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: ks_matrix
       TYPE(qs_rho_type),          POINTER                :: rho
       TYPE(qs_energy_type),       POINTER                :: energy
       LOGICAL,            INTENT(IN)                     :: calculate_forces, just_energy
       TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_tau_rspace
       INTEGER,            INTENT(IN)                     :: ace_rebuild_frequency
       TYPE(section_vals_type), OPTIONAL, POINTER         :: ext_xc_section
 
       CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ace_ks_matrix'
 
       INTEGER                              :: handle, iw, n_rep_hf, nspins, rebuild_freq
       LOGICAL                              :: rebuild_ace, ace_built_now
       TYPE(cp_logger_type),       POINTER  :: logger
       TYPE(dft_control_type),     POINTER  :: dft_control
       TYPE(hfx_type), DIMENSION(:,:), POINTER :: x_data
 
       CALL timeset(routineN, handle)
       NULLIFY(logger, dft_control, x_data)
 
       logger   => cp_get_default_logger()
       iw       =  cp_logger_get_default_io_unit(logger)
 
       CALL get_qs_env(qs_env, x_data=x_data, dft_control=dft_control)
       n_rep_hf = SIZE(x_data, 1)
       nspins   = dft_control%nspins
 
       ! First-implementation restrictions
       IF (n_rep_hf /= 1) CPABORT("ACE: only one &HF section is supported.")
       IF (dft_control%nimages /= 1) &
          CPABORT("ACE: k-points / multiple images are not implemented.")
       IF (dft_control%do_admm) &
          CPABORT("ACE: ADMM is not implemented. Remove &ACE or disable ADMM.")
 
       rebuild_freq = MAX(1, ace_rebuild_frequency)
 
       ! ------------------------------------------------------------------
       ! Bypass A: energy-only call (ks_matrix is not needed by the caller)
       ! ------------------------------------------------------------------
       IF (just_energy) THEN
          IF (DBG_ROUTING .AND. iw > 0) &
             WRITE(iw,'(T2,A)') 'ACE | just_energy=T: full HFX (no matrix update)'
          CALL hfx_call(qs_env, ks_matrix, rho, energy,                  &
                         calculate_forces, just_energy,                    &
                         v_rspace_new, v_tau_rspace, ext_xc_section)
          CALL timestop(handle)
          RETURN
       END IF
 
       ! ------------------------------------------------------------------
       ! Bypass B: ionic forces requested.
       ! The ACE apply path (K_ACE = -W W^T) contains no ionic force
       ! contributions (no derivatives of HFX integrals w.r.t. R_I).
       ! Use full HFX so forces are exact.
       ! ------------------------------------------------------------------
       IF (calculate_forces) THEN
          IF (DBG_ROUTING .AND. iw > 0) &
             WRITE(iw,'(T2,A)') 'ACE | calculate_forces=T: full HFX for exact forces'
          CALL hfx_call(qs_env, ks_matrix, rho, energy,                  &
                         calculate_forces, just_energy,                    &
                         v_rspace_new, v_tau_rspace, ext_xc_section)
          CALL timestop(handle)
          RETURN
       END IF
 
       ! ------------------------------------------------------------------
       ! Rebuild decision
       ! ------------------------------------------------------------------
       rebuild_ace = (.NOT. ace_is_built) .OR. &
                     (MOD(ace_step_counter, rebuild_freq) == 0)
 
       IF (DBG_ROUTING .AND. iw > 0) THEN
          WRITE(iw,'(/,T2,A)')    REPEAT('-', 56)
          WRITE(iw,'(T2,A)')      'ACE | hfx_ace_ks_matrix'
          WRITE(iw,'(T4,A,L1)')   'ace_is_built  = ', ace_is_built
          WRITE(iw,'(T4,A,L1)')   'rebuild_ace   = ', rebuild_ace
          WRITE(iw,'(T4,A,I6)')   'step_counter  = ', ace_step_counter
          WRITE(iw,'(T4,A,I6)')   'rebuild_freq  = ', rebuild_freq
          WRITE(iw,'(T4,A,I4)')   'nspins        = ', nspins
          WRITE(iw,'(T4,A)')      MERGE('-> BUILD', '-> APPLY', rebuild_ace)
          WRITE(iw,'(T2,A)')      REPEAT('-', 56)
       END IF
 
       IF (rebuild_ace) THEN
 
          ace_built_now = .FALSE.
          CALL hfx_ace_build_projector(qs_env, ks_matrix, rho, energy,   &
                                        just_energy,                       &
                                        v_rspace_new, v_tau_rspace,        &
                                        nspins, iw, ace_built_now,         &
                                        ext_xc_section)
          IF (ace_built_now) THEN
             ace_is_built     = .TRUE.
             ace_step_counter = 1
             IF (DBG_ROUTING .AND. iw > 0) &
                WRITE(iw,'(T4,A)') 'ACE | W built. Projector live from next step.'
          ELSE
             ! mo_coeff was zero: full HFX (from step 1 of the build) is
             ! already in ks_matrix.  Keep ace_is_built=.FALSE. so we retry.
             ace_is_built     = .FALSE.
             ace_step_counter = 0
             IF (DBG_ROUTING .AND. iw > 0) &
                WRITE(iw,'(T4,A)') 'ACE | Build deferred (C_occ=0). Full HFX in ks_matrix.'
          END IF
 
       ELSE
 
          CALL hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy, nspins, iw)
          ace_step_counter = ace_step_counter + 1
 
       END IF
 
       IF (DBG_ROUTING .AND. iw > 0) THEN
          WRITE(iw,'(T4,A,F20.10)') 'energy%ex on exit = ', energy%ex
          WRITE(iw,'(T4,A,I6)')     'step_counter now  = ', ace_step_counter
       END IF
 
       CALL timestop(handle)
 
    END SUBROUTINE hfx_ace_ks_matrix
 
 ! **************************************************************************************************
 !> \brief Build the ACE projector W.  Sets build_succeeded=.FALSE. and
 !>        returns early (full HFX result in ks_matrix) when mo_coeff is zero.
 ! **************************************************************************************************
    SUBROUTINE hfx_ace_build_projector(qs_env, ks_matrix, rho, energy,    &
                                        just_energy,                        &
                                        v_rspace_new, v_tau_rspace,         &
                                        nspins, iw, build_succeeded,        &
                                        ext_xc_section)
 
       TYPE(qs_environment_type),  POINTER                :: qs_env
       TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: ks_matrix
       TYPE(qs_rho_type),          POINTER                :: rho
       TYPE(qs_energy_type),       POINTER                :: energy
       LOGICAL,            INTENT(IN)                     :: just_energy
       TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_tau_rspace
       INTEGER,            INTENT(IN)                     :: nspins, iw
       LOGICAL,            INTENT(OUT)                    :: build_succeeded
       TYPE(section_vals_type), OPTIONAL, POINTER         :: ext_xc_section
 
       CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ace_build_projector'
 
       INTEGER                              :: handle, info_chol, ispin, nao, nmo, nocc
       REAL(dp)                             :: ehfx_full, frob
       TYPE(cp_blacs_env_type),    POINTER  :: blacs_env
       TYPE(cp_fm_struct_type),    POINTER  :: fmstruct
       TYPE(cp_fm_type)                     :: K_fm, C_occ_fm, xi_fm, M_fm
       TYPE(cp_fm_type),           POINTER  :: mo_coeff
       TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER :: matrix_h, rho_ao
       TYPE(dbcsr_type)                     :: K_ao_dbcsr
       TYPE(mo_set_type), DIMENSION(:), POINTER    :: mos
       TYPE(mp_para_env_type),     POINTER  :: para_env
       REAL(dp), DIMENSION(:),     POINTER  :: occ_nums
 
       CALL timeset(routineN, handle)
       NULLIFY(blacs_env, para_env, mos, mo_coeff, matrix_h, occ_nums, fmstruct, rho_ao)
 
       build_succeeded = .FALSE.
 
       CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env,   &
                       mos=mos, matrix_h_kp=matrix_h)
 
       ! ----------------------------------------------------------------
       ! Step 1: full HFX.
       !   After return: ks_matrix(ispin,1) = K_HFX(ispin) + H_core
       !                 energy%ex           = E_x(full HFX)
       ! forces=.FALSE.: the caller already handled the forces bypass.
       ! ----------------------------------------------------------------
       CALL hfx_call(qs_env, ks_matrix, rho, energy,                     &
                     .FALSE., just_energy,                                 &
                     v_rspace_new, v_tau_rspace, ext_xc_section)
       ehfx_full = energy%ex

       CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
 
       IF (DBG_BUILD .AND. iw > 0) &
          WRITE(iw,'(/,T2,A,F20.10)') 'ACE BUILD | E_x(full HFX) = ', ehfx_full
 
       ! Allocate / resize module storage
       IF (ALLOCATED(ace_W)) THEN
          IF (SIZE(ace_W, 2) /= nspins) CALL hfx_ace_release()
       END IF
       IF (.NOT. ALLOCATED(ace_W)) ALLOCATE(ace_W(1, nspins))
 
       ! ----------------------------------------------------------------
       ! Per-spin build loop
       ! ----------------------------------------------------------------
       DO ispin = 1, nspins
 
          ! --------------------------------------------------------------
          ! Sync mo_coeff (fm) from mo_coeff_b (dbcsr) when needed.
          !
          ! OT uses use_mo_coeff_b=.TRUE.: the live orbitals are stored in
          ! the sparse mo_coeff_b after each ot_scf_mini step; mo_coeff
          ! (fm) may be stale.  copy_dbcsr_to_fm makes mo_coeff current.
          !
          ! Standard diagonalisation uses use_mo_coeff_b=.FALSE.: mo_coeff
          ! (fm) is updated directly after each Hamiltonian diagonalisation.
          ! No action needed.
          ! --------------------------------------------------------------
          IF (mos(ispin)%use_mo_coeff_b) &
             CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff)
 
          CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, nao=nao, nmo=nmo, &
                          homo=nocc, occupation_numbers=occ_nums)
 
          IF (nocc <= 0) CPABORT("ACE: homo <= 0.")
          IF (nocc > nmo) CPABORT("ACE: homo > nmo.")
          CPASSERT(ASSOCIATED(mo_coeff))
 
          ! --------------------------------------------------------------
          ! Guard: mo_coeff may be zero at the first SCF call.
          !
          ! qs_scf_initialization::calculate_first_density_matrix builds
          ! rho_ao from atomic superposition without diagonalising, so
          ! mo_coeff is never populated on that call.  It becomes non-zero
          ! only after the first diagonalisation / ot_scf_mini.
          ! If zero here: full HFX from step 1 is already the correct
          ! ks_matrix for this step.  Signal the caller to retry BUILD.
          ! --------------------------------------------------------------
          CALL cp_fm_trace(mo_coeff, mo_coeff, frob)  ! frob = ||mo_coeff||_F^2
 
          IF (DBG_BUILD .AND. iw > 0) THEN
             WRITE(iw,'(/,T2,A,I4)')   'ACE BUILD | ispin = ', ispin
             WRITE(iw,'(T4,A,I8)')     'nao  = ', nao
             WRITE(iw,'(T4,A,I8)')     'nmo  = ', nmo
             WRITE(iw,'(T4,A,I8)')     'nocc = ', nocc
             WRITE(iw,'(T4,A,L1)')     'use_mo_coeff_b = ', mos(ispin)%use_mo_coeff_b
             WRITE(iw,'(T4,A,ES12.4)') '||mo_coeff||_F = ', SQRT(MAX(frob, 0.0_dp))
          END IF
 
          IF (frob < 1.0e-20_dp) THEN
             IF (DBG_BUILD .AND. iw > 0) &
                WRITE(iw,'(T4,A)') 'mo_coeff=0: build deferred to next step.'
             CALL timestop(handle)
             RETURN    ! build_succeeded stays .FALSE.
          END IF
 
          ! --------------------------------------------------------------
          ! Step 2: K_AO = ks_matrix - H_core  (sparse, then dense nao x nao)
          !
          ! ks_matrix(ispin,1) = K_HFX + H_core  (from step 1)
          ! K_AO = K_HFX, the exchange contribution (negative semidefinite).
          ! --------------------------------------------------------------
          CALL dbcsr_create(K_ao_dbcsr, template=ks_matrix(ispin,1)%matrix, name="K_AO")
          CALL dbcsr_copy(K_ao_dbcsr, ks_matrix(ispin,1)%matrix)
          CALL dbcsr_add(K_ao_dbcsr, matrix_h(1,1)%matrix, 1.0_dp, -1.0_dp)
 
          NULLIFY(fmstruct)
          CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                    nrow_global=nao, ncol_global=nao)
          CALL cp_fm_create(K_fm, fmstruct, name="K_AO_dense")
          CALL cp_fm_struct_release(fmstruct)
 
          CALL copy_dbcsr_to_fm(K_ao_dbcsr, K_fm)
          CALL dbcsr_release(K_ao_dbcsr)
 
          IF (DBG_BUILD .AND. iw > 0) THEN
             CALL cp_fm_trace(K_fm, K_fm, frob)
             WRITE(iw,'(T4,A,ES12.4)') '||K_AO||_F = ', SQRT(MAX(frob, 0.0_dp))
          END IF
 
          ! --------------------------------------------------------------
          ! Step 3: C_occ (nao x nocc) — redistribute mo_coeff.
          !
          ! mo_coeff (nao x nmo) was created by the SCF ao_mo pool with its
          ! own block-cyclic layout; K_fm was created with blacs_env defaults.
          ! ScaLAPACK PDGEMM silently returns zero if the three operands use
          ! incompatible distributions (see module header note (a)).
          !
          ! cp_fm_to_fm redistributes between arbitrary cp_fm_struct via BLACS
          ! and copies only the first nocc columns (occupied MOs), trimming
          ! the nmo-wide source down to the nocc-wide target automatically.
          !
          ! C_occ_fm and xi_fm share the same fmstruct (nao x nocc), so all
          ! subsequent GEMMs and the trace in APPLY are distribution-safe.
          ! --------------------------------------------------------------
          NULLIFY(fmstruct)
          CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                    nrow_global=nao, ncol_global=nocc)
          CALL cp_fm_create(C_occ_fm, fmstruct, name="C_occ")
          CALL cp_fm_create(xi_fm,    fmstruct, name="xi")
          CALL cp_fm_struct_release(fmstruct)
 
          CALL cp_fm_to_fm(mo_coeff, C_occ_fm)   ! redistribute + copy first nocc cols
          
           ! --------------------------------------------------------------
            ! Diagnostic: compare current rho_ao with density reconstructed
            ! from the same C_occ used by ACE.
            !
            ! P_mo = C_occ * Occ * C_occ^T
            !
            ! We avoid forming diff = P_rho - P_mo explicitly.
            ! Instead:
            !
            !        ||P_rho - P_mo||_F^2 =
            !           ||P_rho||_F^2 + ||P_mo||_F^2 - 2 Tr(P_rho^T P_mo)
            !        --------------------------------------------------------------
            BLOCK
               TYPE(cp_fm_type) :: P_rho_fm, P_mo_fm, C_scaled_fm
               TYPE(cp_fm_struct_type), POINTER :: dens_fmstruct
               REAL(dp), ALLOCATABLE :: occ_sqrt(:)
               REAL(dp) :: norm_rho_sq, norm_mo_sq, cross_trace
               REAL(dp) :: norm_diff_sq, norm_rho, norm_mo, norm_diff
               INTEGER :: iocc
            
               NULLIFY(dens_fmstruct)
            
               CALL cp_fm_struct_create(dens_fmstruct, context=blacs_env, para_env=para_env, &
                                        nrow_global=nao, ncol_global=nao)
            
               CALL cp_fm_create(P_rho_fm, dens_fmstruct, name="P_rho_current")
               CALL cp_fm_create(P_mo_fm,  dens_fmstruct, name="P_from_Cocc")
            
               CALL cp_fm_struct_release(dens_fmstruct)
            
               ! Current CP2K density matrix used by HFX
               CALL copy_dbcsr_to_fm(rho_ao(ispin,1)%matrix, P_rho_fm)
            
               ! Build occupation-weighted C_occ
               ALLOCATE(occ_sqrt(nocc))
               DO iocc = 1, nocc
                  occ_sqrt(iocc) = SQRT(MAX(occ_nums(iocc), 0.0_dp))
               END DO
           
               CALL cp_fm_create(C_scaled_fm, C_occ_fm%matrix_struct, name="C_occ_scaled")
               CALL cp_fm_to_fm(C_occ_fm, C_scaled_fm)
               CALL cp_fm_column_scale(C_scaled_fm, occ_sqrt)
           
               ! P_mo = C_scaled * C_scaled^T
               CALL cp_fm_gemm('N', 'T', nao, nao, nocc, &
                               1.0_dp, C_scaled_fm, C_scaled_fm, 0.0_dp, P_mo_fm)
           
               ! Norms and cross term
               CALL cp_fm_trace(P_rho_fm, P_rho_fm, norm_rho_sq)
               CALL cp_fm_trace(P_mo_fm,  P_mo_fm,  norm_mo_sq)
               CALL cp_fm_trace(P_rho_fm, P_mo_fm,  cross_trace)
           
               norm_diff_sq = norm_rho_sq + norm_mo_sq - 2.0_dp * cross_trace
               norm_diff_sq = MAX(norm_diff_sq, 0.0_dp)
           
               norm_rho  = SQRT(MAX(norm_rho_sq,  0.0_dp))
               norm_mo   = SQRT(MAX(norm_mo_sq,   0.0_dp))
               norm_diff = SQRT(norm_diff_sq)
           
               IF (iw > 0) THEN
                  WRITE(iw,'(T4,A,ES14.6)') 'ACE DIAG ||rho_ao||_F        = ', norm_rho
                  WRITE(iw,'(T4,A,ES14.6)') 'ACE DIAG ||P(C_occ)||_F      = ', norm_mo
                  WRITE(iw,'(T4,A,ES14.6)') 'ACE DIAG Tr(rho_ao*P(C))    = ', cross_trace
                  WRITE(iw,'(T4,A,ES14.6)') 'ACE DIAG ||rho_ao-P(C)||_F   = ', norm_diff
                  WRITE(iw,'(T4,A,ES14.6)') 'ACE DIAG relative difference = ', &
                     norm_diff / MAX(norm_rho, 1.0E-14_dp)
                  WRITE(iw,'(T4,A,100F10.5)') 'ACE DIAG occ_nums = ', occ_nums(1:nocc)
               END IF
           
               CALL cp_fm_release(P_rho_fm)
               CALL cp_fm_release(P_mo_fm)
               CALL cp_fm_release(C_scaled_fm)
               DEALLOCATE(occ_sqrt)
            END BLOCK
          ! --------------------------------------------------------------
          ! Step 4: xi = K_AO * C_occ  (nao x nocc)
          ! --------------------------------------------------------------
          CALL cp_fm_gemm('N', 'N', nao, nocc, nao, &
                           1.0_dp, K_fm, C_occ_fm, 0.0_dp, xi_fm)
          CALL cp_fm_release(K_fm)
 
          IF (DBG_BUILD .AND. iw > 0) THEN
             CALL cp_fm_trace(xi_fm, xi_fm, frob)
             WRITE(iw,'(T4,A,ES12.4)') '||xi||_F = ', SQRT(MAX(frob, 0.0_dp))
          END IF
 
          ! --------------------------------------------------------------
          ! Step 5: M = C_occ^T * xi  (nocc x nocc)
          !
          ! M is symmetric and negative definite for a non-degenerate system
          ! (K_AO negative semidefinite, C_occ full column rank).
          ! --------------------------------------------------------------
          NULLIFY(fmstruct)
          CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                    nrow_global=nocc, ncol_global=nocc)
          CALL cp_fm_create(M_fm, fmstruct, name="M")
          CALL cp_fm_struct_release(fmstruct)
 
          CALL cp_fm_gemm('T', 'N', nocc, nocc, nao, &
                           1.0_dp, C_occ_fm, xi_fm, 0.0_dp, M_fm)
          CALL cp_fm_release(C_occ_fm)   ! no longer needed
 
          IF (DBG_BUILD .AND. iw > 0) THEN
             CALL cp_fm_trace(M_fm, M_fm, frob)
             WRITE(iw,'(T4,A,ES12.4)') '||M||_F = ', SQRT(MAX(frob, 0.0_dp))
          END IF
 
          ! --------------------------------------------------------------
          ! Step 6: Cholesky of -M = U^T U  (U upper triangular)
          !
          ! cp_fm_cholesky_decompose calls PDPOTRF('U',...).
          ! On return M_fm holds U in its upper triangle.
          ! -M must be positive definite; if not, K_AO or C_occ is wrong.
          ! --------------------------------------------------------------
          CALL cp_fm_scale(-1.0_dp, M_fm)
          CALL cp_fm_cholesky_decompose(M_fm, n=nocc, info_out=info_chol)
 
          IF (info_chol /= 0) THEN
             IF (iw > 0) THEN
                WRITE(iw,'(T4,A,I6)')   'ACE | Cholesky failed, info = ', info_chol
                WRITE(iw,'(T4,A,F20.10)') 'ACE | E_x(full) = ', ehfx_full
                WRITE(iw,'(T4,A,I8,A,I8)') 'ACE | nao=', nao, '  nocc=', nocc
             END IF
             CPABORT("ACE: Cholesky of -M failed (not positive definite).")
          END IF
 
          IF (DBG_BUILD .AND. iw > 0) &
             WRITE(iw,'(T4,A)') 'Cholesky OK (info=0).'
 
          ! --------------------------------------------------------------
          ! Step 7: W = xi * U^{-1}  (nao x nocc)
          !
          ! cp_fm_triangular_multiply computes  W := W * U^{-1}  in-place:
          !   side='R'          right-multiply by the triangular matrix
          !   uplo_tr='U'       U is stored in the upper triangle of M_fm
          !   invert_tr=.TRUE.  use U^{-1} rather than U
          !   transpose_tr=.FALSE.  no additional transpose
          ! --------------------------------------------------------------
          IF (ASSOCIATED(ace_W(1,ispin)%matrix_struct)) &
             CALL cp_fm_release(ace_W(1,ispin))
 
          CALL cp_fm_create(ace_W(1,ispin), xi_fm%matrix_struct, name="W_ACE")
          CALL cp_fm_to_fm(xi_fm, ace_W(1,ispin))   ! copy xi → W (same struct, no redistrib)
 
          CALL cp_fm_triangular_multiply(M_fm, ace_W(1,ispin),            &
                                          side='R', uplo_tr='U',           &
                                          transpose_tr=.FALSE.,            &
                                          invert_tr=.TRUE.,                &
                                          n_rows=nao, n_cols=nocc,         &
                                          alpha=1.0_dp)
 
          CALL cp_fm_release(xi_fm)
          CALL cp_fm_release(M_fm)
 
          IF (DBG_BUILD .AND. iw > 0) THEN
             CALL cp_fm_trace(ace_W(1,ispin), ace_W(1,ispin), frob)
             WRITE(iw,'(T4,A,I4,A,2I8,A,ES12.4)') &
                'W spin=', ispin, ' shape=', nao, nocc, &
                '  ||W||_F=', SQRT(MAX(frob, 0.0_dp))
          END IF
 
       END DO  ! ispin
 
       build_succeeded = .TRUE.
 
       ! ----------------------------------------------------------------
       ! Step 8: apply immediately.
       ! Sets ks_matrix = H_core - W W^T and energy%ex = E_x(ACE).
       ! Passing ehfx_full_ref triggers the quality-check print.
       ! ----------------------------------------------------------------
       CALL hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy,       &
                                     nspins, iw, ehfx_full_ref=ehfx_full)
 
       CALL timestop(handle)
 
    END SUBROUTINE hfx_ace_build_projector
 
 ! **************************************************************************************************
 !> \brief Apply the stored ACE projector.
 !>
 !>  For each spin:
 !>    ks_matrix(ispin) = H_core - W(ispin) * W(ispin)^T
 !>    E_x += -0.5 * Tr[ W(ispin)^T * P(ispin) * W(ispin) ]
 !>
 !>  The optional argument ehfx_full_ref is used only for the quality-check
 !>  print on BUILD steps (compare E_x(ACE) to E_x(full HFX)).
 ! **************************************************************************************************
    SUBROUTINE hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy,    &
                                        nspins, iw, ehfx_full_ref)
 
       TYPE(qs_environment_type),  POINTER                :: qs_env
       TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: ks_matrix
       TYPE(qs_rho_type),          POINTER                :: rho
       TYPE(qs_energy_type),       POINTER                :: energy
       INTEGER,            INTENT(IN)                     :: nspins, iw
       REAL(dp), OPTIONAL, INTENT(IN)                     :: ehfx_full_ref
 
       CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ace_apply_projector'
 
       INTEGER                              :: handle, ispin, nao, nocc
       REAL(dp)                             :: ehfx_ace, trace_val
       TYPE(cp_blacs_env_type),    POINTER  :: blacs_env
       TYPE(cp_fm_struct_type),    POINTER  :: fmstruct
       TYPE(cp_fm_type)                     :: P_fm, PW_fm
       TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER :: matrix_h, rho_ao
       TYPE(mp_para_env_type),     POINTER  :: para_env
 
       CALL timeset(routineN, handle)
       NULLIFY(blacs_env, para_env, matrix_h, rho_ao, fmstruct)
       CPASSERT(ALLOCATED(ace_W))
 
       CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env,   &
                       matrix_h_kp=matrix_h)
       CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
 
       ! ----------------------------------------------------------------
       ! KS matrix: ks_matrix(ispin) = H_core - W * W^T
       !
       ! cp_dbcsr_plus_fm_fm_t adds  alpha * V * V^T  to a sparse matrix.
       ! With alpha=-1, V=ace_W: sparse_matrix -= W * W^T
       ! keep_sparsity=.TRUE.: only fill elements in the existing pattern.
       ! ----------------------------------------------------------------
       DO ispin = 1, nspins
          CALL cp_fm_get_info(ace_W(1,ispin), nrow_global=nao, ncol_global=nocc)
          CALL dbcsr_set(ks_matrix(ispin,1)%matrix, 0.0_dp)
          CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=ks_matrix(ispin,1)%matrix, &
                                      matrix_v=ace_W(1,ispin), ncol=nocc,      &
                                      alpha=-1.0_dp, keep_sparsity=.TRUE.)
          CALL dbcsr_add(ks_matrix(ispin,1)%matrix, matrix_h(1,1)%matrix,     &
                         1.0_dp, 1.0_dp)
       END DO
 
       ! ----------------------------------------------------------------
       ! Exchange energy: E_x = -0.5 * sum_spin Tr[ W^T * P * W ]
       !
       ! PW = P * W (nao x nocc), then Tr[W^T * PW] = cp_fm_trace(W, PW).
       !
       ! PW_fm is created with ace_W's own matrix_struct so that the
       ! cp_fm_trace call (which requires identical distributions) is safe.
       ! ----------------------------------------------------------------
       ehfx_ace = 0.0_dp
       DO ispin = 1, nspins
          CALL cp_fm_get_info(ace_W(1,ispin), nrow_global=nao, ncol_global=nocc)
 
          NULLIFY(fmstruct)
          CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                    nrow_global=nao, ncol_global=nao)
          CALL cp_fm_create(P_fm, fmstruct, name="P_dense")
          CALL cp_fm_struct_release(fmstruct)
          CALL copy_dbcsr_to_fm(rho_ao(ispin,1)%matrix, P_fm)
 
          ! PW_fm shares the struct of ace_W: same distribution for cp_fm_trace
          CALL cp_fm_create(PW_fm, ace_W(1,ispin)%matrix_struct, name="PW")
          CALL cp_fm_gemm('N', 'N', nao, nocc, nao,                        &
                           1.0_dp, P_fm, ace_W(1,ispin), 0.0_dp, PW_fm)
          CALL cp_fm_trace(ace_W(1,ispin), PW_fm, trace_val)
 
          ehfx_ace = ehfx_ace - 0.5_dp * trace_val
 
          CALL cp_fm_release(P_fm)
          CALL cp_fm_release(PW_fm)
 
          IF (DBG_ENERGY .AND. iw > 0) &
             WRITE(iw,'(T4,A,I4,A,F20.10)') &
                'ispin=', ispin, '  E_x(ACE) += ', -0.5_dp * trace_val
       END DO
 
       energy%ex = ehfx_ace
 
       IF (DBG_ENERGY .AND. iw > 0) THEN
          WRITE(iw,'(T2,A,F20.10)') 'ACE | E_x(ACE)  = ', ehfx_ace
          IF (PRESENT(ehfx_full_ref)) THEN
             WRITE(iw,'(T2,A,F20.10)') 'ACE | E_x(full) = ', ehfx_full_ref
             WRITE(iw,'(T2,A,ES12.4)') 'ACE | |delta|   = ', ABS(ehfx_ace - ehfx_full_ref)
             WRITE(iw,'(T2,A)') '(|delta| should be ~0 on BUILD steps; small is good)'
          END IF
       END IF
 
       CALL timestop(handle)
 
    END SUBROUTINE hfx_ace_apply_projector
 
 ! **************************************************************************************************
 !> \brief Release all ACE storage and reset state flags.
 ! **************************************************************************************************
    SUBROUTINE hfx_ace_release(iw_opt)
 
       INTEGER, OPTIONAL, INTENT(IN) :: iw_opt
 
       INTEGER :: i, j, iw
 
       iw = -1
       IF (PRESENT(iw_opt)) iw = iw_opt
 
       IF (ALLOCATED(ace_W)) THEN
          DO j = 1, SIZE(ace_W, 2)
             DO i = 1, SIZE(ace_W, 1)
                IF (ASSOCIATED(ace_W(i,j)%matrix_struct)) CALL cp_fm_release(ace_W(i,j))
             END DO
          END DO
          DEALLOCATE(ace_W)
       END IF
 
       ace_is_built     = .FALSE.
       ace_step_counter = 0
 
       IF (iw > 0) WRITE(iw,'(T2,A)') 'ACE | storage released, counters reset'
 
    END SUBROUTINE hfx_ace_release
 
 ! **************************************************************************************************
 !> \brief Private helper: call hfx_ks_matrix with or without ext_xc_section.
 !>        Avoids duplicating the PRESENT branch everywhere.
 ! **************************************************************************************************
    SUBROUTINE hfx_call(qs_env, ks_matrix, rho, energy,                   &
                         calculate_forces, just_energy,                     &
                         v_rspace_new, v_tau_rspace, ext_xc_section)
 
       TYPE(qs_environment_type),  POINTER                :: qs_env
       TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: ks_matrix
       TYPE(qs_rho_type),          POINTER                :: rho
       TYPE(qs_energy_type),       POINTER                :: energy
       LOGICAL,            INTENT(IN)                     :: calculate_forces, just_energy
       TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_tau_rspace
       TYPE(section_vals_type), OPTIONAL, POINTER         :: ext_xc_section
 
       IF (PRESENT(ext_xc_section)) THEN
          CALL hfx_ks_matrix(qs_env, ks_matrix, rho, energy,              &
                              calculate_forces, just_energy,                &
                              v_rspace_new, v_tau_rspace,                   &
                              ext_xc_section=ext_xc_section)
       ELSE
          CALL hfx_ks_matrix(qs_env, ks_matrix, rho, energy,              &
                              calculate_forces, just_energy,                &
                              v_rspace_new, v_tau_rspace)
       END IF
 
    END SUBROUTINE hfx_call
 
 END MODULE hfx_ace_methods
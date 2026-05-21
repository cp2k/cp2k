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
!>  Diagnostics (controlled by DBG_STALE / DBG_EXACT_EX module flags):
!>
!>    DIAG A  projector staleness  (cheap, runs every APPLY step)
!>      Computes ||W^T C_occ^current||_F / ||W^T C_occ^BUILD||_F.
!>      Ratio = 1 → W still accurate.  Ratio -> 0 → W is stale.
!>
!>    DIAG B  exact vs ACE exchange energy  (expensive: one full HFX per APPLY)
!>      Calls full HFX with just_energy=.TRUE. to get E_x^exact[P^k] and
!>      compares to E_x^ACE[P^k].  Growing |delta| confirms stale W.
!>      ACE ks_matrix and energy%ex are restored after the diagnostic.
!>
!> \author Ritama Kar
! **************************************************************************************************

MODULE hfx_ace_methods

   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE bibliography,                    ONLY: Lin2016ACE,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_plus_fm_fm_t
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                              cp_fm_scale,&
                                              cp_fm_trace,&
                                              cp_fm_triangular_multiply
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE hfx_admm_utils,                  ONLY: hfx_ks_matrix
   USE hfx_types,                       ONLY: hfx_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ace_methods'

   ! -----------------------------------------------------------------------
   ! Module-level state: persists across SCF steps within one run.
   !
   !   ace_W(1,ispin)   ACE projector W, shape nao x nocc
   !   ace_is_built     .FALSE. until at least one successful BUILD
   !   ace_step_counter counts calls since last BUILD
   !   ace_W_ref_norm   ||W^T C_occ^BUILD||_F stored at BUILD for DIAG A
   ! -----------------------------------------------------------------------
   TYPE(cp_fm_type), ALLOCATABLE, SAVE :: ace_W(:, :)
   LOGICAL, SAVE :: ace_is_built = .FALSE.
   INTEGER, SAVE :: ace_step_counter = 0
   REAL(dp), SAVE :: ace_W_ref_norm = 0.0_dp
   INTEGER, SAVE :: ace_geo_step = 0   ! NEW

   ! -----------------------------------------------------------------------
   ! Debug / diagnostic flags  — set all .FALSE. for production.
   !
   !   DBG_ROUTING   BUILD/APPLY/DEFER decisions, counters
   !   DBG_BUILD     norms during BUILD
   !   DBG_ENERGY    E_x(ACE) at every step; DIAG C on BUILD steps
   !   DBG_STALE     DIAG A: staleness ratio at every APPLY step (cheap)
   !   DBG_EXACT_EX  DIAG B: full HFX energy at every APPLY step (expensive)
   ! -----------------------------------------------------------------------
   LOGICAL, PARAMETER, PRIVATE :: DBG_ROUTING = .FALSE.
   LOGICAL, PARAMETER, PRIVATE :: DBG_BUILD = .FALSE.
   LOGICAL, PARAMETER, PRIVATE :: DBG_ENERGY = .FALSE.
   LOGICAL, PARAMETER, PRIVATE :: DBG_STALE = .FALSE.
   LOGICAL, PARAMETER, PRIVATE :: DBG_EXACT_EX = .FALSE.

   LOGICAL, SAVE :: ace_dynamic_mode = .FALSE.
   ! Set to .TRUE. by hfx_ace_set_dynamic_mode before geo_opt/MD starts.
   ! Stays .FALSE. for ENERGY/ENERGY_FORCE single-point runs.

   PUBLIC :: hfx_ace_ks_matrix, hfx_ace_release, hfx_ace_set_dynamic_mode

CONTAINS

! **************************************************************************************************
!> \brief Main ACE entry point, replacing hfx_ks_matrix in qs_ks_methods.
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \param ace_rebuild_frequency ...
!> \param ext_xc_section ...
! **************************************************************************************************
   SUBROUTINE hfx_ace_ks_matrix(qs_env, ks_matrix, rho, energy, &
                                calculate_forces, just_energy, &
                                v_rspace_new, v_tau_rspace, &
                                ace_rebuild_frequency, ext_xc_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(IN)                                :: calculate_forces, just_energy
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_tau_rspace
      INTEGER, INTENT(IN)                                :: ace_rebuild_frequency
      TYPE(section_vals_type), OPTIONAL, POINTER         :: ext_xc_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ace_ks_matrix'

      INTEGER                                            :: handle, iw, n_rep_hf, nspins, &
                                                            rebuild_freq
      LOGICAL                                            :: ace_built_now, rebuild_ace
      REAL(dp)                                           :: ex_ace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      CALL cite_reference(Lin2016ACE)
      NULLIFY (logger, dft_control, x_data, scf_control)

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, x_data=x_data, dft_control=dft_control)
      n_rep_hf = SIZE(x_data, 1)
      nspins = dft_control%nspins

      IF (n_rep_hf /= 1) CPABORT("ACE: only one &HF section is supported.")
      IF (dft_control%nimages /= 1) &
         CPABORT("ACE: k-points / multiple images are not implemented.")

      ! ACE requires explicit MO coefficients (C_occ) which are only available
      ! with diagonalization-based SCF.  OT never constructs mo_coeff during
      ! the SCF, so the ACE projector build loop would silently get garbage.
      CALL get_qs_env(qs_env, scf_control=scf_control)
      IF (scf_control%use_ot) &
         CPABORT("ACE: OT doesn't work, use diagonalization-based SCF.")

      rebuild_freq = MAX(1, ace_rebuild_frequency)

      ! ------------------------------------------------------------------
      ! Bypass A: energy-only call
      ! ------------------------------------------------------------------
      IF (just_energy) THEN
         IF (DBG_ROUTING .AND. iw > 0) &
            WRITE (iw, '(T2,A)') 'ACE | just_energy=T: full HFX (no matrix update)'
         CALL hfx_call(qs_env, ks_matrix, rho, energy, &
                       calculate_forces, just_energy, &
                       v_rspace_new, v_tau_rspace, ext_xc_section)
         CALL timestop(handle)
         RETURN
      END IF

      ! ------------------------------------------------------------------
      ! Bypass B: ionic forces requested
      ! ------------------------------------------------------------------
      IF (calculate_forces) THEN
         IF (DBG_ROUTING .AND. iw > 0) &
            WRITE (iw, '(T2,A)') 'ACE | calculate_forces=T: full HFX for exact forces'
         CALL hfx_call(qs_env, ks_matrix, rho, energy, &
                       calculate_forces, just_energy, &
                       v_rspace_new, v_tau_rspace, ext_xc_section)
         ace_is_built = .FALSE.
         ace_step_counter = 0
         ace_geo_step = ace_geo_step + 1   ! NEW: step 0 done, ACE active from now
         CALL timestop(handle)
         RETURN
      END IF

      ! ------------------------------------------------------------------
      ! Bypass C: first geometry step.
      !
      ! The ATOMIC initial guess gives C_occ far from self-consistency,
      ! which would produce an inaccurate W.  Running full HFX for the
      ! entire first geometry step ensures that wavefunction extrapolation
      ! delivers a near-converged C_occ to geometry step 1, making the
      ! ACE BUILD there accurate from the first application.
      ! ------------------------------------------------------------------
      IF (ace_geo_step == 0 .AND. ace_dynamic_mode) THEN
         IF (iw > 0) WRITE (iw, '(T2,A)') &
            'ACE | geo_step=0 (MD/GEO_OPT): full HFX for reference wavefunction'
         CALL hfx_call(qs_env, ks_matrix, rho, energy, &
                       .FALSE., just_energy, &
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
         WRITE (iw, '(/,T2,A)') REPEAT('-', 56)
         WRITE (iw, '(T2,A)') 'ACE | hfx_ace_ks_matrix'
         WRITE (iw, '(T4,A,L1)') 'ace_is_built  = ', ace_is_built
         WRITE (iw, '(T4,A,L1)') 'rebuild_ace   = ', rebuild_ace
         WRITE (iw, '(T4,A,I6)') 'step_counter  = ', ace_step_counter
         WRITE (iw, '(T4,A,I6)') 'rebuild_freq  = ', rebuild_freq
         WRITE (iw, '(T4,A,I4)') 'nspins        = ', nspins
         WRITE (iw, '(T4,A)') MERGE('-> BUILD', '-> APPLY', rebuild_ace)
         WRITE (iw, '(T2,A)') REPEAT('-', 56)
      END IF

      IF (rebuild_ace) THEN

         ace_built_now = .FALSE.
         CALL hfx_ace_build_projector(qs_env, ks_matrix, rho, energy, &
                                      just_energy, &
                                      v_rspace_new, v_tau_rspace, &
                                      nspins, iw, ace_built_now, &
                                      ext_xc_section)
         IF (ace_built_now) THEN
            ace_is_built = .TRUE.
            ace_step_counter = 1
            IF (DBG_ROUTING .AND. iw > 0) &
               WRITE (iw, '(T4,A)') 'ACE | W built. Projector live from next step.'
         ELSE
            ace_is_built = .FALSE.
            ace_step_counter = 0
            IF (DBG_ROUTING .AND. iw > 0) &
               WRITE (iw, '(T4,A)') 'ACE | Build deferred (C_occ=0). Full HFX in ks_matrix.'
         END IF

      ELSE

         CALL hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy, nspins, iw)
         ace_step_counter = ace_step_counter + 1

         ! ----------------------------------------------------------------
         ! DIAGNOSTIC B: compare E_x^ACE[P^k] with E_x^exact[P^k].
         !
         ! Calls full HFX (just_energy=.TRUE.) to get the exact exchange
         ! energy at the current ACE-converging density P^k.  Compares
         ! with E_x^ACE[P^k] already stored in energy%ex.
         !
         ! Growing |delta| over the SCF confirms the root cause: W was
         ! built from C_occ^(step 1), which is far from the converged
         ! C_occ, so K_ACE = -WW^T no longer represents K_x accurately.
         !
         ! After the comparison, hfx_ace_apply_projector is called a
         ! second time to restore ks_matrix and energy%ex to ACE values
         ! so the SCF continues correctly.  Cost: +2 full HFX per step.
         ! ----------------------------------------------------------------
         IF (DBG_EXACT_EX) THEN
            ex_ace = energy%ex

            CALL hfx_call(qs_env, ks_matrix, rho, energy, &
                          .FALSE., .TRUE., &
                          v_rspace_new, v_tau_rspace, ext_xc_section)
            IF (iw > 0) THEN
               WRITE (iw, '(/,T2,A)') REPEAT('-', 56)
               WRITE (iw, '(T2,A,I6)') 'ACE DIAG B | ace_step_counter = ', ace_step_counter
               WRITE (iw, '(T4,A,F20.10)') 'E_x(exact, P^k)  = ', energy%ex
               WRITE (iw, '(T4,A,F20.10)') 'E_x(ACE,   P^k)  = ', ex_ace
               WRITE (iw, '(T4,A,ES12.4)') '|delta|          = ', ABS(ex_ace - energy%ex)
               WRITE (iw, '(T4,A)') &
                  '|delta|->0 on BUILD step; growth confirms stale projector'
               WRITE (iw, '(T2,A)') REPEAT('-', 56)
            END IF

            ! Restore ACE ks_matrix and energy%ex
            CALL hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy, nspins, iw)
         END IF

      END IF

      IF (DBG_ROUTING .AND. iw > 0) THEN
         WRITE (iw, '(T4,A,F20.10)') 'energy%ex on exit = ', energy%ex
         WRITE (iw, '(T4,A,I6)') 'step_counter now  = ', ace_step_counter
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_ace_ks_matrix

! **************************************************************************************************
!> \brief Build the ACE projector W.
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param energy ...
!> \param just_energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \param nspins ...
!> \param iw ...
!> \param build_succeeded ...
!> \param ext_xc_section ...
! **************************************************************************************************
   SUBROUTINE hfx_ace_build_projector(qs_env, ks_matrix, rho, energy, &
                                      just_energy, &
                                      v_rspace_new, v_tau_rspace, &
                                      nspins, iw, build_succeeded, &
                                      ext_xc_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(IN)                                :: just_energy
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_tau_rspace
      INTEGER, INTENT(IN)                                :: nspins, iw
      LOGICAL, INTENT(OUT)                               :: build_succeeded
      TYPE(section_vals_type), OPTIONAL, POINTER         :: ext_xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ace_build_projector'

      INTEGER                                            :: handle, info_chol, ispin, nao, nmo, nocc
      LOGICAL                                            :: do_admm
      REAL(dp)                                           :: ehfx_full, frob
      REAL(dp), DIMENSION(:), POINTER                    :: occ_nums
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type)                                   :: A_ref_fm, C_occ_fm, K_fm, M_fm, xi_fm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h
      TYPE(dbcsr_type)                                   :: K_ao_dbcsr
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_for_ace
      TYPE(mp_para_env_type), POINTER                    :: para_env

! A_ref_fm: temporary nocc x nocc scratch for DIAG A reference norm

      CALL timeset(routineN, handle)
      NULLIFY (blacs_env, para_env, mos, mo_coeff, matrix_h, occ_nums, fmstruct)
      NULLIFY (dft_control, admm_env, mos_for_ace, ks_aux_fit)

      build_succeeded = .FALSE.

      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, &
                      mos=mos, matrix_h_kp=matrix_h, &
                      dft_control=dft_control)

      do_admm = dft_control%do_admm

      IF (do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, &
                           matrix_ks_aux_fit=ks_aux_fit, &
                           mos_aux_fit=mos_for_ace)
      ELSE
         mos_for_ace => mos
      END IF

      ! Step 1: full HFX
      CALL hfx_call(qs_env, ks_matrix, rho, energy, &
                    .FALSE., just_energy, &
                    v_rspace_new, v_tau_rspace, ext_xc_section)
      ehfx_full = energy%ex

      IF (DBG_BUILD .AND. iw > 0) &
         WRITE (iw, '(/,T2,A,F20.10)') 'ACE BUILD | E_x(full HFX) = ', ehfx_full

      ! Allocate / resize module storage
      IF (ALLOCATED(ace_W)) THEN
         IF (SIZE(ace_W, 2) /= nspins) CALL hfx_ace_release()
      END IF
      IF (.NOT. ALLOCATED(ace_W)) ALLOCATE (ace_W(1, nspins))

      ! Reset reference norm; accumulated per spin in the loop below
      ace_W_ref_norm = 0.0_dp

      ! ----------------------------------------------------------------
      ! Per-spin build loop
      ! ----------------------------------------------------------------
      DO ispin = 1, nspins

         IF (mos_for_ace(ispin)%use_mo_coeff_b) &
            CALL copy_dbcsr_to_fm(mos_for_ace(ispin)%mo_coeff_b, &
                                  mos_for_ace(ispin)%mo_coeff)

         CALL get_mo_set(mos_for_ace(ispin), mo_coeff=mo_coeff, &
                         nao=nao, nmo=nmo, homo=nocc, &
                         occupation_numbers=occ_nums)

         IF (nocc <= 0) CPABORT("ACE: homo <= 0.")
         IF (nocc > nmo) CPABORT("ACE: homo > nmo.")
         CPASSERT(ASSOCIATED(mo_coeff))

         CALL cp_fm_trace(mo_coeff, mo_coeff, frob)

         IF (DBG_BUILD .AND. iw > 0) THEN
            WRITE (iw, '(/,T2,A,I4)') 'ACE BUILD | ispin = ', ispin
            WRITE (iw, '(T4,A,I8)') 'nao  = ', nao
            WRITE (iw, '(T4,A,I8)') 'nmo  = ', nmo
            WRITE (iw, '(T4,A,I8)') 'nocc = ', nocc
            WRITE (iw, '(T4,A,L1)') 'use_mo_coeff_b = ', mos_for_ace(ispin)%use_mo_coeff_b
            WRITE (iw, '(T4,A,ES12.4)') '||mo_coeff||_F = ', SQRT(MAX(frob, 0.0_dp))
         END IF

         IF (frob < 1.0e-20_dp) THEN
            IF (DBG_BUILD .AND. iw > 0) &
               WRITE (iw, '(T4,A)') 'mo_coeff=0: build deferred to next step.'
            CALL timestop(handle)
            RETURN
         END IF

         ! Step 2: K_AO
         IF (do_admm) THEN
            CALL dbcsr_create(K_ao_dbcsr, template=ks_aux_fit(ispin)%matrix, &
                              name="K_ACE_aux")
            CALL dbcsr_copy(K_ao_dbcsr, ks_aux_fit(ispin)%matrix)
         ELSE
            CALL dbcsr_create(K_ao_dbcsr, template=ks_matrix(ispin, 1)%matrix, &
                              name="K_AO")
            CALL dbcsr_copy(K_ao_dbcsr, ks_matrix(ispin, 1)%matrix)
            CALL dbcsr_add(K_ao_dbcsr, matrix_h(1, 1)%matrix, 1.0_dp, -1.0_dp)
         END IF

         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nao, ncol_global=nao)
         CALL cp_fm_create(K_fm, fmstruct, name="K_dense")
         CALL cp_fm_struct_release(fmstruct)
         CALL copy_dbcsr_to_fm(K_ao_dbcsr, K_fm)
         CALL dbcsr_release(K_ao_dbcsr)

         IF (DBG_BUILD .AND. iw > 0) THEN
            CALL cp_fm_trace(K_fm, K_fm, frob)
            WRITE (iw, '(T4,A,ES12.4)') '||K_AO||_F = ', SQRT(MAX(frob, 0.0_dp))
         END IF

         ! Step 3: C_occ
         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nao, ncol_global=nocc)
         CALL cp_fm_create(C_occ_fm, fmstruct, name="C_occ")
         CALL cp_fm_create(xi_fm, fmstruct, name="xi")
         CALL cp_fm_struct_release(fmstruct)

         CALL cp_fm_to_fm(mo_coeff, C_occ_fm)

         ! Step 4: xi = K_AO * C_occ
         CALL cp_fm_gemm('N', 'N', nao, nocc, nao, &
                         1.0_dp, K_fm, C_occ_fm, 0.0_dp, xi_fm)
         CALL cp_fm_release(K_fm)

         IF (DBG_BUILD .AND. iw > 0) THEN
            CALL cp_fm_trace(xi_fm, xi_fm, frob)
            WRITE (iw, '(T4,A,ES12.4)') '||xi||_F = ', SQRT(MAX(frob, 0.0_dp))
         END IF

         ! Step 5: M = C_occ^T * xi
         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nocc, ncol_global=nocc)
         CALL cp_fm_create(M_fm, fmstruct, name="M")
         CALL cp_fm_struct_release(fmstruct)

         CALL cp_fm_gemm('T', 'N', nocc, nocc, nao, &
                         1.0_dp, C_occ_fm, xi_fm, 0.0_dp, M_fm)
         CALL cp_fm_release(C_occ_fm)

         IF (DBG_BUILD .AND. iw > 0) THEN
            CALL cp_fm_trace(M_fm, M_fm, frob)
            WRITE (iw, '(T4,A,ES12.4)') '||M||_F = ', SQRT(MAX(frob, 0.0_dp))
         END IF

         ! Step 6: Cholesky of -M = U^T U
         CALL cp_fm_scale(-1.0_dp, M_fm)
         CALL cp_fm_cholesky_decompose(M_fm, n=nocc, info_out=info_chol)

         IF (info_chol /= 0) THEN
            IF (iw > 0) THEN
               WRITE (iw, '(T4,A,I6)') 'ACE | Cholesky failed, info = ', info_chol
               WRITE (iw, '(T4,A,F20.10)') 'ACE | E_x(full) = ', ehfx_full
               WRITE (iw, '(T4,A,I8,A,I8)') 'ACE | nao=', nao, '  nocc=', nocc
            END IF
            CPABORT("ACE: Cholesky of -M failed (not positive definite).")
         END IF

         IF (DBG_BUILD .AND. iw > 0) &
            WRITE (iw, '(T4,A)') 'Cholesky OK (info=0).'

         ! Step 7: W = xi * U^{-1}
         IF (ASSOCIATED(ace_W(1, ispin)%matrix_struct)) &
            CALL cp_fm_release(ace_W(1, ispin))

         CALL cp_fm_create(ace_W(1, ispin), xi_fm%matrix_struct, name="W_ACE")
         CALL cp_fm_to_fm(xi_fm, ace_W(1, ispin))

         CALL cp_fm_triangular_multiply(M_fm, ace_W(1, ispin), &
                                        side='R', uplo_tr='U', &
                                        transpose_tr=.FALSE., &
                                        invert_tr=.TRUE., &
                                        n_rows=nao, n_cols=nocc, &
                                        alpha=1.0_dp)

         CALL cp_fm_release(xi_fm)
         CALL cp_fm_release(M_fm)

         IF (DBG_BUILD .AND. iw > 0) THEN
            CALL cp_fm_trace(ace_W(1, ispin), ace_W(1, ispin), frob)
            WRITE (iw, '(T4,A,I4,A,2I8,A,ES12.4)') &
               'W spin=', ispin, ' shape=', nao, nocc, &
               '  ||W||_F=', SQRT(MAX(frob, 0.0_dp))
         END IF

         ! ----------------------------------------------------------------
         ! DIAG A reference: compute ||W^T C_occ^BUILD||_F for this spin.
         !
         ! C_occ_fm was released after step 5, but mo_coeff is still valid
         ! (it is a pointer into mos_for_ace, not allocated here).
         ! We re-create C_occ_fm from mo_coeff.
         !
         ! Theory: W^T C_occ^BUILD = U^{-T}(-U^T U) = -U  →  norm = ||U||_F
         ! (computed directly rather than storing U).
         ! ----------------------------------------------------------------
         IF (DBG_STALE) THEN
            NULLIFY (fmstruct)
            CALL cp_fm_struct_create(fmstruct, context=blacs_env, &
                                     para_env=para_env, &
                                     nrow_global=nao, ncol_global=nocc)
            CALL cp_fm_create(C_occ_fm, fmstruct, name="C_occ_ref_diag")
            CALL cp_fm_struct_release(fmstruct)
            CALL cp_fm_to_fm(mo_coeff, C_occ_fm)

            NULLIFY (fmstruct)
            CALL cp_fm_struct_create(fmstruct, context=blacs_env, &
                                     para_env=para_env, &
                                     nrow_global=nocc, ncol_global=nocc)
            CALL cp_fm_create(A_ref_fm, fmstruct, name="WtC_ref")
            CALL cp_fm_struct_release(fmstruct)

            CALL cp_fm_gemm('T', 'N', nocc, nocc, nao, &
                            1.0_dp, ace_W(1, ispin), C_occ_fm, 0.0_dp, A_ref_fm)
            CALL cp_fm_trace(A_ref_fm, A_ref_fm, frob)
            ace_W_ref_norm = ace_W_ref_norm + SQRT(MAX(frob, 0.0_dp))

            CALL cp_fm_release(C_occ_fm)
            CALL cp_fm_release(A_ref_fm)
         END IF

      END DO  ! ispin

      IF (DBG_STALE .AND. iw > 0) THEN
         WRITE (iw, '(/,T2,A)') REPEAT('-', 56)
         WRITE (iw, '(T2,A)') 'ACE DIAG A | Reference norm stored at BUILD'
         WRITE (iw, '(T4,A,ES12.4)') &
            '||W^T C_occ^BUILD||_F (sum over spins) = ', ace_W_ref_norm
         WRITE (iw, '(T4,A)') &
            'Staleness ratio = 1.0 at BUILD step; decreasing means W is becoming stale'
         WRITE (iw, '(T2,A)') REPEAT('-', 56)
      END IF

      build_succeeded = .TRUE.

      ! Step 8: apply immediately (DIAG C ratio printed via ehfx_full_ref)
      CALL hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy, &
                                   nspins, iw, ehfx_full_ref=ehfx_full)

      CALL timestop(handle)

   END SUBROUTINE hfx_ace_build_projector

! **************************************************************************************************
!> \brief Apply the stored ACE projector.
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param energy ...
!> \param nspins ...
!> \param iw ...
!> \param ehfx_full_ref ...
! **************************************************************************************************
   SUBROUTINE hfx_ace_apply_projector(qs_env, ks_matrix, rho, energy, &
                                      nspins, iw, ehfx_full_ref)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_energy_type), POINTER                      :: energy
      INTEGER, INTENT(IN)                                :: nspins, iw
      REAL(dp), INTENT(IN), OPTIONAL                     :: ehfx_full_ref

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ace_apply_projector'

      INTEGER                                            :: handle, ispin, nao, nao_d, nmo_d, nocc, &
                                                            nocc_d
      LOGICAL                                            :: do_admm
      REAL(dp)                                           :: ehfx_ace, frob_A, stale_norm, trace_val
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstruct_diag
      TYPE(cp_fm_type)                                   :: A_diag_fm, C_occ_diag, P_fm, PW_fm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_diag
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_aux_fit, ks_aux_fit_hfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos_aux_diag, mos_diag
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit

! ------------------------------------------------------------------
! DIAG A local variables
!   stale_norm    ||W^T C_occ^current||_F summed over spins
!   frob_A        scratch for cp_fm_trace
!   C_occ_diag    current C_occ redistributed to W layout
!   A_diag_fm     nocc x nocc overlap  W^T C_occ^current
!   mos_diag      primary mos (non-ADMM path)
!   mos_aux_diag  auxiliary mos (ADMM path)
!   mo_coeff_diag pointer to the relevant mo_coeff
!   nao_d, nmo_d, nocc_d  dimensions from get_mo_set
! ------------------------------------------------------------------

      CALL timeset(routineN, handle)
      NULLIFY (blacs_env, para_env, matrix_h, rho_ao, fmstruct)
      NULLIFY (dft_control, admm_env, ks_aux_fit, ks_aux_fit_hfx, rho_aux_fit)
      NULLIFY (mos_diag, mos_aux_diag, mo_coeff_diag, fmstruct_diag)
      CPASSERT(ALLOCATED(ace_W))

      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, &
                      matrix_h_kp=matrix_h, dft_control=dft_control)
      do_admm = dft_control%do_admm

      IF (do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, &
                           matrix_ks_aux_fit=ks_aux_fit, &
                           matrix_ks_aux_fit_hfx=ks_aux_fit_hfx, &
                           rho_aux_fit=rho_aux_fit)
         CALL qs_rho_get(rho_aux_fit, rho_ao_kp=rho_ao)
      ELSE
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_get_info(ace_W(1, ispin), nrow_global=nao, ncol_global=nocc)

         IF (do_admm) THEN
            CALL dbcsr_set(ks_matrix(ispin, 1)%matrix, 0.0_dp)
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, &
                           matrix_h(1, 1)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(ks_aux_fit(ispin)%matrix, 0.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t( &
               sparse_matrix=ks_aux_fit(ispin)%matrix, &
               matrix_v=ace_W(1, ispin), &
               ncol=nocc, &
               alpha=-1.0_dp, &
               keep_sparsity=.TRUE.)
            CALL dbcsr_add(ks_aux_fit_hfx(ispin)%matrix, &
                           ks_aux_fit(ispin)%matrix, 0.0_dp, 1.0_dp)
         ELSE
            CALL dbcsr_set(ks_matrix(ispin, 1)%matrix, 0.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t( &
               sparse_matrix=ks_matrix(ispin, 1)%matrix, &
               matrix_v=ace_W(1, ispin), &
               ncol=nocc, &
               alpha=-1.0_dp, &
               keep_sparsity=.TRUE.)
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, &
                           matrix_h(1, 1)%matrix, 1.0_dp, 1.0_dp)
         END IF
      END DO

      ! Exchange energy: E_x = -0.5 * sum_spin Tr[ W^T * P * W ]
      ehfx_ace = 0.0_dp
      DO ispin = 1, nspins
         CALL cp_fm_get_info(ace_W(1, ispin), nrow_global=nao, ncol_global=nocc)

         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nao, ncol_global=nao)
         CALL cp_fm_create(P_fm, fmstruct, name="P_dense")
         CALL cp_fm_struct_release(fmstruct)
         CALL copy_dbcsr_to_fm(rho_ao(ispin, 1)%matrix, P_fm)

         CALL cp_fm_create(PW_fm, ace_W(1, ispin)%matrix_struct, name="PW")
         CALL cp_fm_gemm('N', 'N', nao, nocc, nao, &
                         1.0_dp, P_fm, ace_W(1, ispin), 0.0_dp, PW_fm)
         CALL cp_fm_trace(ace_W(1, ispin), PW_fm, trace_val)
         ehfx_ace = ehfx_ace - 0.5_dp*trace_val

         CALL cp_fm_release(P_fm)
         CALL cp_fm_release(PW_fm)

         IF (DBG_ENERGY .AND. iw > 0) &
            WRITE (iw, '(T4,A,I4,A,F20.10)') &
            'ispin=', ispin, '  E_x(ACE) += ', -0.5_dp*trace_val
      END DO

      energy%ex = ehfx_ace

      ! DIAG C: BUILD-step consistency check (printed when ehfx_full_ref present)
      IF (DBG_ENERGY .AND. iw > 0) THEN
         WRITE (iw, '(T2,A,F20.10)') 'ACE | E_x(ACE)  = ', ehfx_ace
         IF (PRESENT(ehfx_full_ref)) THEN
            WRITE (iw, '(T2,A,F20.10)') 'ACE | E_x(full) = ', ehfx_full_ref
            WRITE (iw, '(T2,A,ES12.4)') 'ACE | |delta|   = ', ABS(ehfx_ace - ehfx_full_ref)
            WRITE (iw, '(T2,A)') '(|delta| should be ~0 on BUILD steps; small is good)'
         END IF
      END IF

      ! ----------------------------------------------------------------
      ! DIAG A: projector staleness check.
      !
      ! Computes ||W^T C_occ^current||_F (summed over spins) and divides
      ! by ace_W_ref_norm = ||W^T C_occ^BUILD||_F stored at BUILD time.
      !
      ! staleness_ratio:
      !   1.0  → C_occ hasn't changed since BUILD; projector is fresh
      !   < 1  → C_occ has rotated; how much depends on the SCF dynamics
      !   → 0  → C_occ is orthogonal to the BUILD-time span; W is useless
      !
      ! For non-ADMM: C_occ comes from primary mos (nao_orb x nocc).
      ! For ADMM:     C_occ comes from mos_aux_fit (nao_aux x nocc_aux),
      !               consistent with ace_W dimensions.
      ! ----------------------------------------------------------------
      IF (DBG_STALE .AND. ace_W_ref_norm > 0.0_dp) THEN
         stale_norm = 0.0_dp
         CALL get_qs_env(qs_env, mos=mos_diag)

         DO ispin = 1, nspins
            CALL cp_fm_get_info(ace_W(1, ispin), nrow_global=nao_d, ncol_global=nocc_d)

            IF (do_admm) THEN
               CALL get_admm_env(admm_env, mos_aux_fit=mos_aux_diag)
               IF (mos_aux_diag(ispin)%use_mo_coeff_b) &
                  CALL copy_dbcsr_to_fm(mos_aux_diag(ispin)%mo_coeff_b, &
                                        mos_aux_diag(ispin)%mo_coeff)
               CALL get_mo_set(mos_aux_diag(ispin), mo_coeff=mo_coeff_diag, &
                               nao=nao_d, nmo=nmo_d, homo=nocc_d)
            ELSE
               IF (mos_diag(ispin)%use_mo_coeff_b) &
                  CALL copy_dbcsr_to_fm(mos_diag(ispin)%mo_coeff_b, &
                                        mos_diag(ispin)%mo_coeff)
               CALL get_mo_set(mos_diag(ispin), mo_coeff=mo_coeff_diag, &
                               nao=nao_d, nmo=nmo_d, homo=nocc_d)
            END IF

            NULLIFY (fmstruct_diag)
            CALL cp_fm_struct_create(fmstruct_diag, context=blacs_env, &
                                     para_env=para_env, &
                                     nrow_global=nao_d, ncol_global=nocc_d)
            CALL cp_fm_create(C_occ_diag, fmstruct_diag, name="C_stale")
            CALL cp_fm_struct_release(fmstruct_diag)
            CALL cp_fm_to_fm(mo_coeff_diag, C_occ_diag)

            NULLIFY (fmstruct_diag)
            CALL cp_fm_struct_create(fmstruct_diag, context=blacs_env, &
                                     para_env=para_env, &
                                     nrow_global=nocc_d, ncol_global=nocc_d)
            CALL cp_fm_create(A_diag_fm, fmstruct_diag, name="WtC_stale")
            CALL cp_fm_struct_release(fmstruct_diag)

            CALL cp_fm_gemm('T', 'N', nocc_d, nocc_d, nao_d, &
                            1.0_dp, ace_W(1, ispin), C_occ_diag, 0.0_dp, A_diag_fm)
            CALL cp_fm_trace(A_diag_fm, A_diag_fm, frob_A)
            stale_norm = stale_norm + SQRT(MAX(frob_A, 0.0_dp))

            CALL cp_fm_release(C_occ_diag)
            CALL cp_fm_release(A_diag_fm)
         END DO
         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') REPEAT('-', 56)
            WRITE (iw, '(T2,A,I6)') 'ACE DIAG A | ace_step_counter      = ', ace_step_counter
            WRITE (iw, '(T4,A,ES12.4)') '||W^T C_occ^current||_F          = ', stale_norm
            WRITE (iw, '(T4,A,ES12.4)') '||W^T C_occ^BUILD||_F  (ref)     = ', ace_W_ref_norm
            WRITE (iw, '(T4,A,F10.6)') 'staleness ratio (1=fresh, 0=stale) = ', &
               stale_norm/MAX(ace_W_ref_norm, 1.0e-30_dp)
            WRITE (iw, '(T2,A)') REPEAT('-', 56)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_ace_apply_projector

! **************************************************************************************************
!> \brief Release all ACE storage and reset state flags.
!> \param iw_opt ...
! **************************************************************************************************
   SUBROUTINE hfx_ace_release(iw_opt)

      INTEGER, INTENT(IN), OPTIONAL                      :: iw_opt

      INTEGER                                            :: i, iw, j

      iw = -1
      IF (PRESENT(iw_opt)) iw = iw_opt

      IF (ALLOCATED(ace_W)) THEN
         DO j = 1, SIZE(ace_W, 2)
            DO i = 1, SIZE(ace_W, 1)
               IF (ASSOCIATED(ace_W(i, j)%matrix_struct)) CALL cp_fm_release(ace_W(i, j))
            END DO
         END DO
         DEALLOCATE (ace_W)
      END IF

      ace_is_built = .FALSE.
      ace_step_counter = 0
      ace_W_ref_norm = 0.0_dp
      ace_geo_step = 0
      ace_dynamic_mode = .FALSE.   ! ADDED: reset dynamic mode on release, so it must be explicitly re-enabled for GEO_OPT/MD runs

      IF (iw > 0) WRITE (iw, '(T2,A)') 'ACE | storage released, counters reset'

   END SUBROUTINE hfx_ace_release

   ! **************************************************************************************************
   !> \brief Mark this run as dynamic (GEO_OPT/MD) so Bypass C fires for geo step 0.
   !>        Call this once from the geo_opt or MD driver before the first SCF.
   !> \param is_dynamic .TRUE. for GEO_OPT/MD, .FALSE. to reset.
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param is_dynamic ...
! **************************************************************************************************
   SUBROUTINE hfx_ace_set_dynamic_mode(is_dynamic)
      LOGICAL, INTENT(IN)                                :: is_dynamic

      ace_dynamic_mode = is_dynamic
   END SUBROUTINE hfx_ace_set_dynamic_mode

   ! **************************************************************************************************
   !> \brief Private helper: call hfx_ks_matrix with or without ext_xc_section.
   !> \param qs_env ...
   !> \param ks_matrix ...
   !> \param rho ...
   !> \param energy ...
   !> \param calculate_forces ...
   !> \param just_energy ...
   !> \param v_rspace_new ...
   !> \param v_tau_rspace ...
   !> \param ext_xc_section ...
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \param ext_xc_section ...
! **************************************************************************************************
   SUBROUTINE hfx_call(qs_env, ks_matrix, rho, energy, &
                       calculate_forces, just_energy, &
                       v_rspace_new, v_tau_rspace, ext_xc_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(IN)                                :: calculate_forces, just_energy
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_tau_rspace
      TYPE(section_vals_type), OPTIONAL, POINTER         :: ext_xc_section

      IF (PRESENT(ext_xc_section)) THEN
         CALL hfx_ks_matrix(qs_env, ks_matrix, rho, energy, &
                            calculate_forces, just_energy, &
                            v_rspace_new, v_tau_rspace, &
                            ext_xc_section=ext_xc_section)
      ELSE
         CALL hfx_ks_matrix(qs_env, ks_matrix, rho, energy, &
                            calculate_forces, just_energy, &
                            v_rspace_new, v_tau_rspace)
      END IF

   END SUBROUTINE hfx_call

END MODULE hfx_ace_methods

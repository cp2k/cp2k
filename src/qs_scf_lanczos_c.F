!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex block-Lanczos refinement kernel for k-point SCF
!>        diagonalization. Processes unconverged columns in blocks of
!>        NBLOCK with two-pass DGKS reorthogonalization against the
!>        working basis and deflation prefix, BLAS-3 CholeskyQR2 batch
!>        orthonormalization (with CGS2 fallback on rank-deficient
!>        batches), and an explicit Hermitian Rayleigh-Ritz projection
!>        T = V^dag(A V). Scratch matrices and local LAPACK buffers live
!>        in a caller-owned workspace reused across sweeps.
! **************************************************************************************************
MODULE qs_scf_lanczos_c
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_gemm,&
                                              cp_cfm_triangular_multiply,&
                                              cp_cfm_vectorsnorm
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_mone,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_lanczos_c'
   INTEGER, PARAMETER, PRIVATE :: shared_budget_factor = 4

   PUBLIC :: lanczos_refinement_2v_c
   PUBLIC :: lanczos_channel_begin, lanczos_sweep_suffix, lanczos_channel_end
   PUBLIC :: lanczos_workspace_create, lanczos_workspace_release

! **************************************************************************************************
!> \brief Caller-owned scratch of the refinement kernel: every matrix and
!>        array that the sweeps of one channel reuse. The driver creates
!>        the workspace once per diagonalization call and hands it to
!>        every kernel call. The widths follow the worst case of the run.
!>        The basis budget holds the larger of the NKRYLOV block budget
!>        and four times the widest single block, and every per-call
!>        clamp stays below it. The small Rayleigh-Ritz pair keeps its
!>        per-block right-sized allocation, because the solver stack
!>        diagonalizes whole matrices and knows no leading-block argument
! **************************************************************************************************
   TYPE, PUBLIC :: lanczos_workspace_type
      ! entry Rayleigh-Ritz, its rotation and the shared scratch
      TYPE(cp_cfm_type) :: hmo, mo_gram, evec, rot_scr
      ! working basis and its A-images
      TYPE(cp_cfm_type) :: vmat, hvmat
      ! candidate block, within-batch scratch column, batch Gram, deflation and
      ! projection-coefficient scratch
      TYPE(cp_cfm_type) :: cand_blk, one_col, batch_gram, proj_coeff, defl_coeff
      REAL(dp), ALLOCATABLE :: theta(:), nrm(:), rank_scale(:), pair_res(:), entry_res(:)
      LOGICAL, ALLOCATABLE :: ldone(:)
      COMPLEX(dp), ALLOCATABLE :: eval_scale(:)
   END TYPE lanczos_workspace_type

CONTAINS

! **************************************************************************************************
!> \brief Refine the non-converged columns of c0 on the standardized operator ks.
!>        The entry Rayleigh-Ritz pass updates the convergence statistics.
!>        Every non-converged column block then runs a single-pass
!>        Lanczos refinement. A full Cholesky sweep closes the call.
!> \param ks standardized operator A = U^-dag H U^-1, nao x nao
!> \param c0 current MOS in the standardized basis, refined in place.
!>        The matrix dimensions of the call follow c0
!> \param eval eigenvalues of c0, updated in place by the entry Rayleigh-Ritz
!> \param eps_col per-column convergence thresholds on the residual norms:
!>        the driver assigns the relaxed empty threshold to the unoccupied
!>        columns. The occupation decides, not the column index. This
!>        mirrors the Davidson k-point driver
!> \param nkrylov number of block expansion steps in the block regime.
!>        The kernel uses it as a floor on the budget only when NBLOCK
!>        covers every pending column
!> \param nblock number of columns refined jointly, one or more
!> \param ws caller-owned scratch of the call, created by the driver at
!>        the worst-case widths of the run
!> \param max_res_norm largest residual column norm on exit: the entry
!>        Rayleigh-Ritz values, or the wanted-pair residuals of the
!>        refinement blocks when a refinement ran
!> \param min_res_norm smallest residual column norm on exit (same basis)
!> \param nmo_conv number of leading converged columns on exit
!> \param nmo_nc columns outside the leading converged prefix: every
!>        still-unconverged column, plus converged columns beyond a gap.
!>        The count stays conservative for the not-converged reporting
!> \param check_moconv_only only run the entry Rayleigh-Ritz and the
!>        convergence count, skip the refinement
! **************************************************************************************************
   SUBROUTINE lanczos_refinement_2v_c(ks, c0, eval, eps_col, nkrylov, nblock, ws, &
                                      max_res_norm, min_res_norm, nmo_conv, nmo_nc, check_moconv_only)

      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      INTEGER, INTENT(OUT)                               :: nmo_conv, nmo_nc
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_moconv_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_refinement_2v_c'

      INTEGER                                            :: handle
      LOGICAL                                            :: my_check_moconv_only, next_suffix_rr, &
                                                            refined

      CALL timeset(routineN, handle)

      my_check_moconv_only = .FALSE.
      IF (PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

      CALL lanczos_channel_begin(ks, c0, eval, eps_col, ws, max_res_norm, min_res_norm, &
                                 nmo_conv, nmo_nc, check_moconv_only=my_check_moconv_only)

      IF (.NOT. my_check_moconv_only) THEN
         refined = .FALSE.
         IF (nmo_nc > 0) THEN
            CALL lanczos_sweep_suffix(ks, c0, eval, eps_col, nkrylov, nblock, .TRUE., .FALSE., &
                                      ws, max_res_norm, min_res_norm, nmo_conv, nmo_nc, &
                                      refined, next_suffix_rr)
         END IF
         IF (refined) THEN
            CALL lanczos_channel_end(c0, ws)
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE lanczos_refinement_2v_c

! **************************************************************************************************
!> \brief Entry preparation of one (k-point, spin) channel before its inner
!>        Lanczos sweeps: orthonormalizes the input MOS (skipped for a
!>        read-only convergence probe), runs the full-width entry
!>        Rayleigh-Ritz pass, and computes the initial leading converged
!>        prefix and residual bounds.
!> \param ks standardized operator A = U^-dag H U^-1, nao x nao
!> \param c0 current MOS in the standardized basis, rotated in place
!> \param eval eigenvalues of c0, updated in place
!> \param eps_col per-column convergence thresholds
!> \param ws caller-owned channel workspace
!> \param max_res_norm largest entry residual norm
!> \param min_res_norm smallest entry residual norm
!> \param nmo_conv number of leading converged columns on entry
!> \param nmo_nc number of columns outside the leading converged prefix
!> \param check_moconv_only skip the entry Cholesky orthonormalization for
!>        a read-only probe
! **************************************************************************************************
   SUBROUTINE lanczos_channel_begin(ks, c0, eval, eps_col, ws, max_res_norm, min_res_norm, &
                                    nmo_conv, nmo_nc, check_moconv_only)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      INTEGER, INTENT(OUT)                               :: nmo_conv, nmo_nc
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_moconv_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_channel_begin'

      INTEGER                                            :: handle, info, iwork, nao, nmo
      LOGICAL                                            :: my_check_moconv_only

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(c0, nrow_global=nao, ncol_global=nmo)

      my_check_moconv_only = .FALSE.
      IF (PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

      IF (.NOT. my_check_moconv_only) THEN
         CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, ws%mo_gram)
         CALL cp_cfm_cholesky_decompose(ws%mo_gram, nmo, info)
         IF (info /= 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "The Gram matrix of the entry MOS lost positive definiteness.")
         END IF
         CALL cp_cfm_triangular_multiply(ws%mo_gram, c0, 'R', invert_tr=.TRUE.)
      END IF

      CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, ks, c0, z_zero, ws%hmo)
      CALL rayleigh_ritz(ws%hmo, c0, nmo, nmo, eval, ws%mo_gram, ws%evec, &
                         ws%eval_scale, ws%rot_scr, 1, ws%entry_res)

      nmo_conv = 0
      nmo_nc = 0
      max_res_norm = MAXVAL(ws%entry_res(1:nmo))
      min_res_norm = MINVAL(ws%entry_res(1:nmo))
      DO iwork = 1, nmo
         IF (ws%entry_res(iwork) <= eps_col(iwork)) THEN
            nmo_conv = nmo_conv + 1
         ELSE
            nmo_nc = nmo - nmo_conv
            EXIT
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE lanczos_channel_begin

! **************************************************************************************************
!> \brief Run one inner Lanczos sweep on the unconverged suffix of c0.
!>        When the previous sweep ran in multi-block mode or had a
!>        truncated block, a local suffix Rayleigh-Ritz pass first
!>        decouples the unconverged tail of width nmo_nc without touching
!>        the locked prefix 1..nmo_conv.
!> \param ks standardized operator, nao x nao
!> \param c0 MOS, refined in place on columns nmo_conv+1..nmo
!> \param eval eigenvalues, updated in place on the unconverged suffix
!> \param eps_col per-column convergence thresholds
!> \param nkrylov number of block expansion steps in the block regime
!> \param nblock number of columns refined jointly
!> \param first_sweep true on the first sweep right after lanczos_channel_begin
!> \param need_suffix_rr whether to run a local suffix Rayleigh-Ritz before
!>        refining
!> \param ws caller-owned channel workspace
!> \param max_res_norm largest residual norm of the active suffix on exit
!> \param min_res_norm smallest residual norm of the active suffix on exit
!> \param nmo_conv leading converged column count, advanced in place
!> \param nmo_nc unconverged suffix width, updated in place
!> \param step_refined whether this sweep modified c0
!> \param next_need_suffix_rr whether the next sweep requires a suffix RR
! **************************************************************************************************
   SUBROUTINE lanczos_sweep_suffix(ks, c0, eval, eps_col, nkrylov, nblock, first_sweep, &
                                   need_suffix_rr, ws, max_res_norm, min_res_norm, &
                                   nmo_conv, nmo_nc, step_refined, next_need_suffix_rr)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock
      LOGICAL, INTENT(IN)                                :: first_sweep, need_suffix_rr
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws
      REAL(dp), INTENT(INOUT)                            :: max_res_norm, min_res_norm
      INTEGER, INTENT(INOUT)                             :: nmo_conv, nmo_nc
      LOGICAL, INTENT(OUT)                               :: step_refined, next_need_suffix_rr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_sweep_suffix'

      INTEGER                                            :: handle, n_conv_cols, nao, nmo_converged, &
                                                            nmo_nc_loc
      LOGICAL                                            :: refined, skip_first_defl
      REAL(dp)                                           :: prev_max, prev_min

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(c0, nrow_global=nao)

      step_refined = .FALSE.
      next_need_suffix_rr = .FALSE.

      IF (need_suffix_rr .AND. nmo_nc > 0) THEN
         CALL suffix_rayleigh_ritz(ks, c0, eval, eps_col, ws, nmo_conv, nmo_nc, &
                                   max_res_norm, min_res_norm)
         step_refined = .TRUE.
         IF (nmo_nc == 0) THEN
            CALL timestop(handle)
            RETURN
         END IF
      END IF

      IF (nmo_nc > 0) THEN
         prev_max = max_res_norm
         prev_min = min_res_norm
         nmo_converged = nmo_conv
         nmo_nc_loc = nmo_nc
         skip_first_defl = first_sweep .OR. need_suffix_rr
         CALL lanczos_refine(ks, c0, eval, nao, eps_col, nkrylov, nblock, &
                             nmo_converged, nmo_nc_loc, skip_first_defl, n_conv_cols, &
                             max_res_norm, min_res_norm, refined, next_need_suffix_rr, ws)
         IF (.NOT. refined) THEN
            max_res_norm = prev_max
            min_res_norm = prev_min
         ELSE
            step_refined = .TRUE.
         END IF
         nmo_conv = nmo_converged + n_conv_cols
         nmo_nc = nmo_nc_loc - n_conv_cols
      END IF

      CALL timestop(handle)
   END SUBROUTINE lanczos_sweep_suffix

! **************************************************************************************************
!> \brief Defensive full-width Cholesky orthonormalization at the end of a
!>        channel's refinement sweeps, executed at most once per channel.
!> \param c0 refined MOS in the standardized basis
!> \param ws caller-owned channel workspace
! **************************************************************************************************
   SUBROUTINE lanczos_channel_end(c0, ws)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_channel_end'

      INTEGER                                            :: handle, info, nao, nmo

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(c0, nrow_global=nao, ncol_global=nmo)
      CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, ws%mo_gram)
      CALL cp_cfm_cholesky_decompose(ws%mo_gram, nmo, info)
      IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "The Gram matrix of the refined set lost positive definiteness.")
      END IF
      CALL cp_cfm_triangular_multiply(ws%mo_gram, c0, 'R', invert_tr=.TRUE.)

      CALL timestop(handle)
   END SUBROUTINE lanczos_channel_end

! **************************************************************************************************
!> \brief Orthonormalize and run a local Rayleigh-Ritz pass on the
!>        unconverged suffix c0(:, nmo_conv+1:nmo) after a multi-block or
!>        truncated sweep, advancing nmo_conv if newly decoupled pairs meet
!>        eps_col.
!> \param ks standardized operator, nao x nao
!> \param c0 MOS, updated in place on the suffix nmo_conv+1..nmo
!> \param eval eigenvalues, updated in place on the suffix
!> \param eps_col per-column convergence thresholds
!> \param ws caller-owned channel workspace
!> \param nmo_conv leading converged count, advanced in place
!> \param nmo_nc unconverged suffix width, updated in place
!> \param max_res_norm largest residual norm of the suffix
!> \param min_res_norm smallest residual norm of the suffix
! **************************************************************************************************
   SUBROUTINE suffix_rayleigh_ritz(ks, c0, eval, eps_col, ws, nmo_conv, nmo_nc, &
                                   max_res_norm, min_res_norm)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws
      INTEGER, INTENT(INOUT)                             :: nmo_conv, nmo_nc
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm

      INTEGER                                            :: col_lo, info, ipass, iwork, m, nao, nmo
      TYPE(cp_cfm_type)                                  :: esolve, hsolve
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp

      NULLIFY (fm_struct_tmp)
      CALL cp_cfm_get_info(c0, nrow_global=nao, ncol_global=nmo)

      col_lo = nmo_conv + 1
      m = nmo_nc

      CALL cp_cfm_to_cfm(c0, ws%rot_scr, m, col_lo, 1)
      IF (nmo_conv > 0) THEN
         DO ipass = 1, 2
            CALL cp_cfm_gemm('C', 'N', nmo_conv, m, nao, z_one, c0, ws%rot_scr, z_zero, ws%mo_gram)
            CALL cp_cfm_gemm('N', 'N', nao, m, nmo_conv, z_mone, c0, ws%mo_gram, z_one, ws%rot_scr)
         END DO
      END IF
      CALL cp_cfm_gemm('C', 'N', m, m, nao, z_one, ws%rot_scr, ws%rot_scr, z_zero, ws%mo_gram)
      CALL cp_cfm_cholesky_decompose(ws%mo_gram, m, info)
      IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "The Gram matrix of the suffix MOS lost positive definiteness.")
      END IF
      CALL cp_cfm_triangular_multiply(ws%mo_gram, ws%rot_scr, side='R', invert_tr=.TRUE., &
                                      n_rows=nao, n_cols=m)
      CALL cp_cfm_to_cfm(ws%rot_scr, c0, m, 1, col_lo)

      CALL cp_cfm_gemm('N', 'N', nao, m, nao, z_one, ks, c0, z_zero, ws%hmo, b_first_col=col_lo)
      IF (m == nmo) THEN
         CALL rayleigh_ritz(ws%hmo, c0, m, m, eval(col_lo:nmo), ws%mo_gram, ws%evec, &
                            ws%eval_scale, ws%rot_scr, col_lo, ws%entry_res, bas_col=col_lo)
      ELSE
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=m, ncol_global=m, &
                                  para_env=c0%matrix_struct%para_env, &
                                  context=c0%matrix_struct%context)
         CALL cp_cfm_create(hsolve, matrix_struct=fm_struct_tmp, name="lan_sh", set_zero=.TRUE.)
         CALL cp_cfm_create(esolve, matrix_struct=fm_struct_tmp, name="lan_se", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL rayleigh_ritz(ws%hmo, c0, m, m, eval(col_lo:nmo), hsolve, esolve, &
                            ws%eval_scale, ws%rot_scr, col_lo, ws%entry_res, bas_col=col_lo)
         CALL cp_cfm_release(hsolve)
         CALL cp_cfm_release(esolve)
      END IF

      max_res_norm = MAXVAL(ws%entry_res(1:m))
      min_res_norm = MINVAL(ws%entry_res(1:m))
      DO iwork = 1, m
         IF (ws%entry_res(iwork) <= eps_col(col_lo + iwork - 1)) THEN
            nmo_conv = nmo_conv + 1
         ELSE
            EXIT
         END IF
      END DO
      nmo_nc = nmo - nmo_conv
   END SUBROUTINE suffix_rayleigh_ritz

! **************************************************************************************************
!> \brief Per-block single-pass Lanczos refinement on the unconverged tail.
!>        Processes non-converged columns in blocks of nblock, expanding
!>        each seed block by A-images with two-pass DGKS reorthogonalization
!>        against the working basis and deflation prefix, and solving the
!>        projected Hermitian Rayleigh-Ritz problem T = V^dag(A V).
!> \param ks standardized operator, nao x nao
!> \param c0 MOS, refined in place block by block
!> \param eval eigenvalues, updated per converged column
!> \param nao number of atomic orbitals
!> \param eps_col per-column convergence thresholds
!> \param nkrylov number of block expansion steps in the block regime
!> \param nblock number of columns refined jointly
!> \param nmo_converged number of leading locked columns
!> \param nmo_nc_loc number of non-converged tail columns
!> \param skip_first_defl skip deflation on the first seed block when the
!>        suffix just underwent a Rayleigh-Ritz rotation
!> \param n_conv_cols length of the leading converged run across blocks
!> \param max_res_norm largest residual norm among the wanted pairs
!> \param min_res_norm smallest residual norm among the wanted pairs
!> \param refined whether any block completed a Rayleigh-Ritz pass
!> \param need_suffix_rr whether a subsequent sweep must run a suffix RR
!>        (multi-block run or any block truncated below its target width)
!> \param ws caller-owned scratch of the call. Only the small
!>        Rayleigh-Ritz pair stays local to a block, at its right size,
!>        because the solver stack diagonalizes whole matrices
! **************************************************************************************************
   SUBROUTINE lanczos_refine(ks, c0, eval, nao, eps_col, nkrylov, nblock, &
                             nmo_converged, nmo_nc_loc, skip_first_defl, n_conv_cols, &
                             max_res_norm, min_res_norm, refined, need_suffix_rr, ws)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      INTEGER, INTENT(IN)                                :: nao
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock, nmo_converged, &
                                                            nmo_nc_loc
      LOGICAL, INTENT(IN)                                :: skip_first_defl
      INTEGER, INTENT(OUT)                               :: n_conv_cols
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      LOGICAL, INTENT(OUT)                               :: refined, need_suffix_rr
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lanczos_refine'

      INTEGER                                            :: col_hi, col_lo, handle, ib, j, nbasis, &
                                                            nblock_e, nbmax, nmob_c, nmobt, nres, &
                                                            nseed, nstep
      LOGICAL                                            :: any_truncated, prefix_open, single_block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: esolve, hsolve
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct_tmp, blacs_env, para_env)
      para_env => c0%matrix_struct%para_env
      blacs_env => c0%matrix_struct%context

      nblock_e = MIN(MAX(nblock, 1), MAX(nmo_nc_loc, 1))
      single_block = nblock_e >= nmo_nc_loc
      IF (single_block) THEN
         ! one block covers every pending column: the Krylov budget holds
         ! about four times the pending width, floored at NKRYLOV and at
         ! the block width. The shared subspace holds every target pair
         ! plus a threefold expansion room. NKRYLOV plays no role here
         nbmax = MIN(nao - 2, MAX(shared_budget_factor*nmo_nc_loc, nkrylov, nblock_e + 2))
      ELSE
         nbmax = nkrylov*nblock_e
      END IF
      IF (nbmax < 2) THEN
         ! no working basis beyond the seed exists (nao below four, or a
         ! degenerate NKRYLOV/NBLOCK combination): leave the columns at
         ! their entry state for the caller to report, instead of aborting
         n_conv_cols = 0
         max_res_norm = 0.0_dp
         min_res_norm = HUGE(1.0_dp)
         refined = .FALSE.
         need_suffix_rr = .FALSE.
         CPWARN("Krylov budget below two columns, Lanczos refinement skipped")
         CALL timestop(handle)
         RETURN
      END IF
      ! the budget bounds the block width: a wider block would append
      ! past the end of vmat at its first seed call (nao - 2 caps nbmax
      ! below the pending width on minimal bases). The pending set then
      ! splits into narrower blocks, and no column is dropped
      IF (nblock_e > nbmax - 2) THEN
         nblock_e = MAX(nbmax - 2, 1)
         single_block = nblock_e >= nmo_nc_loc
      END IF

      n_conv_cols = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      prefix_open = .TRUE.
      refined = .FALSE.
      any_truncated = .FALSE.

      block_loop: DO ib = 1, nmo_nc_loc, nblock_e
         col_lo = nmo_converged + ib
         col_hi = MIN(col_lo + nblock_e - 1, nmo_converged + nmo_nc_loc)
         nmobt = col_hi - col_lo + 1
         ws%ldone(1:nmobt) = .FALSE.
         nres = 0

         ! ---------- seed block: the pending columns, purified and
         ! orthonormalized into the working basis
         CALL cp_cfm_to_cfm(c0, ws%cand_blk, nmobt, col_lo, 1)
         CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)
         ws%rank_scale(1:nmobt) = ws%nrm(1:nmobt)
         nbasis = 0
         CALL purify_and_append(c0, col_lo - 1, ib == 1 .AND. skip_first_defl, ws, &
                                MIN(nmobt, nbmax), nbasis, nseed)
         IF (nseed == 0) THEN
            ! no seed survived the purification: no column of this block
            ! converged, and the leading prefix cannot extend past it
            prefix_open = .FALSE.
            any_truncated = .TRUE.
            CYCLE block_loop
         END IF
         nstep = nseed

         ! ---------- expansion: extend the basis with the A-images of the
         ! last batch until the budget is exhausted
         DO WHILE (nbasis < nbmax .AND. nstep > 0)
            CALL cp_cfm_gemm('N', 'N', nao, nstep, nao, z_one, ks, ws%vmat, z_zero, ws%hvmat, &
                             b_first_col=nbasis - nstep + 1, c_first_col=nbasis - nstep + 1)
            CALL cp_cfm_to_cfm(ws%hvmat, ws%cand_blk, nstep, nbasis - nstep + 1, 1)
            CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)
            ws%rank_scale(1:nstep) = ws%nrm(1:nstep)
            nmob_c = MIN(nstep, nbmax - nbasis)
            CALL purify_and_append(c0, col_lo - 1, .FALSE., ws, nmob_c, nbasis, nseed)
            nstep = nseed
         END DO
         ! the last batch has no successor: its A-images must still be
         ! formed explicitly, or the projection reads the stale tail of
         ! the reused product buffer and silently degrades the Ritz pairs
         IF (nstep > 0) THEN
            CALL cp_cfm_gemm('N', 'N', nao, nstep, nao, z_one, ks, ws%vmat, z_zero, ws%hvmat, &
                             b_first_col=nbasis - nstep + 1, c_first_col=nbasis - nstep + 1)
         END IF

         ! ---------- Rayleigh-Ritz: T = V^dag (A V) as a distributed
         ! problem on a right-sized matrix, through the shared solver
         ! stack
         NULLIFY (fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbasis, ncol_global=nbasis, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_cfm_create(hsolve, matrix_struct=fm_struct_tmp, name="lan_h", set_zero=.TRUE.)
         CALL cp_cfm_create(esolve, matrix_struct=fm_struct_tmp, name="lan_e", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct_tmp)
         nres = MIN(nmobt, nbasis)
         IF (nres < nmobt) any_truncated = .TRUE.
         ! one shared Rayleigh-Ritz pass per block, on the stored
         ! A-images of the working basis: rotates the leading nres Ritz
         ! vectors into c0(:, col_lo:col_lo+nres-1) and forms their exact
         ! residual norms in ws%nrm
         CALL rayleigh_ritz(ws%hvmat, ws%vmat, nbasis, nres, ws%theta, hsolve, esolve, &
                            ws%eval_scale, ws%cand_blk, col_lo, ws%nrm, c_dst=c0)
         ws%pair_res(1:nres) = ws%nrm(1:nres)

         ! ---------- write-back metadata: no pair-to-column identity ever
         ! latches across calls, heevd re-sorts ascending. Columns beyond
         ! nres keep their entry vectors and fall to the end-of-call
         ! fallback
         DO j = 1, nres
            eval(col_lo + j - 1) = ws%theta(j)
            ws%ldone(j) = ws%pair_res(j) < eps_col(col_lo + j - 1)
         END DO
         refined = .TRUE.
         CALL cp_cfm_release(hsolve)
         CALL cp_cfm_release(esolve)

         ! exit statistics: the residuals describe the state the columns
         ! are left in, not the worst value seen during the call. The
         ! converged count is the global leading converged run
         DO j = 1, nres
            max_res_norm = MAX(max_res_norm, ws%pair_res(j))
            min_res_norm = MIN(min_res_norm, ws%pair_res(j))
         END DO
         IF (prefix_open) THEN
            DO j = 1, nmobt
               IF (.NOT. ws%ldone(j)) EXIT
               n_conv_cols = n_conv_cols + 1
            END DO
            prefix_open = ALL(ws%ldone(1:nmobt))
         END IF
      END DO block_loop

      need_suffix_rr = (.NOT. single_block) .OR. any_truncated

      CALL timestop(handle)
   END SUBROUTINE lanczos_refine

! **************************************************************************************************
!> \brief Purify a block of candidate directions and append the survivors to
!>        the working basis: projection passes against the working basis and
!>        the leading ndefl columns of c0, with the second pass decided by
!>        the DGKS norm-ratio test, modified Gram-Schmidt within the block,
!>        then a per-column rank safeguard relative to rank_scale,
!>        normalization and append. The count of appended columns is
!>        returned in napp. The kernel advances nbasis accordingly.
!> \param c0 MOS holding the deflation columns
!> \param ndefl number of leading c0 columns to deflate against
!> \param skip_defl skip deflation when the seed block is already orthogonal
!>        to the deflation prefix (first block after a full-width RR)
!> \param ws kernel workspace: cand_blk holds the candidate directions, nmob
!>        of them live. vmat is the working basis appended in place.
!>        rank_scale carries the per-column rank-safeguard scale, nrm the
!>        column norms, proj_coeff the projection coefficients, defl_coeff the
!>        deflation-projection coefficients and one_col the scratch column
!>        of the within-batch orthogonalization
!> \param nmob number of live candidate columns
!> \param nbasis number of basis columns on entry, advanced on exit
!> \param napp number of appended columns
! **************************************************************************************************
   SUBROUTINE purify_and_append(c0, ndefl, skip_defl, ws, nmob, nbasis, napp)
      TYPE(cp_cfm_type), INTENT(IN)                      :: c0
      INTEGER, INTENT(IN)                                :: ndefl
      LOGICAL, INTENT(IN)                                :: skip_defl
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws
      INTEGER, INTENT(IN)                                :: nmob
      INTEGER, INTENT(INOUT)                             :: nbasis
      INTEGER, INTENT(OUT)                               :: napp

      REAL(dp), PARAMETER :: eta_dgks = 0.70710678118654752_dp

      COMPLEX(dp)                                        :: col_scale1(1)
      INTEGER                                            :: icol, info, ipass, iqr, kept, nao
      LOGICAL                                            :: chol_ok, do_defl, do_purge

      CALL cp_cfm_get_info(ws%cand_blk, nrow_global=nao)

      ! batch purge: two passes of block gemms project the candidate set
      ! against the working basis and the deflation prefix; the first seed
      ! block after a full-width rotation skips it, later seed blocks must
      ! deflate against the columns earlier blocks wrote
      do_defl = (ndefl > 0) .AND. (.NOT. skip_defl)
      do_purge = (nbasis > 0) .OR. do_defl
      IF (do_purge) THEN
         DO ipass = 1, 2
            IF (nbasis > 0) THEN
               CALL cp_cfm_gemm('C', 'N', nbasis, nmob, nao, z_one, ws%vmat, ws%cand_blk, z_zero, ws%proj_coeff)
               CALL cp_cfm_gemm('N', 'N', nao, nmob, nbasis, z_mone, ws%vmat, ws%proj_coeff, z_one, ws%cand_blk)
            END IF
            IF (do_defl) THEN
               CALL cp_cfm_gemm('C', 'N', ndefl, nmob, nao, z_one, c0, ws%cand_blk, z_zero, ws%defl_coeff)
               CALL cp_cfm_gemm('N', 'N', nao, nmob, ndefl, z_mone, c0, ws%defl_coeff, z_one, ws%cand_blk)
            END IF
            CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)
            IF (ipass == 1) THEN
               ! DGKS norm-ratio test (Daniel, Gragg, Kaufman and Stewart): a
               ! first pass that left every candidate norm within eta_dgks of
               ! its pre-purge norm already orthogonalized the batch to the
               ! projected subspace to working precision
               IF (ALL(ws%nrm(1:nmob) >= eta_dgks*ws%rank_scale(1:nmob))) EXIT
            END IF
         END DO
      END IF

      ! rank safeguard: once the Krylov space saturates the available
      ! dimensions, the orthogonal complement is pure cancellation noise
      ! and a residual of 1e-16 passes a plain positivity test
      IF (ALL(ws%nrm(1:nmob) <= 1.0E-12_dp*ws%rank_scale(1:nmob))) THEN
         napp = 0
         RETURN
      END IF

      ! Path A: single-column fast path (norm guard already satisfied above)
      IF (nmob == 1) THEN
         ws%eval_scale(1) = CMPLX(1.0_dp/ws%nrm(1), 0.0_dp, KIND=dp)
         IF (SIZE(ws%eval_scale) > 1) ws%eval_scale(2:) = z_zero
         CALL cp_cfm_column_scale(ws%cand_blk, ws%eval_scale)
         CALL cp_cfm_to_cfm(ws%cand_blk, ws%vmat, 1, 1, nbasis + 1)
         napp = 1
         nbasis = nbasis + 1
         RETURN
      END IF

      ! Path B: healthy multi-column batch fast path via two-pass BLAS-3 CholeskyQR2
      IF (ALL(ws%nrm(1:nmob) > 1.0E-12_dp*ws%rank_scale(1:nmob))) THEN
         CALL cp_cfm_to_cfm(ws%cand_blk, ws%vmat, nmob, 1, nbasis + 1)
         chol_ok = .TRUE.
         DO iqr = 1, 2
            CALL cp_cfm_gemm('C', 'N', nmob, nmob, nao, z_one, ws%cand_blk, ws%cand_blk, &
                             z_zero, ws%batch_gram)
            CALL cp_cfm_cholesky_decompose(ws%batch_gram, nmob, info)
            IF (info /= 0) THEN
               chol_ok = .FALSE.
               EXIT
            END IF
            CALL cp_cfm_triangular_multiply(ws%batch_gram, ws%cand_blk, side='R', &
                                            invert_tr=.TRUE., n_rows=nao, n_cols=nmob)
         END DO
         IF (chol_ok) THEN
            CALL cp_cfm_to_cfm(ws%cand_blk, ws%vmat, nmob, 1, nbasis + 1)
            napp = nmob
            nbasis = nbasis + nmob
            RETURN
         END IF
         CALL cp_cfm_to_cfm(ws%vmat, ws%cand_blk, nmob, nbasis + 1, 1)
      END IF

      ! Path C: rank-deficient / near-dependent fallback by two-pass CGS2
      ! with per-column norm retirement
      kept = 0
      DO icol = 1, nmob
         CALL cp_cfm_to_cfm(ws%cand_blk, ws%one_col, 1, icol, 1)
         DO ipass = 1, 2
            IF (kept > 0) THEN
               CALL cp_cfm_gemm('C', 'N', kept, 1, nao, z_one, ws%cand_blk, &
                                ws%one_col, z_zero, ws%proj_coeff)
               CALL cp_cfm_gemm('N', 'N', nao, 1, kept, z_mone, ws%cand_blk, &
                                ws%proj_coeff, z_one, ws%one_col)
            END IF
         END DO
         CALL cp_cfm_vectorsnorm(ws%one_col, ws%nrm)
         IF (ws%nrm(1) > 1.0E-12_dp*ws%rank_scale(icol)) THEN
            col_scale1(1) = CMPLX(1.0_dp/ws%nrm(1), 0.0_dp, KIND=dp)
            CALL cp_cfm_column_scale(ws%one_col, col_scale1)
            kept = kept + 1
            CALL cp_cfm_to_cfm(ws%one_col, ws%cand_blk, 1, 1, kept)
         END IF
      END DO
      napp = kept
      IF (napp > 0) CALL cp_cfm_to_cfm(ws%cand_blk, ws%vmat, napp, 1, nbasis + 1)
      nbasis = nbasis + napp
   END SUBROUTINE purify_and_append

! **************************************************************************************************
!> \brief One Rayleigh-Ritz pass on a precomputed basis pair: project
!>        the operator onto the basis, solve the small Hermitian
!>        problem, rotate the leading Ritz vectors into the target
!>        matrix window, and return their exact residual norms.
!> \param abas A-images of the basis columns, A*B
!> \param basis orthonormal basis B; updated in place when c_dst is absent
!> \param nbas live basis width
!> \param npairs number of leading Ritz pairs to write back and form residuals
!> \param ritz_val eigenvalues of the projected problem, at least nbas
!> \param small_t scratch holding the projected matrix on entry, right
!>        sized to nbas
!> \param small_z scratch receiving the eigenvectors, right sized to nbas
!> \param eval_scale complex scaling array spanning every global column
!>        of work_blk
!> \param work_blk scratch receiving the rotated Ritz vectors and then
!>        the residual columns
!> \param dst_col starting destination column in c_dst (or basis)
!> \param res_norm residual norms, spans every global column of work_blk
!> \param c_dst optional target matrix receiving the rotated Ritz vectors
!> \param bas_col optional starting column of basis (defaults to 1)
! **************************************************************************************************
   SUBROUTINE rayleigh_ritz(abas, basis, nbas, npairs, ritz_val, small_t, small_z, &
                            eval_scale, work_blk, dst_col, res_norm, c_dst, bas_col)
      TYPE(cp_cfm_type), INTENT(IN)                      :: abas
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: basis
      INTEGER, INTENT(IN)                                :: nbas, npairs
      REAL(dp), DIMENSION(:)                             :: ritz_val
      TYPE(cp_cfm_type)                                  :: small_t, small_z
      COMPLEX(dp), DIMENSION(:)                          :: eval_scale
      TYPE(cp_cfm_type)                                  :: work_blk
      INTEGER, INTENT(IN)                                :: dst_col
      REAL(dp), DIMENSION(:)                             :: res_norm
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: c_dst
      INTEGER, INTENT(IN), OPTIONAL                      :: bas_col

      INTEGER                                            :: my_bas_col, nao, ncol_work

      my_bas_col = 1
      IF (PRESENT(bas_col)) my_bas_col = bas_col

      CALL cp_cfm_get_info(work_blk, nrow_global=nao, ncol_global=ncol_work)

      CALL cp_cfm_gemm('C', 'N', nbas, nbas, nao, z_one, basis, abas, z_zero, small_t, &
                       a_first_col=my_bas_col)
      CALL cp_cfm_heevd(small_t, small_z, ritz_val(1:nbas))

      CALL cp_cfm_gemm('N', 'N', nao, npairs, nbas, z_one, basis, small_z, z_zero, work_blk, &
                       a_first_col=my_bas_col)
      IF (PRESENT(c_dst)) THEN
         CALL cp_cfm_to_cfm(work_blk, c_dst, npairs, 1, dst_col)
      ELSE
         CALL cp_cfm_to_cfm(work_blk, basis, npairs, 1, dst_col)
      END IF

      eval_scale(1:npairs) = CMPLX(-ritz_val(1:npairs), 0.0_dp, KIND=dp)
      IF (npairs < ncol_work) eval_scale(npairs + 1:ncol_work) = z_zero
      CALL cp_cfm_column_scale(work_blk, eval_scale)
      CALL cp_cfm_gemm('N', 'N', nao, npairs, nbas, z_one, abas, small_z, z_one, work_blk)
      CALL cp_cfm_vectorsnorm(work_blk, res_norm)
   END SUBROUTINE rayleigh_ritz

! **************************************************************************************************
!> \brief Create the kernel workspace at the worst-case widths of the run.
!>        The basis budget holds the larger of the NKRYLOV block budget
!>        and four times the widest single block, and every per-call
!>        clamp of the kernel stays below it. A single block covers
!>        every pending column only when NBLOCK spans them, so four
!>        times MIN(NBLOCK, nmo) bounds that regime. The candidate and
!>        Ritz widths hold the widest block a run can form
!> \param ws workspace to create
!> \param mo_struct structure of the channel MOS, nao x nmo
!> \param ao_struct structure of the channel operators, nao x nao
!> \param nkrylov number of block expansion steps of the block regime
!> \param nblock number of columns refined jointly
! **************************************************************************************************
   SUBROUTINE lanczos_workspace_create(ws, mo_struct, ao_struct, nkrylov, nblock)
      TYPE(lanczos_workspace_type), INTENT(OUT)          :: ws
      TYPE(cp_fm_struct_type), POINTER                   :: mo_struct, ao_struct
      INTEGER, INTENT(IN)                                :: nkrylov, nblock

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_workspace_create'

      INTEGER                                            :: handle, nao, nbe, nmo, wdth
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp

      CALL timeset(routineN, handle)
      NULLIFY (fm_struct_tmp)

      nao = ao_struct%nrow_global
      nmo = mo_struct%ncol_global
      nbe = MIN(MAX(nblock, 1), nmo)
      ! a single block covering every pending column does not multiply
      ! the budget by the block count
      IF (MAX(nblock, 1) >= nmo) THEN
         wdth = MAX(shared_budget_factor*nbe, nkrylov, nbe + 2)
      ELSE
         wdth = MAX(nkrylov*nbe, shared_budget_factor*nbe, nbe + 2)
      END IF

      CALL cp_cfm_create(ws%hmo, mo_struct, name="lan_hmo")
      CALL cp_cfm_create(ws%rot_scr, mo_struct, name="lan_rot_scr")

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nmo, &
                               para_env=mo_struct%para_env, context=mo_struct%context)
      CALL cp_cfm_create(ws%mo_gram, matrix_struct=fm_struct_tmp, name="lan_mo_gram")
      CALL cp_cfm_create(ws%evec, matrix_struct=fm_struct_tmp, name="lan_evec")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=wdth, &
                               para_env=ao_struct%para_env, context=ao_struct%context)
      CALL cp_cfm_create(ws%vmat, matrix_struct=fm_struct_tmp, name="lan_v", set_zero=.TRUE.)
      CALL cp_cfm_create(ws%hvmat, matrix_struct=fm_struct_tmp, name="lan_hv", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nbe, &
                               para_env=ao_struct%para_env, context=ao_struct%context)
      CALL cp_cfm_create(ws%cand_blk, matrix_struct=fm_struct_tmp, name="lan_r", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=1, &
                               para_env=ao_struct%para_env, context=ao_struct%context)
      CALL cp_cfm_create(ws%one_col, matrix_struct=fm_struct_tmp, name="lan_w", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbe, ncol_global=nbe, &
                               para_env=ao_struct%para_env, context=ao_struct%context)
      CALL cp_cfm_create(ws%batch_gram, matrix_struct=fm_struct_tmp, name="lan_bgram", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=wdth, ncol_global=nbe, &
                               para_env=ao_struct%para_env, context=ao_struct%context)
      CALL cp_cfm_create(ws%proj_coeff, matrix_struct=fm_struct_tmp, name="lan_pcoeff", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nbe, &
                               para_env=mo_struct%para_env, context=mo_struct%context)
      CALL cp_cfm_create(ws%defl_coeff, matrix_struct=fm_struct_tmp, name="lan_pcb", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      ALLOCATE (ws%theta(wdth), ws%nrm(nbe), ws%rank_scale(nbe), ws%pair_res(nbe), ws%entry_res(nmo))
      ALLOCATE (ws%ldone(nbe))
      ALLOCATE (ws%eval_scale(MAX(nmo, nbe)))

      CALL timestop(handle)
   END SUBROUTINE lanczos_workspace_create

! **************************************************************************************************
!> \brief Release the kernel workspace. Safe on a workspace that was
!>        never created: the matrix components default to disassociated
!>        and the arrays to unallocated
!> \param ws workspace to release
! **************************************************************************************************
   SUBROUTINE lanczos_workspace_release(ws)
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws

      CALL cp_cfm_release(ws%hmo)
      CALL cp_cfm_release(ws%mo_gram)
      CALL cp_cfm_release(ws%evec)
      CALL cp_cfm_release(ws%rot_scr)
      CALL cp_cfm_release(ws%vmat)
      CALL cp_cfm_release(ws%hvmat)
      CALL cp_cfm_release(ws%cand_blk)
      CALL cp_cfm_release(ws%one_col)
      CALL cp_cfm_release(ws%batch_gram)
      CALL cp_cfm_release(ws%proj_coeff)
      CALL cp_cfm_release(ws%defl_coeff)
      IF (ALLOCATED(ws%theta)) DEALLOCATE (ws%theta)
      IF (ALLOCATED(ws%nrm)) DEALLOCATE (ws%nrm)
      IF (ALLOCATED(ws%rank_scale)) DEALLOCATE (ws%rank_scale)
      IF (ALLOCATED(ws%pair_res)) DEALLOCATE (ws%pair_res)
      IF (ALLOCATED(ws%entry_res)) DEALLOCATE (ws%entry_res)
      IF (ALLOCATED(ws%ldone)) DEALLOCATE (ws%ldone)
      IF (ALLOCATED(ws%eval_scale)) DEALLOCATE (ws%eval_scale)
   END SUBROUTINE lanczos_workspace_release

END MODULE qs_scf_lanczos_c

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex block-Lanczos refinement kernel for the k-point SCF
!>        diagonalization. The kernel processes the non-converged
!>        columns in blocks of NBLOCK. Each seed block enters the
!>        working basis after purification, and the basis grows by the
!>        A-images of the last appended batch, with full
!>        reorthogonalization (two passes) against the working basis,
!>        the columns of c0 ahead of the block, and the other columns
!>        of the same block. The budget holds NKRYLOV block expansion
!>        steps in the block regime, and about four times the pending
!>        width in the shared-subspace regime. The kernel rebuilds the
!>        projected matrix explicitly from the stored A-images, because
!>        full reorthogonalization invalidates any three-term
!>        coefficient bookkeeping. The kernel solves the small
!>        Rayleigh-Ritz problem as complex Hermitian. The leading Ritz
!>        pairs overwrite the target columns and refresh their
!>        eigenvalues at every call. The kernel keeps no state of its own:
!>        every scratch matrix and array lives in a caller-owned workspace,
!>        created once per diagonalization call and reused by every sweep.
!>        The workspace matrices are wider than a single call needs, and
!>        every read stays inside the prefix the call wrote, so the reused
!>        scratch needs no zeroing
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

   PUBLIC :: lanczos_refinement_2v_c
   PUBLIC :: lanczos_workspace_create, lanczos_workspace_release

! **************************************************************************************************
!> \brief Caller-owned scratch of the refinement kernel: every matrix and
!>        array that the sweeps of one channel reuse. The driver creates
!>        the workspace once per diagonalization call and hands it to
!>        every kernel call. The widths follow the worst case of the run:
!>        the basis budget holds the larger of the NKRYLOV block budget
!>        and four times the widest single block, and every per-call
!>        clamp stays below it. The small Rayleigh-Ritz pair keeps its
!>        per-block
!>        right-sized allocation, because the solver stack diagonalizes
!>        whole matrices and knows no leading-block argument
! **************************************************************************************************
   TYPE, PUBLIC :: lanczos_workspace_type
      ! entry Rayleigh-Ritz, its rotation and the shared scratch column
      TYPE(cp_cfm_type) :: entry_rblk, hmo, mo_gram, evec, rot_scr
      ! entry-pass Ritz vector copies at the full MO width
      TYPE(cp_cfm_type) :: entry_ritz_vec, entry_ritz_scaled
      ! working basis and its A-images
      TYPE(cp_cfm_type) :: vmat, hvmat
      ! candidate block, single-column fallback, deflation and
      ! projection-coefficient scratch
      TYPE(cp_cfm_type) :: cand_blk, one_col, proj_coeff, defl_coeff
      ! batch Gram scratch of the within-block orthogonalization
      TYPE(cp_cfm_type) :: batch_gram
      REAL(dp), ALLOCATABLE :: theta(:), nrm(:), rank_scale(:), pair_res(:), entry_res(:), gram_eval(:)
      LOGICAL, ALLOCATABLE :: ldone(:)
      COMPLEX(dp), ALLOCATABLE :: ctheta(:), col_scale(:), eval_scale(:), gram_scale(:)
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

      INTEGER                                            :: handle, info, ipass, iwork, n_conv_cols, &
                                                            nao, nmo, nmo_converged, nmo_nc_loc
      LOGICAL                                            :: my_check_moconv_only, refined
      REAL(dp)                                           :: gfloor, max_norm, min_norm

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(c0, nrow_global=nao, ncol_global=nmo)

      nmo_conv = 0
      nmo_nc = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      my_check_moconv_only = .FALSE.
      IF (PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

      ! entry Rayleigh-Ritz pass on the current MOS: A-images, then the
      ! shared projection, solve and exact residuals
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, ks, c0, z_zero, ws%hmo)
      CALL rayleigh_ritz(ws%hmo, c0, nmo, nmo, eval, ws%mo_gram, ws%evec, &
                         ws%entry_ritz_vec, ws%entry_ritz_scaled, ws%eval_scale, &
                         ws%entry_rblk, ws%entry_res)

      ! Rotate the C vectors
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, c0, ws%evec, z_zero, ws%rot_scr)
      CALL cp_cfm_to_cfm(ws%rot_scr, c0, nmo, 1, 1)

      nmo_converged = 0
      nmo_nc_loc = 0
      max_norm = MAXVAL(ws%entry_res)
      min_norm = MINVAL(ws%entry_res)
      DO iwork = 1, nmo
         IF (ws%entry_res(iwork) <= eps_col(iwork)) THEN
            nmo_converged = nmo_converged + 1
         ELSE
            nmo_nc_loc = nmo - nmo_converged
            EXIT
         END IF
      END DO

      nmo_conv = nmo_converged
      nmo_nc = nmo_nc_loc
      max_res_norm = max_norm
      min_res_norm = min_norm

      IF (.NOT. my_check_moconv_only) THEN
         IF (nmo_nc_loc > 0) THEN
            CALL lanczos_refine(ks, c0, eval, nao, nmo, eps_col, nkrylov, nblock, &
                                nmo_converged, nmo_nc_loc, n_conv_cols, max_res_norm, min_res_norm, &
                                refined, ws)
            IF (.NOT. refined) THEN
               ! no refinement pass ran (degenerate Krylov budget, or every
               ! seed block was rank-deficient): keep the entry statistics
               ! instead of reporting zeros for unconverged columns
               max_res_norm = max_norm
               min_res_norm = min_norm
            END IF
            nmo_conv = nmo_converged + n_conv_cols
            nmo_nc = nmo_nc_loc - n_conv_cols
         END IF

         ! full re-orthonormalization of the refined set, defensive: the
         ! refinement leaves an orthonormal, deflated block behind. But
         ! columns beyond the pass's nres (rank-deficient seeds) keep
         ! their entry vectors and mix with the written block. The Gram
         ! of such a full-rank set can be numerically indefinite, and
         ! the Cholesky factorization can fail. The kernel then falls
         ! back to a spectral reorthonormalization (clamp the Gram
         ! spectrum at a floor and apply the inverse square root)
         ! instead of aborting. The pass cannot repair a genuinely
         ! rank-deficient set: a right transformation preserves the
         ! rank, and the clamped inverse leaves a missing direction
         ! near zero norm instead of creating it. The construction
         ! avoids that case: the entry MOS are orthonormal, the written
         ! block is orthonormal and deflated against the leading
         ! prefix, and mixing the two groups loses no direction. The
         ! next sweep also re-purifies the columns from fresh seeds
         CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, ws%mo_gram)
         CALL cp_cfm_cholesky_decompose(ws%mo_gram, nmo, info)
         IF (info == 0) THEN
            CALL cp_cfm_triangular_multiply(ws%mo_gram, c0, 'R', invert_tr=.TRUE.)
         ELSE
            ! spectral fallback: clamp the Gram spectrum at a floor and
            ! apply the inverse square root. A clamped direction keeps a
            ! residual off-diagonal overlap of order sqrt(g/gfloor) after
            ! one pass. A second pass squares it down. Per-column
            ! renormalization is deliberately avoided: it would rescale a
            ! clamped direction back to unit norm and reimport the very
            ! overlap the clamp bounds. Normalization is not
            ! re-orthogonalization
            DO ipass = 1, 2
               CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, ws%mo_gram)
               CALL cp_cfm_heevd(ws%mo_gram, ws%evec, ws%gram_eval)
               gfloor = 1.0E-10_dp*MAX(MAXVAL(ABS(ws%gram_eval)), TINY(1.0_dp))
               ws%gram_scale(1:nmo) = CMPLX(1.0_dp/SQRT(MAX(ws%gram_eval, gfloor)), 0.0_dp, KIND=dp)
               CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, c0, ws%evec, z_zero, ws%rot_scr)
               CALL cp_cfm_column_scale(ws%rot_scr, ws%gram_scale)
               CALL cp_cfm_gemm('N', 'C', nao, nmo, nmo, z_one, ws%rot_scr, ws%evec, z_zero, c0)
            END DO
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE lanczos_refinement_2v_c

! **************************************************************************************************
!> \brief Per-block single-pass Lanczos refinement.
!>        The kernel processes the non-converged columns in blocks of
!>        nblock columns. When one block covers every pending column
!>        (NBLOCK at least the number of pending MOS), the budget holds
!>        about four times the pending width. The working basis then
!>        becomes the single subspace shared by all columns. The basis
!>        grows by the A-images of the last appended batch, with full
!>        reorthogonalization (two passes) against the working basis,
!>        the columns of c0 ahead of the block, and the other columns
!>        of the same block. The kernel rebuilds the projection
!>        T = V^dag(A V) explicitly from the stored A-images HV = A V,
!>        and diagonalizes it as complex Hermitian. The targets of a
!>        block are the leading pairs of the operator deflated by the
!>        columns ahead of it, one pair per pending column of the
!>        block. The pass rewrites the whole target block with the
!>        current Ritz pairs and refreshes the eigenvalues: the
!>        pair-to-column identity is re-judged at every call and never
!>        latched, because the ascending re-sort can shift pair
!>        identities as soon as an expansion captures a lower state than
!>        the subspace held before.
!> \param ks standardized operator, nao x nao
!> \param c0 MOS, refined in place block by block
!> \param eval eigenvalues, updated per converged column
!> \param nao number of atomic orbitals
!> \param nmo number of molecular orbitals
!> \param eps_col per-column convergence thresholds (see lanczos_refinement_2v_c)
!> \param nkrylov number of block expansion steps in the block regime.
!>        The kernel uses it as a floor on the budget only when nblock
!>        covers every pending column
!> \param nblock number of columns refined jointly
!> \param nmo_converged number of leading locked columns
!> \param nmo_nc_loc number of non-converged tail columns
!> \param n_conv_cols length of the leading converged run across all
!>        blocks: once a block ends short of full convergence, later
!>        blocks no longer contribute even when all their own columns
!>        converged (their columns sit beyond the gap)
!> \param max_res_norm largest residual norm among the wanted pairs of
!>        the pass of each block (exit state, not worst-seen)
!> \param min_res_norm smallest such residual norm. HUGE when no
!>        refinement pass ran at all
!> \param refined whether any block completed a Rayleigh-Ritz pass
!> \param ws caller-owned scratch of the call. Only the small
!>        Rayleigh-Ritz pair stays local to a block, at its right size,
!>        because the solver stack diagonalizes whole matrices
! **************************************************************************************************
   SUBROUTINE lanczos_refine(ks, c0, eval, nao, nmo, eps_col, nkrylov, nblock, &
                             nmo_converged, nmo_nc_loc, n_conv_cols, max_res_norm, min_res_norm, &
                             refined, ws)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      INTEGER, INTENT(IN)                                :: nao, nmo
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock, nmo_converged, &
                                                            nmo_nc_loc
      INTEGER, INTENT(OUT)                               :: n_conv_cols
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      LOGICAL, INTENT(OUT)                               :: refined
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lanczos_refine'

      INTEGER                                            :: col_hi, col_lo, handle, ib, j, nbasis, &
                                                            nblock_e, nbmax, nmob_c, nmobt, nres, &
                                                            nseed, nstep
      LOGICAL                                            :: prefix_open, single_block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: esolve, esolve2, hsolve
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
         nbmax = MIN(nao - 2, MAX(4*nmo_nc_loc, nkrylov, nblock_e + 2))
      ELSE
         nbmax = nkrylov*nblock_e
      END IF
      IF (nbmax < 3) THEN
         ! no valid working basis exists (minimal basis, or a degenerate
         ! NKRYLOV/NBLOCK combination): leave the columns at their entry
         ! state for the caller to report, instead of aborting the run
         n_conv_cols = 0
         max_res_norm = 0.0_dp
         min_res_norm = HUGE(1.0_dp)
         refined = .FALSE.
         CPWARN("Krylov budget below three columns, Lanczos refinement skipped")
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

      block_loop: DO ib = 1, nmo_nc_loc, nblock_e
         col_lo = nmo_converged + ib
         col_hi = MIN(col_lo + nblock_e - 1, nmo_converged + nmo_nc_loc)
         nmobt = col_hi - col_lo + 1
         ws%ldone(1:nmobt) = .FALSE.
         nres = 0

         ! ---------- seed block: the pending columns, purified and
         ! orthonormalized into the working basis (the clamp is defensive:
         ! the block width is already bounded by the budget above)
         CALL cp_cfm_to_cfm(c0, ws%cand_blk, nmobt, col_lo, 1)
         CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)
         ws%rank_scale(1:nmobt) = ws%nrm(1:nmobt)
         nbasis = 0
         CALL purify_and_append(c0, col_lo - 1, ws, MIN(nmobt, nbmax), nbasis, nseed)
         IF (nseed == 0) THEN
            ! defensive branch, not reachable from orthonormal entry MOS:
            ! no seed survived the purification. No column of this block
            ! converged, and the leading prefix cannot extend past it
            prefix_open = .FALSE.
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
            CALL purify_and_append(c0, col_lo - 1, ws, nmob_c, nbasis, nseed)
            nstep = nseed
         END DO
         ! the last batch has no successor: its A-images must still be
         ! formed explicitly, or the projection reads the stale tail of
         ! the reused product buffer and silently degrades the Ritz pairs
         IF (nstep > 0) THEN
            CALL cp_cfm_gemm('N', 'N', nao, nstep, nao, z_one, ks, ws%vmat, z_zero, ws%hvmat, &
                             b_first_col=nbasis - nstep + 1, c_first_col=nbasis - nstep + 1)
         END IF

         ! ---------- Rayleigh-Ritz: T = V^dag (A V), solved as a
         ! distributed problem over the group. The projection forms as a
         ! distributed product on a right-sized matrix and decomposes
         ! through the shared solver stack (ELPA when the run selected
         ! it, ScaLAPACK otherwise, the serial LAPACK driver on a
         ! one-rank context). Only the eigenvector columns that the
         ! consumers read are gathered back. The earlier design
         ! replicated T, decomposed it on the group root and rebroadcast
         ! the spectrum: the serial tridiagonal reduction then idled the
         ! remaining ranks of the group for most of the kernel time
         NULLIFY (fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbasis, ncol_global=nbasis, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_cfm_create(hsolve, matrix_struct=fm_struct_tmp, name="lan_h", set_zero=.TRUE.)
         CALL cp_cfm_create(esolve, matrix_struct=fm_struct_tmp, name="lan_e", set_zero=.TRUE.)
         CALL cp_cfm_create(esolve2, matrix_struct=fm_struct_tmp, name="lan_e2", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct_tmp)
         nres = MIN(nmobt, nbasis)
         ! one shared Rayleigh-Ritz pass per block, on the stored
         ! A-images of the working basis. The Ritz scratch stays at the
         ! right size of the block, so every copy of the eigenvectors
         ! crosses no layouts
         CALL rayleigh_ritz(ws%hvmat, ws%vmat, nbasis, nres, ws%theta, hsolve, esolve, &
                            esolve, esolve2, ws%ctheta, ws%cand_blk, ws%nrm)
         ws%pair_res(1:nres) = ws%nrm(1:nres)

         ! ---------- write-back: the pass rotates the whole target block
         ! into place and refreshes its eigenvalues. No pair-to-column
         ! identity ever latches across calls: heevd re-sorts ascending,
         ! and an expansion that captures a lower state than the subspace
         ! held before shifts every index above it. This write-back is
         ! variationally safe: the seed block contains the purified entry
         ! columns, and their projected matrix equals diag(entry evals).
         ! Courant-Fischer therefore bounds every refreshed eigenvalue
         ! from above by its entry value. Columns beyond nres (basis
         ! saturated below the pending width, rank-deficient seeds) keep
         ! their entry vectors and fall to the end-of-call fallback
         CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_one, ws%vmat, esolve, z_zero, c0, &
                          c_first_col=col_lo)
         DO j = 1, nres
            eval(col_lo + j - 1) = ws%theta(j)
            ws%ldone(j) = ws%pair_res(j) < eps_col(col_lo + j - 1)
         END DO
         refined = .TRUE.
         CALL cp_cfm_release(hsolve)
         CALL cp_cfm_release(esolve)
         CALL cp_cfm_release(esolve2)

         ! exit statistics: the residual numbers describe the state the
         ! columns are left in (the pass's wanted pairs, not the worst
         ! value seen during the call, and not columns beyond nres, which
         ! keep their entry vectors). The converged count is the global
         ! leading converged run. Once a block ends short of full
         ! convergence, no later column can extend the prefix. Later
         ! blocks then contribute nothing, even when all of their own
         ! columns converged
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

      CALL timestop(handle)
   END SUBROUTINE lanczos_refine

! **************************************************************************************************
!> \brief Purify a block of candidate directions and append the survivors to
!>        the working basis: two full projection passes against the working
!>        basis and the leading ndefl columns of c0, modified Gram-Schmidt
!>        within the block, then a per-column rank safeguard relative to
!>        rank_scale, normalization and append. The count of appended columns
!>        is returned in napp. The kernel advances nbasis accordingly.
!> \param c0 MOS holding the deflation columns
!> \param ndefl number of leading c0 columns to deflate against
!> \param ws kernel workspace: cand_blk holds the candidate directions, nmob
!>        of them live. vmat is the working basis appended in place.
!>        rank_scale carries the per-column rank-safeguard scale, nrm the
!>        column norms, ritz_vec the projection coefficients, defl_coeff the
!>        deflation-projection coefficients, batch_gram the batch Gram and
!>        one_col the scratch column of the column-at-a-time fallback
!> \param nmob number of live candidate columns
!> \param nbasis number of basis columns on entry, advanced on exit
!> \param napp number of appended columns
! **************************************************************************************************
   SUBROUTINE purify_and_append(c0, ndefl, ws, nmob, nbasis, napp)
      TYPE(cp_cfm_type), INTENT(IN)                      :: c0
      INTEGER, INTENT(IN)                                :: ndefl
      TYPE(lanczos_workspace_type), INTENT(INOUT)        :: ws
      INTEGER, INTENT(IN)                                :: nmob
      INTEGER, INTENT(INOUT)                             :: nbasis
      INTEGER, INTENT(OUT)                               :: napp

      COMPLEX(dp)                                        :: col_scale1(1)
      INTEGER                                            :: i, icol, info, ipass, mlive, nao, ncolg
      LOGICAL                                            :: lseq

      CALL cp_cfm_get_info(ws%cand_blk, nrow_global=nao, ncol_global=ncolg)

      ! batch purge: two full passes of block gemms project the whole
      ! candidate set against the working basis and the deflation prefix.
      ! This gives the same double-pass quality as the single-column path
      DO ipass = 1, 2
         IF (nbasis > 0) THEN
            CALL cp_cfm_gemm('C', 'N', nbasis, nmob, nao, z_one, ws%vmat, ws%cand_blk, z_zero, ws%proj_coeff)
            CALL cp_cfm_gemm('N', 'N', nao, nmob, nbasis, z_mone, ws%vmat, ws%proj_coeff, z_one, ws%cand_blk)
         END IF
         IF (ndefl > 0) THEN
            CALL cp_cfm_gemm('C', 'N', ndefl, nmob, nao, z_one, c0, ws%cand_blk, z_zero, ws%defl_coeff)
            CALL cp_cfm_gemm('N', 'N', nao, nmob, ndefl, z_mone, c0, ws%defl_coeff, z_one, ws%cand_blk)
         END IF
      END DO
      CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)

      ! rank safeguard: once the Krylov space saturates the available
      ! dimensions, the orthogonal complement is pure cancellation noise.
      ! A residual of 1e-16 passes a plain positivity test, as
      ! He-minimal-basis tests with NKRYLOV larger than the basis showed
      IF (ALL(ws%nrm(1:nmob) <= 1.0E-12_dp*ws%rank_scale(1:nmob))) THEN
         napp = 0
         RETURN
      END IF

      ! within-batch orthogonalization by CholeskyQR2 (two rounds of a
      ! Cholesky factorization of the batch Gram): near-parallel candidates
      ! pass the per-column norm guard, but they square the condition of
      ! the first round into its roundoff, and the second round removes
      ! that error again. A batch of nearly-spanned candidates (the usual
      ! late-SCF expansion) carries a few exhausted columns. Only those
      ! columns are retired: the healthy remainder is compacted to the
      ! front and keeps the batched treatment. A single deficient
      ! direction therefore cannot demote the whole batch to the
      ! column-at-a-time loop below. Only a factorization that reports
      ! rank deficiency on the healthy set demotes the batch
      mlive = 0
      DO i = 1, nmob
         IF (ws%nrm(i) > 1.0E-12_dp*ws%rank_scale(i)) THEN
            mlive = mlive + 1
            IF (mlive /= i) THEN
               CALL cp_cfm_to_cfm(ws%cand_blk, ws%one_col, 1, i, 1)
               CALL cp_cfm_to_cfm(ws%one_col, ws%cand_blk, 1, 1, mlive)
               ws%rank_scale(mlive) = ws%rank_scale(i)
            END IF
         END IF
      END DO
      lseq = .FALSE.
      DO ipass = 1, 2
         CALL cp_cfm_gemm('C', 'N', mlive, mlive, nao, z_one, ws%cand_blk, ws%cand_blk, z_zero, ws%batch_gram)
         CALL cp_cfm_cholesky_decompose(ws%batch_gram, mlive, info)
         IF (info /= 0) THEN
            lseq = .TRUE.
            EXIT
         END IF
         ! n_cols restricts the transform to the live columns: the Gram
         ! scratch is wider than the batch and its stale triangle must
         ! not touch the candidates
         CALL cp_cfm_triangular_multiply(ws%batch_gram, ws%cand_blk, "R", invert_tr=.TRUE., n_cols=mlive)
      END DO
      IF (.NOT. lseq) THEN
         CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)
         ! the scaling array spans every global column of cand_blk: the
         ! parallel branch of cp_cfm_column_scale indexes it with the
         ! global column id. A shorter array reads out of bounds on
         ! ranks owning columns above the live count. The retired tail
         ! scales by one
         ws%col_scale(1:ncolg) = z_one
         DO i = 1, mlive
            ws%col_scale(i) = CMPLX(1.0_dp/ws%nrm(i), 0.0_dp, KIND=dp)
         END DO
         CALL cp_cfm_column_scale(ws%cand_blk, ws%col_scale)
         napp = mlive
         CALL cp_cfm_to_cfm(ws%cand_blk, ws%vmat, mlive, 1, nbasis + 1)
         nbasis = nbasis + napp
         RETURN
      END IF

      ! column-at-a-time path: each purge sees the basis columns of the
      ! previous candidates, which makes the within-batch
      ! orthogonalization a byproduct. A surviving first CholeskyQR2
      ! round already renormalized the batch. The rank test must
      ! therefore reference the current norms, not the pre-purge scales
      ! the caller passed
      CALL cp_cfm_vectorsnorm(ws%cand_blk, ws%nrm)
      ws%rank_scale(1:mlive) = ws%nrm(1:mlive)
      napp = 0
      DO icol = 1, mlive
         CALL cp_cfm_to_cfm(ws%cand_blk, ws%one_col, 1, icol, 1)
         DO ipass = 1, 2
            IF (nbasis + napp > 0) THEN
               CALL cp_cfm_gemm('C', 'N', nbasis + napp, 1, nao, z_one, ws%vmat, ws%one_col, z_zero, ws%proj_coeff)
               CALL cp_cfm_gemm('N', 'N', nao, 1, nbasis + napp, z_mone, ws%vmat, ws%proj_coeff, z_one, ws%one_col)
            END IF
            IF (ndefl > 0) THEN
               CALL cp_cfm_gemm('C', 'N', ndefl, 1, nao, z_one, c0, ws%one_col, z_zero, ws%defl_coeff)
               CALL cp_cfm_gemm('N', 'N', nao, 1, ndefl, z_mone, c0, ws%defl_coeff, z_one, ws%one_col)
            END IF
         END DO
         CALL cp_cfm_vectorsnorm(ws%one_col, ws%nrm)
         IF (ws%nrm(1) > 1.0E-12_dp*ws%rank_scale(icol)) THEN
            col_scale1(1) = CMPLX(1.0_dp/ws%nrm(1), 0.0_dp, KIND=dp)
            CALL cp_cfm_column_scale(ws%one_col, col_scale1)
            CALL cp_cfm_to_cfm(ws%one_col, ws%vmat, 1, 1, nbasis + napp + 1)
            napp = napp + 1
         END IF
      END DO
      nbasis = nbasis + napp
   END SUBROUTINE purify_and_append

! **************************************************************************************************
!> \brief One Rayleigh-Ritz pass on a precomputed basis pair: project
!>        the operator onto the basis, solve the small Hermitian
!>        problem, and return the exact residual norms of the leading
!>        Ritz pairs. The entry pass and every refinement block call
!>        this routine. The scratch is caller-owned, because the entry
!>        works at the full MO width while the refinement blocks stay
!>        at their budget width
!> \param abas A-images of the basis columns, A*B
!> \param basis orthonormal basis B
!> \param nbas live basis width
!> \param npairs number of leading Ritz pairs whose residuals are formed
!> \param ritz_val eigenvalues of the projected problem, at least nbas
!> \param small_t scratch holding the projected matrix on entry, right
!>        sized to nbas
!> \param small_z scratch receiving the eigenvectors, right sized to
!>        nbas. The caller may read it after the pass
!> \param ritz_vec receives the leading npairs eigenvector columns, in
!>        basis coordinates, for the caller's write-back. Same shape as
!>        small_z: the copy between them crosses no layouts
!> \param ritz_scaled receives the eigenvalue-scaled copy of those
!>        columns. Same shape as small_z
!> \param theta_scale complex scaling array spanning every global column
!>        of ritz_scaled. The dead tail beyond npairs scales by one
!> \param res_blk scratch receiving the residual columns
!> \param res_norm residual norms, spans every global column of res_blk
! **************************************************************************************************
   SUBROUTINE rayleigh_ritz(abas, basis, nbas, npairs, ritz_val, small_t, small_z, &
                            ritz_vec, ritz_scaled, theta_scale, res_blk, res_norm)
      TYPE(cp_cfm_type), INTENT(IN)                      :: abas, basis
      INTEGER, INTENT(IN)                                :: nbas, npairs
      REAL(dp), DIMENSION(:)                             :: ritz_val
      TYPE(cp_cfm_type)                                  :: small_t, small_z, ritz_vec, ritz_scaled
      COMPLEX(dp), DIMENSION(:)                          :: theta_scale
      TYPE(cp_cfm_type)                                  :: res_blk
      REAL(dp), DIMENSION(:)                             :: res_norm

      INTEGER                                            :: nao

      CALL cp_cfm_get_info(abas, nrow_global=nao)

      CALL cp_cfm_gemm('C', 'N', nbas, nbas, nao, z_one, basis, abas, z_zero, small_t)
      CALL cp_cfm_heevd(small_t, small_z, ritz_val(1:nbas))

      CALL cp_cfm_to_cfm(small_z, ritz_vec, npairs, 1, 1)
      CALL cp_cfm_to_cfm(small_z, ritz_scaled, npairs, 1, 1)
      theta_scale(:) = z_one
      theta_scale(1:npairs) = CMPLX(ritz_val(1:npairs), 0.0_dp, KIND=dp)
      CALL cp_cfm_column_scale(ritz_scaled, theta_scale)

      CALL cp_cfm_gemm('N', 'N', nao, npairs, nbas, z_one, abas, ritz_vec, z_zero, res_blk)
      CALL cp_cfm_gemm('N', 'N', nao, npairs, nbas, z_mone, basis, ritz_scaled, z_one, res_blk)
      CALL cp_cfm_vectorsnorm(res_blk, res_norm)
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
      wdth = MAX(nkrylov*MAX(nblock, 1), 4*MIN(MAX(nblock, 1), nmo))
      nbe = MIN(MAX(nblock, 1), nmo)

      CALL cp_cfm_create(ws%entry_rblk, mo_struct, name="lan_ctmp")
      CALL cp_cfm_create(ws%hmo, mo_struct, name="lan_hmo")
      CALL cp_cfm_create(ws%rot_scr, mo_struct, name="lan_rot_scr")

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nmo, &
                               para_env=mo_struct%para_env, context=mo_struct%context)
      CALL cp_cfm_create(ws%mo_gram, matrix_struct=fm_struct_tmp, name="lan_mo_gram")
      CALL cp_cfm_create(ws%evec, matrix_struct=fm_struct_tmp, name="lan_evec")
      CALL cp_cfm_create(ws%entry_ritz_vec, matrix_struct=fm_struct_tmp, name="lan_ev1")
      CALL cp_cfm_create(ws%entry_ritz_scaled, matrix_struct=fm_struct_tmp, name="lan_ev2")
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

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=wdth, ncol_global=wdth, &
                               para_env=ao_struct%para_env, context=ao_struct%context)
      CALL cp_cfm_create(ws%batch_gram, matrix_struct=fm_struct_tmp, name="lan_t", set_zero=.TRUE.)
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
      ALLOCATE (ws%ctheta(wdth), ws%col_scale(nbe))
      ALLOCATE (ws%eval_scale(nmo), ws%gram_eval(nmo), ws%gram_scale(nmo))

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

      CALL cp_cfm_release(ws%entry_rblk)
      CALL cp_cfm_release(ws%hmo)
      CALL cp_cfm_release(ws%mo_gram)
      CALL cp_cfm_release(ws%evec)
      CALL cp_cfm_release(ws%entry_ritz_vec)
      CALL cp_cfm_release(ws%entry_ritz_scaled)
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
      IF (ALLOCATED(ws%ctheta)) DEALLOCATE (ws%ctheta)
      IF (ALLOCATED(ws%col_scale)) DEALLOCATE (ws%col_scale)
      IF (ALLOCATED(ws%eval_scale)) DEALLOCATE (ws%eval_scale)
      IF (ALLOCATED(ws%gram_eval)) DEALLOCATE (ws%gram_eval)
      IF (ALLOCATED(ws%gram_scale)) DEALLOCATE (ws%gram_scale)
   END SUBROUTINE lanczos_workspace_release

END MODULE qs_scf_lanczos_c

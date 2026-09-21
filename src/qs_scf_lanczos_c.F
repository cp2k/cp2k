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
!>        eigenvalues at every call. All scratch matrices are local to
!>        the call, so the kernel stays free of CP2K state and
!>        independently testable.
! **************************************************************************************************
MODULE qs_scf_lanczos_c
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_triangular_multiply,&
                                              cp_cfm_vectorsnorm
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_submatrix,&
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

CONTAINS

! **************************************************************************************************
!> \brief Refine the non-converged columns of c0 on the standardized operator ks.
!>        The entry Rayleigh-Ritz pass updates the convergence statistics.
!>        Every non-converged column block then runs a single-pass
!>        Lanczos refinement. A full Cholesky sweep closes the call.
!> \param ks standardized operator A = U^-dag H U^-1, nao x nao
!> \param c0 current MOS in the standardized basis, nao x nmo, refined in place
!> \param c1 scratch of the same shape as c0
!> \param eval eigenvalues of c0, updated in place by the entry Rayleigh-Ritz
!> \param nao number of atomic orbitals
!> \param nmo number of molecular orbitals
!> \param eps_col per-column convergence thresholds on the residual norms:
!>        the driver assigns the relaxed empty threshold to the unoccupied
!>        columns. The occupation decides, not the column index. This
!>        mirrors the Davidson k-point driver
!> \param nkrylov number of block expansion steps in the block regime.
!>        The kernel uses it as a floor on the budget only when NBLOCK
!>        covers every pending column
!> \param nblock number of columns refined jointly, one or more
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
   SUBROUTINE lanczos_refinement_2v_c(ks, c0, c1, eval, nao, nmo, eps_col, nkrylov, nblock, &
                                      max_res_norm, min_res_norm, nmo_conv, nmo_nc, check_moconv_only)

      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0, c1
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      INTEGER, INTENT(IN)                                :: nao, nmo
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      INTEGER, INTENT(OUT)                               :: nmo_conv, nmo_nc
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_moconv_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_refinement_2v_c'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: c_eval_c, g_scale
      INTEGER                                            :: handle, info, ipass, iwork, n_conv_cols, &
                                                            nmo_converged, nmo_nc_loc
      LOGICAL                                            :: my_check_moconv_only, refined
      REAL(dp)                                           :: gfloor, max_norm, min_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: c_res, g_eval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: c_tmp, chc, evec, hc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct_tmp, blacs_env, para_env)

      para_env => c0%matrix_struct%para_env
      blacs_env => c0%matrix_struct%context

      nmo_conv = 0
      nmo_nc = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      my_check_moconv_only = .FALSE.
      IF (PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_cfm_create(c_tmp, matrix_struct=fm_struct_tmp, name="c_tmp")
      CALL cp_cfm_create(hc, matrix_struct=fm_struct_tmp, name="hc")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_cfm_create(chc, matrix_struct=fm_struct_tmp, name="chc")
      CALL cp_cfm_create(evec, matrix_struct=fm_struct_tmp, name="evec")
      CALL cp_fm_struct_release(fm_struct_tmp)

      ! Compute (C^dag)HC
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, ks, c0, z_zero, hc)
      CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, hc, z_zero, chc)

      ! Diagonalize (C^dag)HC
      CALL cp_cfm_heevd(chc, evec, eval)

      ! Rotate the C vectors
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, c0, evec, z_zero, c1)

      ! Check for converged states: residual = H*C*Y - C*Y*diag(eval)
      CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, hc, evec, z_zero, c_tmp)
      CALL cp_cfm_to_cfm(c1, c0, nmo, 1, 1)
      ALLOCATE (c_eval_c(nmo))
      c_eval_c(:) = CMPLX(eval(1:nmo), 0.0_dp, KIND=dp)
      CALL cp_cfm_column_scale(c1, c_eval_c)
      DEALLOCATE (c_eval_c)
      CALL cp_cfm_scale_and_add(z_one, c_tmp, z_mone, c1)
      ALLOCATE (c_res(nmo))
      CALL cp_cfm_vectorsnorm(c_tmp, c_res)

      nmo_converged = 0
      nmo_nc_loc = 0
      max_norm = MAXVAL(c_res)
      min_norm = MINVAL(c_res)
      DO iwork = 1, nmo
         IF (c_res(iwork) <= eps_col(iwork)) THEN
            nmo_converged = nmo_converged + 1
         ELSE
            nmo_nc_loc = nmo - nmo_converged
            EXIT
         END IF
      END DO
      DEALLOCATE (c_res)

      CALL cp_cfm_release(c_tmp)
      CALL cp_cfm_release(hc)

      nmo_conv = nmo_converged
      nmo_nc = nmo_nc_loc
      max_res_norm = max_norm
      min_res_norm = min_norm

      IF (.NOT. my_check_moconv_only) THEN
         IF (nmo_nc_loc > 0) THEN
            CALL lanczos_refine(ks, c0, eval, nao, nmo, eps_col, nkrylov, nblock, &
                                nmo_converged, nmo_nc_loc, n_conv_cols, max_res_norm, min_res_norm, &
                                refined)
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
         ! of such a set can be numerically indefinite, and the Cholesky
         ! factorization can fail. The kernel then falls back to a
         ! spectral reorthonormalization (clamp the Gram spectrum at a
         ! floor and apply the inverse square root) instead of aborting
         CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, chc)
         CALL cp_cfm_cholesky_decompose(chc, nmo, info)
         IF (info == 0) THEN
            CALL cp_cfm_triangular_multiply(chc, c0, 'R', invert_tr=.TRUE.)
         ELSE
            ! spectral fallback: clamp the Gram spectrum at a floor and
            ! apply the inverse square root. A clamped direction keeps a
            ! residual off-diagonal overlap of order sqrt(g/gfloor) after
            ! one pass. A second pass squares it down. Per-column
            ! renormalization is deliberately avoided: it would rescale a
            ! clamped direction back to unit norm and reimport the very
            ! overlap the clamp bounds. Normalization is not
            ! re-orthogonalization
            ALLOCATE (g_eval(nmo), g_scale(nmo))
            DO ipass = 1, 2
               CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, chc)
               CALL cp_cfm_heevd(chc, evec, g_eval)
               gfloor = 1.0E-10_dp*MAX(MAXVAL(ABS(g_eval)), TINY(1.0_dp))
               g_scale(:) = CMPLX(1.0_dp/SQRT(MAX(g_eval, gfloor)), 0.0_dp, KIND=dp)
               CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, c0, evec, z_zero, c1)
               CALL cp_cfm_column_scale(c1, g_scale)
               CALL cp_cfm_gemm('N', 'C', nao, nmo, nmo, z_one, c1, evec, z_zero, c0)
            END DO
            DEALLOCATE (g_eval, g_scale)
         END IF
      END IF

      CALL cp_cfm_release(chc)
      CALL cp_cfm_release(evec)

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
! **************************************************************************************************
   SUBROUTINE lanczos_refine(ks, c0, eval, nao, nmo, eps_col, nkrylov, nblock, &
                             nmo_converged, nmo_nc_loc, n_conv_cols, max_res_norm, min_res_norm, &
                             refined)
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lanczos_refine'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: ctheta
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: zvec
      INTEGER                                            :: col_hi, col_lo, handle, ib, j, nbasis, &
                                                            nblock_e, nbmax, nmob_c, nmobt, nres, &
                                                            nseed, nstep
      LOGICAL                                            :: prefix_open, single_block
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: ldone
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: gscale, nrm, rhot, theta
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: esolve, hsolve, hvmat, proj_cb, rblk, &
                                                            tmatrix, vmat, wvec, yw1, yw2
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

      ! workspaces, all local to the call
      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nbmax, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(vmat, matrix_struct=fm_struct_tmp, name="lan_v", set_zero=.TRUE.)
      CALL cp_cfm_create(hvmat, matrix_struct=fm_struct_tmp, name="lan_hv", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nblock_e, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(rblk, matrix_struct=fm_struct_tmp, name="lan_r", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=1, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(wvec, matrix_struct=fm_struct_tmp, name="lan_w", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbmax, ncol_global=nbmax, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(tmatrix, matrix_struct=fm_struct_tmp, name="lan_t", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbmax, ncol_global=nblock_e, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(yw1, matrix_struct=fm_struct_tmp, name="lan_yw1", set_zero=.TRUE.)
      CALL cp_cfm_create(yw2, matrix_struct=fm_struct_tmp, name="lan_yw2", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nblock_e, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(proj_cb, matrix_struct=fm_struct_tmp, name="lan_pcb", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      ALLOCATE (theta(nbmax), nrm(nblock_e), gscale(nblock_e), rhot(nblock_e), ldone(nblock_e))
      ALLOCATE (ctheta(nblock_e))
      ALLOCATE (zvec(nbmax, MAX(nmo_nc_loc, 1)))

      n_conv_cols = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      prefix_open = .TRUE.
      refined = .FALSE.

      block_loop: DO ib = 1, nmo_nc_loc, nblock_e
         col_lo = nmo_converged + ib
         col_hi = MIN(col_lo + nblock_e - 1, nmo_converged + nmo_nc_loc)
         nmobt = col_hi - col_lo + 1
         ldone(1:nmobt) = .FALSE.
         nres = 0

         ! ---------- seed block: the pending columns, purified and
         ! orthonormalized into the working basis (the clamp is defensive:
         ! the block width is already bounded by the budget above)
         CALL cp_cfm_to_cfm(c0, rblk, nmobt, col_lo, 1)
         CALL cp_cfm_vectorsnorm(rblk, nrm)
         gscale(1:nmobt) = nrm(1:nmobt)
         nbasis = 0
         CALL purify_and_append(c0, col_lo - 1, rblk, MIN(nmobt, nbmax), gscale, vmat, nbasis, nseed, &
                                yw1, proj_cb, tmatrix, wvec, nrm)
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
            CALL cp_cfm_gemm('N', 'N', nao, nstep, nao, z_one, ks, vmat, z_zero, hvmat, &
                             b_first_col=nbasis - nstep + 1, c_first_col=nbasis - nstep + 1)
            CALL cp_cfm_to_cfm(hvmat, rblk, nstep, nbasis - nstep + 1, 1)
            CALL cp_cfm_vectorsnorm(rblk, nrm)
            gscale(1:nstep) = nrm(1:nstep)
            nmob_c = MIN(nstep, nbmax - nbasis)
            CALL purify_and_append(c0, col_lo - 1, rblk, nmob_c, gscale, vmat, nbasis, nseed, &
                                   yw1, proj_cb, tmatrix, wvec, nrm)
            nstep = nseed
         END DO
         ! the last batch has no successor: its A-images must still be
         ! formed explicitly, or the projection reads the zeroed tail of
         ! the product buffer and silently degrades the Ritz pairs
         IF (nstep > 0) THEN
            CALL cp_cfm_gemm('N', 'N', nao, nstep, nao, z_one, ks, vmat, z_zero, hvmat, &
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
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_cfm_gemm('C', 'N', nbasis, nbasis, nao, z_one, vmat, hvmat, z_zero, hsolve)
         CALL cp_cfm_heevd(hsolve, esolve, theta(1:nbasis))
         nres = MIN(nmobt, nbasis)
         ! only the leading nres pairs, one per pending column, are ever
         ! read, so that block alone is gathered back
         CALL cp_cfm_get_submatrix(esolve, zvec(1:nbasis, 1:nres), 1, 1, nbasis, nres)
         CALL cp_cfm_release(hsolve)
         CALL cp_cfm_release(esolve)

         ! ---------- exact residuals of the target pairs: the leading
         ! Ritz pairs, one per pending column
         CALL cp_cfm_set_submatrix(yw1, zvec(1:nbasis, 1:nres), 1, 1, nbasis, nres)
         CALL cp_cfm_set_submatrix(yw2, zvec(1:nbasis, 1:nres), 1, 1, nbasis, nres)
         ! the scaling array spans every global column of yw2: the
         ! parallel branch of cp_cfm_column_scale indexes it with the
         ! global column id, and a rank owning columns above nres would
         ! read past a shorter array. The dead tail scales by one
         ctheta(:) = z_one
         ctheta(1:nres) = CMPLX(theta(1:nres), 0.0_dp, KIND=dp)
         CALL cp_cfm_column_scale(yw2, ctheta)
         CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_one, hvmat, yw1, z_zero, rblk)
         CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_mone, vmat, yw2, z_one, rblk)
         CALL cp_cfm_vectorsnorm(rblk, nrm)
         rhot(1:nres) = nrm(1:nres)

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
         CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_one, vmat, yw1, z_zero, c0, &
                          c_first_col=col_lo)
         DO j = 1, nres
            eval(col_lo + j - 1) = theta(j)
            ldone(j) = rhot(j) < eps_col(col_lo + j - 1)
         END DO
         refined = .TRUE.

         ! exit statistics: the residual numbers describe the state the
         ! columns are left in (the pass's wanted pairs, not the worst
         ! value seen during the call, and not columns beyond nres, which
         ! keep their entry vectors). The converged count is the global
         ! leading converged run. Once a block ends short of full
         ! convergence, no later column can extend the prefix. Later
         ! blocks then contribute nothing, even when all of their own
         ! columns converged
         DO j = 1, nres
            max_res_norm = MAX(max_res_norm, rhot(j))
            min_res_norm = MIN(min_res_norm, rhot(j))
         END DO
         IF (prefix_open) THEN
            DO j = 1, nmobt
               IF (.NOT. ldone(j)) EXIT
               n_conv_cols = n_conv_cols + 1
            END DO
            prefix_open = ALL(ldone(1:nmobt))
         END IF
      END DO block_loop

      DEALLOCATE (theta, nrm, gscale, rhot, ldone, ctheta, zvec)
      CALL cp_cfm_release(vmat)
      CALL cp_cfm_release(hvmat)
      CALL cp_cfm_release(rblk)
      CALL cp_cfm_release(wvec)
      CALL cp_cfm_release(proj_cb)
      CALL cp_cfm_release(tmatrix)
      CALL cp_cfm_release(yw1)
      CALL cp_cfm_release(yw2)

      CALL timestop(handle)
   END SUBROUTINE lanczos_refine

! **************************************************************************************************
!> \brief Purify a block of candidate directions and append the survivors to
!>        the working basis: two full projection passes against the working
!>        basis and the leading ndefl columns of c0, modified Gram-Schmidt
!>        within the block, then a per-column rank safeguard relative to
!>        gscale, normalization and append. The count of appended columns
!>        is returned in napp. The kernel advances nbasis accordingly.
!> \param c0 MOS holding the deflation columns
!> \param ndefl number of leading c0 columns to deflate against
!> \param rblk candidate directions, nao columns x nblock_e, nmob live
!> \param nmob number of live candidate columns
!> \param gscale per-column rank-safeguard scale (pre-purge norms)
!> \param vmat working basis, appended in place
!> \param nbasis number of basis columns on entry, advanced on exit
!> \param napp number of appended columns
!> \param coef scratch for the projection coefficients, nbmax x nblock_e
!> \param proj_cb deflation-projection coefficients, nmo x nblock_e
!> \param gram scratch for the batch Gram, nbmax x nbmax
!> \param wvec scratch column of the column-at-a-time fallback
!> \param nrm column norms, at least nblock_e
! **************************************************************************************************
   SUBROUTINE purify_and_append(c0, ndefl, rblk, nmob, gscale, vmat, nbasis, napp, &
                                coef, proj_cb, gram, wvec, nrm)
      TYPE(cp_cfm_type), INTENT(IN)                      :: c0
      INTEGER, INTENT(IN)                                :: ndefl
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: rblk
      INTEGER, INTENT(IN)                                :: nmob
      REAL(dp), DIMENSION(:)                             :: gscale
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: vmat
      INTEGER, INTENT(INOUT)                             :: nbasis
      INTEGER, INTENT(OUT)                               :: napp
      TYPE(cp_cfm_type)                                  :: coef
      TYPE(cp_cfm_type), INTENT(IN)                      :: proj_cb
      TYPE(cp_cfm_type)                                  :: gram
      TYPE(cp_cfm_type), INTENT(IN)                      :: wvec
      REAL(dp), DIMENSION(:)                             :: nrm

      COMPLEX(dp)                                        :: zscl1(1)
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: zscl
      INTEGER                                            :: i, icol, info, ipass, mlive, nao, ncolg
      LOGICAL                                            :: lseq

      CALL cp_cfm_get_info(rblk, nrow_global=nao, ncol_global=ncolg)

      ! batch purge: two full passes of block gemms project the whole
      ! candidate set against the working basis and the deflation prefix.
      ! This gives the same double-pass quality as the single-column path
      DO ipass = 1, 2
         IF (nbasis > 0) THEN
            CALL cp_cfm_gemm('C', 'N', nbasis, nmob, nao, z_one, vmat, rblk, z_zero, coef)
            CALL cp_cfm_gemm('N', 'N', nao, nmob, nbasis, z_mone, vmat, coef, z_one, rblk)
         END IF
         IF (ndefl > 0) THEN
            CALL cp_cfm_gemm('C', 'N', ndefl, nmob, nao, z_one, c0, rblk, z_zero, proj_cb)
            CALL cp_cfm_gemm('N', 'N', nao, nmob, ndefl, z_mone, c0, proj_cb, z_one, rblk)
         END IF
      END DO
      CALL cp_cfm_vectorsnorm(rblk, nrm)

      ! rank safeguard: once the Krylov space saturates the available
      ! dimensions, the orthogonal complement is pure cancellation noise.
      ! A residual of 1e-16 passes a plain positivity test, as
      ! He-minimal-basis tests with NKRYLOV larger than the basis showed
      IF (ALL(nrm(1:nmob) <= 1.0E-12_dp*gscale(1:nmob))) THEN
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
         IF (nrm(i) > 1.0E-12_dp*gscale(i)) THEN
            mlive = mlive + 1
            IF (mlive /= i) THEN
               CALL cp_cfm_to_cfm(rblk, wvec, 1, i, 1)
               CALL cp_cfm_to_cfm(wvec, rblk, 1, 1, mlive)
               gscale(mlive) = gscale(i)
            END IF
         END IF
      END DO
      lseq = .FALSE.
      DO ipass = 1, 2
         CALL cp_cfm_gemm('C', 'N', mlive, mlive, nao, z_one, rblk, rblk, z_zero, gram)
         CALL cp_cfm_cholesky_decompose(gram, mlive, info)
         IF (info /= 0) THEN
            lseq = .TRUE.
            EXIT
         END IF
         ! n_cols restricts the transform to the live columns: the Gram
         ! scratch is wider than the batch and its stale triangle must
         ! not touch the candidates
         CALL cp_cfm_triangular_multiply(gram, rblk, "R", invert_tr=.TRUE., n_cols=mlive)
      END DO
      IF (.NOT. lseq) THEN
         CALL cp_cfm_vectorsnorm(rblk, nrm)
         ! the scaling array spans every global column of rblk: the
         ! parallel branch of cp_cfm_column_scale indexes it with the
         ! global column id. A shorter array reads out of bounds on
         ! ranks owning columns above the live count. The retired tail
         ! scales by one
         ALLOCATE (zscl(ncolg))
         zscl(:) = z_one
         DO i = 1, mlive
            zscl(i) = CMPLX(1.0_dp/nrm(i), 0.0_dp, KIND=dp)
         END DO
         CALL cp_cfm_column_scale(rblk, zscl)
         DEALLOCATE (zscl)
         napp = mlive
         CALL cp_cfm_to_cfm(rblk, vmat, mlive, 1, nbasis + 1)
         nbasis = nbasis + napp
         RETURN
      END IF

      ! column-at-a-time path: each purge sees the basis columns of the
      ! previous candidates, which makes the within-batch
      ! orthogonalization a byproduct. A surviving first CholeskyQR2
      ! round already renormalized the batch. The rank test must
      ! therefore reference the current norms, not the pre-purge scales
      ! the caller passed
      CALL cp_cfm_vectorsnorm(rblk, nrm)
      gscale(1:mlive) = nrm(1:mlive)
      napp = 0
      DO icol = 1, mlive
         CALL cp_cfm_to_cfm(rblk, wvec, 1, icol, 1)
         DO ipass = 1, 2
            IF (nbasis + napp > 0) THEN
               CALL cp_cfm_gemm('C', 'N', nbasis + napp, 1, nao, z_one, vmat, wvec, z_zero, coef)
               CALL cp_cfm_gemm('N', 'N', nao, 1, nbasis + napp, z_mone, vmat, coef, z_one, wvec)
            END IF
            IF (ndefl > 0) THEN
               CALL cp_cfm_gemm('C', 'N', ndefl, 1, nao, z_one, c0, wvec, z_zero, proj_cb)
               CALL cp_cfm_gemm('N', 'N', nao, 1, ndefl, z_mone, c0, proj_cb, z_one, wvec)
            END IF
         END DO
         CALL cp_cfm_vectorsnorm(wvec, nrm)
         IF (nrm(1) > 1.0E-12_dp*gscale(icol)) THEN
            zscl1(1) = CMPLX(1.0_dp/nrm(1), 0.0_dp, KIND=dp)
            CALL cp_cfm_column_scale(wvec, zscl1)
            CALL cp_cfm_to_cfm(wvec, vmat, 1, 1, nbasis + napp + 1)
            napp = napp + 1
         END IF
      END DO
      nbasis = nbasis + napp
   END SUBROUTINE purify_and_append

END MODULE qs_scf_lanczos_c

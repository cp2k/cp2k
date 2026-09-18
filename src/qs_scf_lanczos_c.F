!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex thick-restart block-Lanczos (TRLan) refinement kernel for the k-point
!>        SCF diagonalization, after Wu & Simon, SIAM J. Matrix Anal. Appl. 22,
!>        602 (2000). The non-converged columns are processed in blocks of
!>        NBLOCK: the working Krylov basis of a block persists across restart
!>        cycles, and each restart keeps the Ritz-optimal vectors together
!>        with their A-images, so the subspace accumulates information
!>        instead of freezing onto the current columns. The projected matrix
!>        is rebuilt explicitly from the stored A-images (full
!>        reorthogonalization invalidates any three-term coefficient
!>        bookkeeping), and the small Rayleigh-Ritz problem is solved as
!>        complex Hermitian. All scratch matrices are local to the call, so
!>        the kernel stays free of CP2K state and independently testable.
! **************************************************************************************************
MODULE qs_scf_lanczos_c
   USE bibliography,                    ONLY: WuSimon2000,&
                                              cite_reference
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
!>        The entry Rayleigh-Ritz pass updates the convergence statistics; every
!>        non-converged column block then runs its own thick-restart Lanczos
!>        cycle loop; a full Cholesky sweep closes the call.
!> \param ks standardized operator A = U^-dag H U^-1, nao x nao
!> \param c0 current MOS in the standardized basis, nao x nmo, refined in place
!> \param c1 scratch of the same shape as c0
!> \param eval eigenvalues of c0, updated in place by the entry Rayleigh-Ritz
!> \param nao number of atomic orbitals
!> \param nmo number of molecular orbitals
!> \param eps_col per-column convergence thresholds on the residual norms:
!>        the driver assigns the relaxed empty threshold to the unoccupied
!>        columns (the occupation decides, not the column index), mirroring
!>        the Davidson k-point driver
!> \param nkrylov number of block expansion steps per cycle in the block
!>        regime; a floor on the budget only once NBLOCK covers every
!>        pending column
!> \param nblock number of columns refined jointly, one or more
!> \param max_cycles TRLan cycles allowed per column block
!> \param max_res_norm largest residual column norm on exit: the entry
!>        Rayleigh-Ritz values, or the final-cycle wanted-pair residuals
!>        of the refinement blocks when a refinement ran
!> \param min_res_norm smallest residual column norm on exit (same basis)
!> \param nmo_conv number of leading converged columns on exit
!> \param nmo_nc columns outside the leading converged prefix: every
!>        still-unconverged column, plus converged columns beyond a gap;
!>        conservative for the not-converged reporting
!> \param check_moconv_only only run the entry Rayleigh-Ritz and the
!>        convergence count, skip the refinement
!> \param ncycles number of thick-restart cycles consumed by the
!>        refinement blocks
! **************************************************************************************************
   SUBROUTINE lanczos_refinement_2v_c(ks, c0, c1, eval, nao, nmo, eps_col, nkrylov, nblock, max_cycles, &
                                      max_res_norm, min_res_norm, nmo_conv, nmo_nc, check_moconv_only, ncycles)

      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0, c1
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      INTEGER, INTENT(IN)                                :: nao, nmo
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock, max_cycles
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      INTEGER, INTENT(OUT)                               :: nmo_conv, nmo_nc
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_moconv_only
      INTEGER, INTENT(OUT), OPTIONAL                     :: ncycles

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_refinement_2v_c'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: c_eval_c, g_scale
      INTEGER                                            :: handle, info, ipass, iwork, my_ncycles, &
                                                            n_conv_cols, nmo_converged, nmo_nc_loc
      LOGICAL                                            :: my_check_moconv_only
      REAL(dp)                                           :: gfloor, max_norm, min_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: c_res, g_eval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: c_tmp, chc, evec, hc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL cite_reference(WuSimon2000)
      NULLIFY (fm_struct_tmp, blacs_env, para_env)

      para_env => c0%matrix_struct%para_env
      blacs_env => c0%matrix_struct%context

      nmo_conv = 0
      nmo_nc = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      my_check_moconv_only = .FALSE.
      IF (PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only
      IF (PRESENT(ncycles)) ncycles = 0

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
            CALL trlan_refine(ks, c0, eval, nao, nmo, eps_col, nkrylov, nblock, max_cycles, &
                              nmo_converged, nmo_nc_loc, n_conv_cols, max_res_norm, min_res_norm, &
                              my_ncycles)
            IF (my_ncycles == 0) THEN
               ! no refinement cycle ran (degenerate Krylov budget, or every
               ! seed block was rank-deficient): keep the entry statistics
               ! instead of reporting zeros for unconverged columns
               max_res_norm = max_norm
               min_res_norm = min_norm
            END IF
            IF (PRESENT(ncycles)) ncycles = my_ncycles
            nmo_conv = nmo_converged + n_conv_cols
            nmo_nc = nmo_nc_loc - n_conv_cols
         END IF

         ! full re-orthonormalization of the refined set, defensive: the
         ! cycles leave an orthonormal, deflated block behind, but columns
         ! beyond the last cycle's nres (rank-deficient seeds) keep their
         ! entry vectors and mix with the written block; the Gram of such a
         ! set can be numerically indefinite and the Cholesky factorization
         ! can fail; fall back to a spectral reorthonormalization (clamp
         ! the Gram spectrum at a floor and apply the inverse square root)
         ! instead of aborting
         CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, c0, c0, z_zero, chc)
         CALL cp_cfm_cholesky_decompose(chc, nmo, info)
         IF (info == 0) THEN
            CALL cp_cfm_triangular_multiply(chc, c0, 'R', invert_tr=.TRUE.)
         ELSE
            ! spectral fallback: clamp the Gram spectrum at a floor and
            ! apply the inverse square root. A clamped direction keeps a
            ! residual off-diagonal overlap of order sqrt(g/gfloor) after
            ! one pass; a second pass squares it down. Per-column
            ! renormalization is deliberately avoided: it would rescale a
            ! clamped direction back to unit norm and reimport the very
            ! overlap the clamp bounds, and normalization is not
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
!> \brief Per-block thick-restart Lanczos refinement (Wu & Simon 2000).
!>        The non-converged columns are processed in blocks of nblock
!>        columns; when one block covers every pending column (NBLOCK at
!>        least the number of pending MOS), the budget follows the m ~ 2k
!>        convention of the thick-restart Lanczos method and the restart
!>        keeps every target pair, which turns the working basis into the
!>        single subspace shared by all columns. A block runs its own cycle loop: the Krylov basis grows
!>        to nkrylov blocks of columns with full reorthogonalization (two
!>        passes) against the working basis, the columns of c0 ahead of the
!>        block, and the other columns of the same block; the projection
!>        T = V^dag(A V) is rebuilt explicitly from the stored A-images
!>        HV = A V and diagonalized as complex Hermitian. The targets of a
!>        block are the leading pairs of the operator deflated by the
!>        columns ahead of it, one pair per pending column of the block.
!>        Every cycle rewrites the whole target block with the current
!>        Ritz pairs and refreshes the eigenvalues: the pair-to-column
!>        identity is re-judged each cycle and never latched, because the
!>        ascending re-sort can shift pair identities as soon as an
!>        expansion captures a lower state than the subspace held before;
!>        a cycle that is not done restarts keeping the lowest Ritz pairs
!>        (vectors and A-images, floored so no target pair leaves the
!>        basis) and seeding with the residual vectors of the
!>        still-pending pairs.
!> \param ks standardized operator, nao x nao
!> \param c0 MOS, refined in place block by block
!> \param eval eigenvalues, updated per converged column
!> \param nao number of atomic orbitals
!> \param nmo number of molecular orbitals
!> \param eps_col per-column convergence thresholds (see lanczos_refinement_2v_c)
!> \param nkrylov number of block expansion steps per cycle in the block
!>        regime; a floor on the budget only once nblock covers every
!>        pending column
!> \param nblock number of columns refined jointly
!> \param max_cycles TRLan cycles allowed per column block
!> \param nmo_converged number of leading locked columns
!> \param nmo_nc_loc number of non-converged tail columns
!> \param n_conv_cols length of the leading converged run across all
!>        blocks: once a block ends short of full convergence, later
!>        blocks no longer contribute even when all their own columns
!>        converged (their columns sit beyond the gap)
!> \param max_res_norm largest residual norm among the wanted pairs of
!>        the final cycle of each block (exit state, not worst-seen)
!> \param min_res_norm smallest such residual norm; HUGE when no
!>        refinement cycle ran at all
!> \param ncycles number of thick-restart cycles consumed across the
!>        column blocks
! **************************************************************************************************
   SUBROUTINE trlan_refine(ks, c0, eval, nao, nmo, eps_col, nkrylov, nblock, max_cycles, &
                           nmo_converged, nmo_nc_loc, n_conv_cols, max_res_norm, min_res_norm, &
                           ncycles)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ks
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c0
      REAL(dp), DIMENSION(:), POINTER                    :: eval
      INTEGER, INTENT(IN)                                :: nao, nmo
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: eps_col
      INTEGER, INTENT(IN)                                :: nkrylov, nblock, max_cycles, &
                                                            nmo_converged, nmo_nc_loc
      INTEGER, INTENT(OUT)                               :: n_conv_cols
      REAL(dp), INTENT(OUT)                              :: max_res_norm, min_res_norm
      INTEGER, INTENT(OUT)                               :: ncycles

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trlan_refine'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: ctheta
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: zt, zvec
      INTEGER                                            :: col_hi, col_lo, cycle_, handle, ib, j, &
                                                            k_save, nbasis, nblock_e, nbmax, nbuf, &
                                                            nmob_c, nmobt, nres, nseed, nstep
      LOGICAL                                            :: prefix_open, single_block
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: ldone
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: gscale, nrm, rhot, theta
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: esolve, hsolve, hvmat, pblk, proj_cb, &
                                                            rblk, restart_buf, tmatrix, vmat, &
                                                            wvec, yw1, yw2
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
         ! about twice the pending width (the restarting heuristic of Wu
         ! and Simon 2000), floored at NKRYLOV and at the block width, so
         ! the shared subspace holds every target pair plus an equally
         ! wide expansion room; NKRYLOV plays no role here
         nbmax = MIN(nao - 2, MAX(2*nmo_nc_loc, nkrylov, nblock_e + 2))
      ELSE
         nbmax = nkrylov*nblock_e
      END IF
      IF (nbmax < 3) THEN
         ! no valid working basis exists (minimal basis, or a degenerate
         ! NKRYLOV/NBLOCK combination): leave the columns at their entry
         ! state for the caller to report, instead of aborting the run
         n_conv_cols = 0
         ncycles = 0
         max_res_norm = 0.0_dp
         min_res_norm = HUGE(1.0_dp)
         CPWARN("Krylov budget below three columns, Lanczos refinement skipped")
         CALL timestop(handle)
         RETURN
      END IF
      ! the budget bounds the block width: a wider block would append
      ! past the end of vmat at its first seed call (nao - 2 caps nbmax
      ! below the pending width on minimal bases); the pending set splits
      ! into narrower blocks instead, no column is dropped
      IF (nblock_e > nbmax - 2) THEN
         nblock_e = MAX(nbmax - 2, 1)
         single_block = nblock_e >= nmo_nc_loc
      END IF

      ! workspaces, all local to the call
      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nbmax, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(vmat, matrix_struct=fm_struct_tmp, name="trl_v", set_zero=.TRUE.)
      CALL cp_cfm_create(hvmat, matrix_struct=fm_struct_tmp, name="trl_hv", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nblock_e, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(pblk, matrix_struct=fm_struct_tmp, name="trl_p", set_zero=.TRUE.)
      CALL cp_cfm_create(rblk, matrix_struct=fm_struct_tmp, name="trl_r", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=1, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(wvec, matrix_struct=fm_struct_tmp, name="trl_w", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbmax, ncol_global=nbmax, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(tmatrix, matrix_struct=fm_struct_tmp, name="trl_t", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbmax, ncol_global=nblock_e, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(yw1, matrix_struct=fm_struct_tmp, name="trl_yw1", set_zero=.TRUE.)
      CALL cp_cfm_create(yw2, matrix_struct=fm_struct_tmp, name="trl_yw2", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      ! restart buffer, twice the largest restart width the retention rule
      ! can produce: k_save is capped by half the budget or by the block
      ! width (the retention floor never exceeds the target count, which
      ! the block width bounds), so the buffer covers whichever cap is
      ! larger, and by nbmax - 2 in any case
      nbuf = 2*MIN(nbmax - 2, MAX(MAX(nbmax/2, 1), nblock_e))
      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nbuf, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(restart_buf, matrix_struct=fm_struct_tmp, name="trl_rb", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nblock_e, &
                               para_env=para_env, context=blacs_env)
      CALL cp_cfm_create(proj_cb, matrix_struct=fm_struct_tmp, name="trl_pcb", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct_tmp)

      ALLOCATE (theta(nbmax), nrm(nblock_e), gscale(nblock_e), rhot(nblock_e), ldone(nblock_e))
      ALLOCATE (ctheta(nblock_e))
      ALLOCATE (zt(nbmax, nbmax))
      ALLOCATE (zvec(nbmax, MAX(nmo_nc_loc, 1)))

      n_conv_cols = 0
      ncycles = 0
      max_res_norm = 0.0_dp
      min_res_norm = HUGE(1.0_dp)
      prefix_open = .TRUE.

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
         IF (nseed == 0) CYCLE block_loop
         nstep = nseed

         cycle_loop: DO cycle_ = 1, max_cycles
            ncycles = ncycles + 1

            ! ---------- expansion: extend the basis with the A-images of the
            ! last block until the budget is exhausted
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
            ! the last block has no successor: its A-images must still be
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
            ! one-rank context); only the eigenvector columns the consumers
            ! read are gathered back. The earlier design replicated T,
            ! decomposed it on the group root and rebroadcast the spectrum:
            ! the serial tridiagonal reduction then idled the remaining
            ! ranks of the group for most of the kernel time
            NULLIFY (fm_struct_tmp)
            CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nbasis, ncol_global=nbasis, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_cfm_create(hsolve, matrix_struct=fm_struct_tmp, name="trl_h", set_zero=.TRUE.)
            CALL cp_cfm_create(esolve, matrix_struct=fm_struct_tmp, name="trl_e", set_zero=.TRUE.)
            CALL cp_fm_struct_release(fm_struct_tmp)
            CALL cp_cfm_gemm('C', 'N', nbasis, nbasis, nao, z_one, vmat, hvmat, z_zero, hsolve)
            CALL cp_cfm_heevd(hsolve, esolve, theta(1:nbasis))
            nres = MIN(nmobt, nbasis)
            ! the single-block regime never reads beyond the leading nres
            ! pairs, so it gathers that block alone; the block regime may
            ! retain more than nres pairs at the restart (its k_save is
            ! capped by the budget, not by the pending width) and gathers
            ! the full spectrum
            IF (single_block) THEN
               CALL cp_cfm_get_submatrix(esolve, zvec(1:nbasis, 1:nres), 1, 1, nbasis, nres)
            ELSE
               CALL cp_cfm_get_submatrix(esolve, zt(1:nbasis, 1:nbasis), 1, 1, nbasis, nbasis)
            END IF
            CALL cp_cfm_release(hsolve)
            CALL cp_cfm_release(esolve)

            ! ---------- exact residuals of the target pairs: the leading
            ! Ritz pairs, one per pending column
            IF (single_block) THEN
               CALL cp_cfm_set_submatrix(yw1, zvec(1:nbasis, 1:nres), 1, 1, nbasis, nres)
               CALL cp_cfm_set_submatrix(yw2, zvec(1:nbasis, 1:nres), 1, 1, nbasis, nres)
            ELSE
               CALL cp_cfm_set_submatrix(yw1, zt(1:nbasis, 1:nres), 1, 1, nbasis, nres)
               CALL cp_cfm_set_submatrix(yw2, zt(1:nbasis, 1:nres), 1, 1, nbasis, nres)
            END IF
            ! the scaling array spans every global column of yw2: the
            ! parallel branch of cp_cfm_column_scale indexes it with the
            ! global column id, a rank owning columns above nres would
            ! read past a shorter array; the dead tail scales by one
            ctheta(:) = z_one
            ctheta(1:nres) = CMPLX(theta(1:nres), 0.0_dp, KIND=dp)
            CALL cp_cfm_column_scale(yw2, ctheta)
            CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_one, hvmat, yw1, z_zero, rblk)
            CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_mone, vmat, yw2, z_one, rblk)
            CALL cp_cfm_vectorsnorm(rblk, nrm)
            rhot(1:nres) = nrm(1:nres)

            ! ---------- re-judged write-back: every cycle rotates the whole
            ! target block into place and refreshes its eigenvalues, so no
            ! pair-to-column identity ever latches across cycles (heevd
            ! re-sorts ascending, and an expansion that captures a lower
            ! state than the subspace held before shifts every index above
            ! it). This is variationally safe: the seed block contains the
            ! purified entry columns, their projected matrix equals
            ! diag(entry evals), so Courant-Fischer bounds every refreshed
            ! eigenvalue from above by its entry value. The staging stays
            ! off rblk: a restart below seeds from the residual vectors
            ! held there. Columns beyond nres (basis saturated below the
            ! pending width, rank-deficient seeds) keep their entry vectors
            ! and fall to the end-of-call fallback
            CALL cp_cfm_gemm('N', 'N', nao, nres, nbasis, z_one, vmat, yw1, z_zero, pblk)
            CALL cp_cfm_to_cfm(pblk, c0, nres, 1, col_lo)
            DO j = 1, nres
               eval(col_lo + j - 1) = theta(j)
               ldone(j) = rhot(j) < eps_col(col_lo + j - 1)
            END DO
            IF (ALL(ldone(1:nmobt))) EXIT cycle_loop

            ! ---------- thick restart: keep the lowest k_save Ritz pairs,
            ! buffered through c1 because the products read the whole basis
            IF (single_block) THEN
               k_save = MIN(nmobt, nbmax - 2, nbasis - 2)
            ELSE
               k_save = MIN(MAX(nbmax/2, 1), nbmax - 2, nbasis - 2)
            END IF
            ! a target pair never leaves the working basis at the restart:
            ! a direction the block still needs would have to be re-derived
            ! from the seeds by a later cycle
            k_save = MAX(k_save, MIN(nres, nbasis - 2))
            IF (k_save < 1) EXIT cycle_loop
            ! the eigenvector block is routed through tmatrix (it already
            ! lives in the broadcast buffer) and both products form as single gemms
            IF (single_block) THEN
               CALL cp_cfm_set_submatrix(tmatrix, zvec(1:nbasis, 1:k_save), 1, 1, nbasis, k_save)
            ELSE
               CALL cp_cfm_set_submatrix(tmatrix, zt(1:nbasis, 1:k_save), 1, 1, nbasis, k_save)
            END IF
            CALL cp_cfm_gemm('N', 'N', nao, k_save, nbasis, z_one, vmat, tmatrix, z_zero, restart_buf, &
                             c_first_col=1)
            CALL cp_cfm_gemm('N', 'N', nao, k_save, nbasis, z_one, hvmat, tmatrix, z_zero, restart_buf, &
                             c_first_col=k_save + 1)
            CALL cp_cfm_to_cfm(restart_buf, vmat, k_save, 1, 1)
            CALL cp_cfm_to_cfm(restart_buf, hvmat, k_save, k_save + 1, 1)
            nbasis = k_save

            ! ---------- new seed block: the residual vectors of the
            ! still-pending pairs, compacted to the front of rblk
            nstep = 0
            DO j = 1, nres
               IF (.NOT. ldone(j)) THEN
                  nstep = nstep + 1
                  gscale(nstep) = rhot(j)
                  IF (nstep /= j) THEN
                     CALL cp_cfm_to_cfm(rblk, wvec, 1, j, 1)
                     CALL cp_cfm_to_cfm(wvec, rblk, 1, 1, nstep)
                  END IF
               END IF
            END DO
            IF (nstep == 0) EXIT cycle_loop
            ! the seed batch is clamped to the budget like the expansion: in
            ! the single-block regime the restart keeps every pending pair,
            ! and those pairs plus their seeds can exceed the basis budget
            nmob_c = MIN(nstep, nbmax - nbasis)
            CALL purify_and_append(c0, col_lo - 1, rblk, nmob_c, gscale, vmat, nbasis, nseed, &
                                   yw1, proj_cb, tmatrix, wvec, nrm)
            IF (nseed == 0) EXIT cycle_loop
            nstep = nseed
         END DO cycle_loop

         ! exit statistics: the residual numbers describe the state the
         ! columns are left in (the final cycle's wanted pairs, not the
         ! worst value seen during the call, and not columns beyond the
         ! final nres, which keep their entry vectors). The converged
         ! count is the global leading converged run: once a block ends
         ! short of full convergence, no later column can extend the
         ! prefix, so later blocks contribute nothing even when all of
         ! their own columns converged
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

      DEALLOCATE (theta, nrm, gscale, rhot, ldone, ctheta, zt, zvec)
      CALL cp_cfm_release(vmat)
      CALL cp_cfm_release(hvmat)
      CALL cp_cfm_release(restart_buf)
      CALL cp_cfm_release(pblk)
      CALL cp_cfm_release(rblk)
      CALL cp_cfm_release(wvec)
      CALL cp_cfm_release(proj_cb)
      CALL cp_cfm_release(tmatrix)
      CALL cp_cfm_release(yw1)
      CALL cp_cfm_release(yw2)

      CALL timestop(handle)
   END SUBROUTINE trlan_refine

! **************************************************************************************************
!> \brief Purify a block of candidate directions and append the survivors to
!>        the working basis: two full projection passes against the working
!>        basis and the leading ndefl columns of c0, modified Gram-Schmidt
!>        within the block, then a per-column rank safeguard relative to
!>        gscale, normalization and append. The count of appended columns
!>        is returned in napp; nbasis is advanced accordingly.
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
      ! candidate set against the working basis and the deflation prefix,
      ! the same double-pass quality the single-column path provides
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
      ! dimensions the orthogonal complement is pure cancellation noise;
      ! a residual of 1e-16 passes a plain positivity test, as
      ! He-minimal-basis tests with NKRYLOV larger than the basis showed
      IF (ALL(nrm(1:nmob) <= 1.0E-12_dp*gscale(1:nmob))) THEN
         napp = 0
         RETURN
      END IF

      ! within-batch orthogonalization by CholeskyQR2 (two rounds of a
      ! Cholesky factorization of the batch Gram): near-parallel candidates
      ! pass the per-column norm guard yet square the condition of the
      ! first round into its roundoff, and the second round removes that
      ! error again. A batch of nearly-spanned candidates (the usual late-SCF
      ! expansion) carries a few exhausted columns: only those are retired,
      ! the healthy remainder is compacted to the front and keeps the
      ! batched treatment, so a single deficient direction cannot demote
      ! the whole batch to the column-at-a-time loop below. Only a
      ! factorization that reports rank deficiency on the healthy set does
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
         ! global column id, so a shorter array reads out of bounds on
         ! ranks owning columns above the live count; the retired tail
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
      ! previous candidates, which makes the within-batch orthogonalization
      ! a byproduct. A surviving first CholeskyQR2 round already
      ! renormalized the batch, so the rank test must reference the
      ! current norms, not the pre-purge scales the caller passed
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

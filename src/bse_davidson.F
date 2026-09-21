!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block Davidson solvers for the lowest excitations of the Bethe-Salpeter equation on top
!>        of the matrix-free application of A and B
!> \par History
!>      09.2026 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_davidson

   USE bse_full_diag,                   ONLY: create_hermitian_form_of_ABBA
   USE bse_matvec,                      ONLY: bse_matvec_apply,&
                                              bse_matvec_diagonal,&
                                              bse_matvec_env_type,&
                                              bse_matvec_subblock,&
                                              bse_matvec_vector_struct
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_transpose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_release, cp_fm_set_all, &
        cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type
   USE input_constants,                 ONLY: bse_iter_and_cond,&
                                              bse_iter_en_cond,&
                                              bse_iter_or_cond,&
                                              bse_iter_res_cond
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_davidson'

   ! energies closer than this belong to one multiplet (Hartree)
   REAL(KIND=dp), PARAMETER, PRIVATE :: deg_thresh = 1.0E-6_dp
   ! a basis that is not orthonormal to this accuracy makes the residuals meaningless
   REAL(KIND=dp), PARAMETER, PRIVATE :: ortho_tol = 1.0E-8_dp

   ! outcome of an ABBA solver: solved, or A-B not positive definite on the trial space
   INTEGER, PARAMETER, PUBLIC :: abba_ok = 0, abba_indefinite = 1

   PUBLIC :: bse_davidson_tda, bse_davidson_abba_mk, bse_davidson_abba_os, &
             bse_davidson_refcheck_tda, bse_davidson_refcheck_abba

CONTAINS

! **************************************************************************************************
!> \brief Block Davidson with thick restart for the NUM_EXC_EN lowest solutions of
!>        sum_jb A_ia,jb X_jb^n = Ω^n X_ia^n. With the orthonormal basis Z_ia,m the reduced matrix is
!>        G_mn = sum_ia Z_ia,m (A Z)_ia,n with eigenpairs (θ_k, c_mk), X_ia^k = sum_m Z_ia,m c_mk,
!>        and the basis grows by the corrections t_ia,k = r_ia,k/(d_ia - θ_k) of the residuals
!>        r_ia,k = sum_m (A Z)_ia,m c_mk - θ_k X_ia^k, d_ia being the diagonal chosen by PRECONDITIONER
!> \param mv_env ...
!> \param mp2_env ...
!> \param unit_nr ...
!> \param exc_ens Ω^n in Hartree
!> \param fm_X X_ia^n on the process grid of mv_env, column n
! **************************************************************************************************
   SUBROUTINE bse_davidson_tda(mv_env, mp2_env, unit_nr, exc_ens, fm_X)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: exc_ens
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_X

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_davidson_tda'

      INTEGER                                            :: block_size, handle, iter, k, m, m_max, &
                                                            n_act, n_cand, n_dependent, n_dropped, &
                                                            n_kernel, n_ov, n_req, n_restart, &
                                                            n_want, nt
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected
      LOGICAL                                            :: all_conv, has_prev, ok, stalled
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: conv, need
      REAL(KIND=dp)                                      :: dev, t_iter
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, diag, res, theta, theta_prev, &
                                                            theta_sub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C, C_cand, C_ritz, G_new, Q
      TYPE(cp_blacs_env_type), POINTER                   :: diag_blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_AZ, fm_red_A, fm_work, fm_Z
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      para_env => mv_env%para_env
      ! the trial vectors live on an npe x 1 grid, which is a bad one to diagonalize on
      NULLIFY (diag_blacs_env)
      CALL cp_blacs_env_create(diag_blacs_env, para_env)
      n_ov = mv_env%n_ov
      CALL davidson_sizes(mp2_env, n_ov, n_want, n_act, block_size)

      ALLOCATE (diag(n_ov))
      CALL bse_matvec_diagonal(mv_env, mp2_env%bse%preconditioner, diag)

      CALL subspace_ceiling(mp2_env, n_ov, n_act, 2, para_env, unit_nr, m_max)

      CALL bse_matvec_vector_struct(mv_env, m_max, fm_struct)
      CALL cp_fm_create(fm_Z, fm_struct, name="fm_Z_davidson", set_zero=.TRUE.)
      CALL cp_fm_create(fm_AZ, fm_struct, name="fm_AZ_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL bse_matvec_vector_struct(mv_env, MIN(2*n_act, m_max), fm_struct)
      CALL cp_fm_create(fm_work, fm_struct, name="fm_work_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=diag_blacs_env, &
                               nrow_global=m_max, ncol_global=m_max)
      CALL cp_fm_create(fm_red_A, fm_struct, name="fm_red_A_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (C(m_max, n_act), C_cand(m_max, 2*n_act), C_ritz(m_max, n_act))
      ALLOCATE (Q(m_max, 2*n_act))
      ALLOCATE (theta(m_max), theta_prev(n_act), res(n_act), dE(n_act), conv(n_act), need(n_act), &
                selected(n_act))
      C_cand(:, :) = 0.0_dp
      theta_prev(:) = 0.0_dp
      has_prev = .FALSE.
      stalled = .FALSE.
      n_restart = 0
      n_dependent = 0

      IF (mp2_env%bse%num_guess_transitions > 0) THEN
         CALL initial_guess_subblock(mv_env, mp2_env, diag, n_act, m_max, deg_thresh, fm_Z, m, theta_sub, &
                                     unit_nr)
      ELSE
         CALL initial_guess(diag, n_act, m_max, deg_thresh, fm_Z, m)
      END IF
      IF (extension_deviation(fm_Z, fm_Z, 0, m, para_env) > ortho_tol) THEN
         CPABORT("BSE Davidson: the initial guess is not orthonormal")
      END IF
      ! (A Z)_ia,m of the guess and G_mn = sum_ia Z_ia,m (A Z)_ia,n
      CALL bse_matvec_apply(mv_env, fm_Z, 1, m, fm_AZ)
      n_kernel = m
      CALL reduced_gram_blocks(fm_Z, m, fm_AZ, 1, m, block_size, para_env, fm_red_A)

      CALL print_iteration_header('Block Davidson iterations within the TDA:', unit_nr)

      all_conv = .FALSE.
      DO iter = 1, mp2_env%bse%max_iter
         t_iter = m_walltime()

         ! sum_n G_mn c_nk = θ_k c_mk, θ ascending
         CALL solve_reduced(fm_red_A, m, n_act, para_env, diag_blacs_env, theta, C)

         ! thick restart onto the current and the previous Ritz vectors, no kernel application
         IF (m + block_size > m_max .AND. m > n_act .AND. m < n_ov) THEN
            C_cand(1:m, 1:n_act) = C(1:m, 1:n_act)
            n_cand = n_act
            IF (has_prev) n_cand = 2*n_act
            CALL restart_basis(C_cand, m, n_cand, 2*n_act, Q, k)
            CALL rotate_in_place(fm_Z, m, Q, k, fm_work)
            CALL rotate_in_place(fm_AZ, m, Q, k, fm_work)
            CALL reduced_rotate(fm_red_A, m, Q, k, para_env, diag_blacs_env)
            m = k
            n_restart = n_restart + 1
            CALL solve_reduced(fm_red_A, m, n_act, para_env, diag_blacs_env, theta, C)
         END IF

         ! r_ia,k = sum_m (A Z)_ia,m c_mk - θ_k sum_m Z_ia,m c_mk for all tracked roots k
         DO k = 1, n_act
            C_ritz(1:m, k) = C(1:m, k)*theta(k)
         END DO
         CALL subspace_rotate(fm_AZ, m, C(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
         CALL subspace_rotate(fm_Z, m, C_ritz(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
         CALL column_norms(fm_work, 1, n_act, para_env, res)

         dE(:) = ABS(theta(1:n_act) - theta_prev(:))
         CALL convergence_flags(mp2_env, dE, res, conv)
         n_req = last_required_root(theta, n_want, n_act)
         ! a subspace that spans all transitions is exact
         IF (m == n_ov) conv(:) = .TRUE.
         CALL roots_in_need(conv, theta, res, n_req, need)
         all_conv = .NOT. ANY(need)
         CALL print_iteration(iter, m, COUNT(conv(1:n_want)), &
                              MAX(MAXVAL(res(1:n_req)), MAXVAL(res, MASK=need)), &
                              MAXVAL(dE(1:n_req)), theta(1), m_walltime() - t_iter, unit_nr)
         IF (all_conv) EXIT

         IF (iter == mp2_env%bse%max_iter) EXIT

         ! corrections for the lowest unconverged roots, t_ia,k = r_ia,k/(d_ia - θ_k)
         CALL select_roots(mp2_env, need, res, MIN(block_size, m_max - m), selected, nt)
         CALL davidson_corrections(fm_work, selected(1:nt), theta, diag, fm_Z, m + 1)
         CALL project_out(fm_Z, m + 1, nt, fm_Z, m, para_env)
         CALL drop_small_columns(fm_Z, m + 1, nt, para_env)
         ! no correction left: the open roots only wait for their energy criterion, allowed once
         IF (nt == 0 .AND. stalled) THEN
            CPABORT("BSE Davidson: correction vectors vanished before convergence")
         END IF
         stalled = nt == 0
         CALL cholqr2(fm_Z, m + 1, nt, para_env, ok, n_dropped)
         IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")
         n_dependent = n_dependent + n_dropped
         IF (extension_deviation(fm_Z, fm_Z, m, nt, para_env) > ortho_tol) THEN
            CPABORT("BSE Davidson: the basis lost orthonormality")
         END IF

         IF (nt > 0) THEN
            CALL bse_matvec_apply(mv_env, fm_Z, m + 1, nt, fm_AZ)
            n_kernel = n_kernel + nt

            ! new columns G_mn = sum_ia Z_ia,m (A Z)_ia,n of the symmetric reduced matrix
            ALLOCATE (G_new(m + nt, nt))
            CALL subspace_gram(fm_Z, 1, m + nt, fm_AZ, m + 1, nt, para_env, G_new)
            CALL reduced_extend(fm_red_A, G_new, m, nt)
            DEALLOCATE (G_new)
         END IF

         ! the Ritz vectors of this iteration are the second block of restart candidates
         C_cand(:, n_act + 1:2*n_act) = 0.0_dp
         C_cand(1:m, n_act + 1:2*n_act) = C(1:m, 1:n_act)
         has_prev = .TRUE.
         theta_prev(:) = theta(1:n_act)
         m = m + nt
      END DO

      IF (.NOT. all_conv) CALL abort_unconverged(theta, res, conv, need, n_req, unit_nr)
      IF (mp2_env%bse%bse_debug_print) CALL print_tracked_roots(theta, res, n_want, n_act, unit_nr)

      IF (mp2_env%bse%bse_debug_print) THEN
         dev = orthonormality_deviation(fm_Z, fm_Z, m, block_size, para_env)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A10,T13,A,T63,ES18.6)') 'BSE|DEBUG|', &
               'Max deviation of the basis from orthonormality', dev
         END IF
      END IF

      IF (ALLOCATED(theta_sub)) CALL check_guess_bound(theta_sub, theta, n_req, unit_nr)
      CALL print_summary(iter, n_kernel, n_restart, n_want, n_req, n_act, n_ov, unit_nr, n_dependent=n_dependent)

      ALLOCATE (exc_ens(n_want))
      exc_ens(:) = theta(1:n_want)
      CALL bse_matvec_vector_struct(mv_env, n_want, fm_struct)
      CALL cp_fm_create(fm_X, fm_struct, name="fm_X_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      ! X_ia^n = sum_m Z_ia,m c_mn
      CALL subspace_rotate(fm_Z, m, C(1:m, 1:n_want), n_want, fm_X, 1.0_dp, 0.0_dp)

      CALL cp_fm_release(fm_Z)
      CALL cp_fm_release(fm_AZ)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_red_A)
      CALL cp_blacs_env_release(diag_blacs_env)
      DEALLOCATE (C, C_cand, C_ritz, Q, theta, theta_prev, res, dE, conv, need, selected, diag)
      IF (ALLOCATED(theta_sub)) DEALLOCATE (theta_sub)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_tda

! **************************************************************************************************
!> \brief Block Davidson with thick restart for the NUM_EXC_EN lowest solutions of the full problem,
!>        written as sum_jb [(A+B)(A-B)]_ia,jb x_jb^n = (Ω^n)^2 x_ia^n with x^n = X^n - Y^n, which is
!>        symmetric in the inner product <u,v> = sum_ia,jb u_ia (A-B)_ia,jb v_jb. The basis V_ia,m is
!>        orthonormal in that inner product, W_ia,m = sum_jb (A-B)_ia,jb V_jb,m, and the reduced matrix
!>        H_mn = sum_ia,jb W_ia,m (A+B)_ia,jb W_jb,n has the eigenpairs ((Ω^k)^2, c_mk). Then
!>        x_ia^k = (Ω^k)^1/2 sum_m V_ia,m c_mk and y_ia^k = (X^k+Y^k)_ia = (Ω^k)^-1/2 sum_m W_ia,m c_mk
!>        obey sum_ia x_ia^k y_ia^k = 1, the residual is r_ia,k = sum_jb (A+B)_ia,jb y_jb^k - Ω^k x_ia^k
!>        and the correction t_ia,k = r_ia,k/(d_ia^2 - (Ω^k)^2), d_ia chosen by PRECONDITIONER
!> \param mv_env ...
!> \param mp2_env ...
!> \param unit_nr ...
!> \param exc_ens Ω^n in Hartree
!> \param fm_X X_ia^n on the process grid of mv_env, column n; created only for abba_ok
!> \param fm_Y Y_ia^n, as fm_X
!> \param abba_status abba_indefinite if A-B is not positive definite on the trial space
!> \param ab_margin smallest eigenvalue of A-B on the trial space in Hartree
!> \param fm_X_tda converged TDA vectors as initial guess
! **************************************************************************************************
   SUBROUTINE bse_davidson_abba_mk(mv_env, mp2_env, unit_nr, exc_ens, fm_X, fm_Y, abba_status, &
                                   ab_margin, fm_X_tda)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: exc_ens
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_X, fm_Y
      INTEGER, INTENT(OUT)                               :: abba_status
      REAL(KIND=dp), INTENT(OUT)                         :: ab_margin
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_X_tda

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_davidson_abba_mk'

      INTEGER :: block_size, guess_col, handle, iter, k, m, m_max, n_act, n_cand, n_dependent, &
         n_guess_cols, n_kernel, n_ov, n_req, n_restart, n_seed, n_want, nrow_local, nt, nt_guess
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected
      LOGICAL                                            :: all_conv, has_prev, stalled
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: conv, need
      REAL(KIND=dp)                                      :: dev, t_iter
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, diag, omega, omega_prev, omega_sq, &
                                                            res, theta, theta_sub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C, C_cand, C_x, C_x_omega, C_y, H_new, Q
      TYPE(cp_blacs_env_type), POINTER                   :: diag_blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_MW, fm_red_M, fm_S, fm_scratch, &
                                                            fm_Svec, fm_V, fm_W, fm_work
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      abba_status = abba_ok
      ab_margin = 0.0_dp
      para_env => mv_env%para_env
      ! the trial vectors live on an npe x 1 grid, which is a bad one to diagonalize on
      NULLIFY (diag_blacs_env)
      CALL cp_blacs_env_create(diag_blacs_env, para_env)
      n_ov = mv_env%n_ov
      CALL davidson_sizes(mp2_env, n_ov, n_want, n_act, block_size)

      ALLOCATE (diag(n_ov))
      CALL bse_matvec_diagonal(mv_env, mp2_env%bse%preconditioner, diag)

      CALL subspace_ceiling(mp2_env, n_ov, n_act, 3, para_env, unit_nr, m_max)

      CALL bse_matvec_vector_struct(mv_env, m_max, fm_struct)
      CALL cp_fm_create(fm_V, fm_struct, name="fm_V_davidson", set_zero=.TRUE.)
      CALL cp_fm_create(fm_W, fm_struct, name="fm_W_davidson", set_zero=.TRUE.)
      CALL cp_fm_create(fm_MW, fm_struct, name="fm_MW_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL bse_matvec_vector_struct(mv_env, MIN(2*n_act, m_max), fm_struct)
      CALL cp_fm_create(fm_work, fm_struct, name="fm_work_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL bse_matvec_vector_struct(mv_env, block_size, fm_struct)
      CALL cp_fm_create(fm_scratch, fm_struct, name="fm_scratch_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=diag_blacs_env, &
                               nrow_global=m_max, ncol_global=m_max)
      CALL cp_fm_create(fm_red_M, fm_struct, name="fm_red_M_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local)

      ALLOCATE (C(m_max, n_act), C_cand(m_max, 2*n_act), Q(m_max, 2*n_act))
      ALLOCATE (C_x(m_max, n_act), C_y(m_max, n_act), C_x_omega(m_max, n_act))
      ALLOCATE (theta(m_max), omega(n_act), omega_prev(n_act), res(n_act), dE(n_act), conv(n_act), &
                need(n_act), selected(n_act))
      C_cand(:, :) = 0.0_dp
      omega_prev(:) = 0.0_dp
      has_prev = .FALSE.
      stalled = .FALSE.
      all_conv = .FALSE.
      n_restart = 0
      n_kernel = 0
      n_dependent = 0

      ! guess: the TDA vectors if given, then the eigenvectors of the exact guess block or unit
      ! vectors on the lowest entries of the diagonal
      n_seed = 0
      IF (PRESENT(fm_X_tda)) THEN
         CALL cp_fm_get_info(fm_X_tda, ncol_global=n_seed)
         n_seed = MIN(n_seed, n_act)
         CALL cp_fm_to_fm(fm_X_tda, fm_V, n_seed)
      END IF
      IF (mp2_env%bse%num_guess_transitions > 0) THEN
         CALL initial_guess_subblock(mv_env, mp2_env, diag, n_act, m_max, deg_thresh, fm_V, n_guess_cols, theta_sub, &
                                     unit_nr, n_seed, abba_status, ab_margin)
      ELSE
         CALL initial_guess(diag, n_act, m_max, deg_thresh, fm_V, n_guess_cols, n_seed)
      END IF
      ! the corrections divide by d_ia^2 - (Ω^k)^2
      diag(:) = diag(:)**2

      ! the guess enters the basis in blocks of the scratch width
      m = 0
      guess_col = 1
      DO WHILE (guess_col <= n_guess_cols)
         nt_guess = MIN(block_size, n_guess_cols - guess_col + 1)
         nt = nt_guess
         IF (guess_col > m + 1) THEN
            CALL cp_fm_to_fm(fm_V, fm_V, nt, guess_col, m + 1)
         END IF
         CALL extend_basis_mk(mv_env, fm_V, fm_W, fm_MW, fm_scratch, m, nt, n_kernel, n_dependent, abba_status, &
                              ab_margin)
         IF (abba_status /= abba_ok) EXIT
         guess_col = guess_col + nt_guess
         m = m + nt
      END DO

      IF (abba_status == abba_ok) THEN
         fm_V%local_data(1:nrow_local, m + 1:n_guess_cols) = 0.0_dp
         ! H_mn = sum_ia W_ia,m (M W)_ia,n on the guess basis
         CALL reduced_gram_blocks(fm_W, m, fm_MW, 1, m, block_size, para_env, fm_red_M)

         CALL print_iteration_header('Block Davidson iterations for ABBA, (A+B)(A-B) x = E^2 x:', unit_nr)

         DO iter = 1, mp2_env%bse%max_iter
            t_iter = m_walltime()

            ! sum_n H_mn c_nk = θ_k c_mk with θ_k = (Ω^k)^2, ascending
            CALL solve_reduced(fm_red_M, m, n_act, para_env, diag_blacs_env, theta, C)

            ! thick restart: V, W and MW are rotated alike, which keeps V^T W = 1
            IF (m + block_size > m_max .AND. m > n_act .AND. m < n_ov) THEN
               C_cand(1:m, 1:n_act) = C(1:m, 1:n_act)
               n_cand = n_act
               IF (has_prev) n_cand = 2*n_act
               CALL restart_basis(C_cand, m, n_cand, 2*n_act, Q, k)
               CALL rotate_in_place(fm_V, m, Q, k, fm_work)
               CALL rotate_in_place(fm_W, m, Q, k, fm_work)
               CALL rotate_in_place(fm_MW, m, Q, k, fm_work)
               CALL reduced_rotate(fm_red_M, m, Q, k, para_env, diag_blacs_env)
               m = k
               n_restart = n_restart + 1
               CALL solve_reduced(fm_red_M, m, n_act, para_env, diag_blacs_env, theta, C)
            END IF

            ! a Ritz value bounds an eigenvalue of (A+B)(A-B) from above
            IF (theta(1) <= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson (MK_DAVIDSON): negative squared excitation energy. Matrix "// &
                             "(A+B) is not positive definite, or the basis lost its orthonormality.")
            END IF
            omega(:) = SQRT(theta(1:n_act))

            ! r_ia,k = sum_m (M W)_ia,m c_mk (Ω^k)^-1/2 - Ω^k (Ω^k)^1/2 sum_m V_ia,m c_mk
            DO k = 1, n_act
               C_y(1:m, k) = C(1:m, k)/SQRT(omega(k))
               C_x_omega(1:m, k) = C(1:m, k)*omega(k)*SQRT(omega(k))
            END DO
            CALL subspace_rotate(fm_MW, m, C_y(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
            CALL subspace_rotate(fm_V, m, C_x_omega(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
            CALL column_norms(fm_work, 1, n_act, para_env, res)

            dE(:) = ABS(omega(:) - omega_prev(:))
            CALL convergence_flags(mp2_env, dE, res, conv)
            n_req = last_required_root(omega, n_want, n_act)
            IF (m == n_ov) conv(:) = .TRUE.
            CALL roots_in_need(conv, omega, res, n_req, need)
            all_conv = .NOT. ANY(need)
            CALL print_iteration(iter, m, COUNT(conv(1:n_want)), &
                                 MAX(MAXVAL(res(1:n_req)), MAXVAL(res, MASK=need)), &
                                 MAXVAL(dE(1:n_req)), omega(1), m_walltime() - t_iter, unit_nr)
            IF (all_conv) EXIT

            IF (iter == mp2_env%bse%max_iter) EXIT

            ! corrections t_ia,k = r_ia,k/(d_ia^2 - (Ω^k)^2) for the lowest unconverged roots
            CALL select_roots(mp2_env, need, res, MIN(block_size, m_max - m), selected, nt)
            CALL davidson_corrections(fm_work, selected(1:nt), theta, diag, fm_V, m + 1)
            CALL extend_basis_mk(mv_env, fm_V, fm_W, fm_MW, fm_scratch, m, nt, n_kernel, n_dependent, abba_status, &
                                 ab_margin)
            IF (abba_status /= abba_ok) EXIT
            IF (nt == 0 .AND. stalled) THEN
               CPABORT("BSE Davidson: correction vectors vanished before convergence")
            END IF
            stalled = nt == 0

            IF (nt > 0) THEN
               ! new columns H_mn = sum_ia W_ia,m (M W)_ia,n of the symmetric reduced matrix
               ALLOCATE (H_new(m + nt, nt))
               CALL subspace_gram(fm_W, 1, m + nt, fm_MW, m + 1, nt, para_env, H_new)
               CALL reduced_extend(fm_red_M, H_new, m, nt)
               DEALLOCATE (H_new)
            END IF

            ! the Ritz vectors of this iteration are the second block of restart candidates
            C_cand(:, n_act + 1:2*n_act) = 0.0_dp
            C_cand(1:m, n_act + 1:2*n_act) = C(1:m, 1:n_act)
            has_prev = .TRUE.
            omega_prev(:) = omega(:)
            m = m + nt
         END DO
      END IF

      IF (abba_status == abba_ok) THEN
         IF (.NOT. all_conv) CALL abort_unconverged(omega, res, conv, need, n_req, unit_nr)
         IF (mp2_env%bse%bse_debug_print) CALL print_tracked_roots(omega, res, n_want, n_act, unit_nr)

         IF (mp2_env%bse%bse_debug_print) THEN
            dev = orthonormality_deviation(fm_V, fm_W, m, block_size, para_env)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A10,T13,A,T63,ES18.6)') 'BSE|DEBUG|', &
                  'Max deviation of the basis from orthonormality', dev
            END IF
         END IF

         ! with sum_ia V_ia,m W_ia,n = δ_mn the smallest eigenvalue of A-B on the span of V is the
         ! inverse of the largest eigenvalue of S_mn = sum_ia V_ia,m V_ia,n
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=diag_blacs_env, &
                                  nrow_global=m, ncol_global=m)
         CALL cp_fm_create(fm_S, fm_struct, name="bse_basis_metric")
         CALL cp_fm_create(fm_Svec, fm_struct, name="bse_basis_metric_vectors")
         CALL cp_fm_struct_release(fm_struct)
         CALL reduced_gram_blocks(fm_V, m, fm_V, 1, m, block_size, para_env, fm_S)
         ALLOCATE (omega_sq(m))
         omega_sq(:) = 0.0_dp
         CALL choose_eigv_solver(fm_S, fm_Svec, omega_sq)
         ab_margin = 1.0_dp/omega_sq(m)
         DEALLOCATE (omega_sq)
         CALL cp_fm_release(fm_S)
         CALL cp_fm_release(fm_Svec)

         IF (ALLOCATED(theta_sub)) CALL check_guess_bound(theta_sub, omega, n_req, unit_nr)
         CALL print_summary(iter, n_kernel, n_restart, n_want, n_req, n_act, n_ov, unit_nr, ab_margin, n_dependent)

         ! x_ia^n = (Ω^n)^1/2 sum_m V_ia,m c_mn, y_ia^n = (Ω^n)^-1/2 sum_m W_ia,m c_mn,
         ! X = (y + x)/2, Y = (y - x)/2
         ALLOCATE (exc_ens(n_want))
         exc_ens(:) = omega(1:n_want)
         DO k = 1, n_want
            C_x(1:m, k) = C(1:m, k)*SQRT(omega(k))
            C_y(1:m, k) = C(1:m, k)/SQRT(omega(k))
         END DO
         CALL bse_matvec_vector_struct(mv_env, n_want, fm_struct)
         CALL cp_fm_create(fm_X, fm_struct, name="fm_X_davidson", set_zero=.TRUE.)
         CALL cp_fm_create(fm_Y, fm_struct, name="fm_Y_davidson", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct)
         CALL subspace_rotate(fm_W, m, C_y(1:m, 1:n_want), n_want, fm_X, 0.5_dp, 0.0_dp)
         CALL subspace_rotate(fm_V, m, C_x(1:m, 1:n_want), n_want, fm_X, 0.5_dp, 1.0_dp)
         CALL subspace_rotate(fm_W, m, C_y(1:m, 1:n_want), n_want, fm_Y, 0.5_dp, 0.0_dp)
         CALL subspace_rotate(fm_V, m, C_x(1:m, 1:n_want), n_want, fm_Y, -0.5_dp, 1.0_dp)
      END IF

      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_W)
      CALL cp_fm_release(fm_MW)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_scratch)
      CALL cp_fm_release(fm_red_M)
      CALL cp_blacs_env_release(diag_blacs_env)
      DEALLOCATE (C, C_cand, C_x, C_y, C_x_omega, Q, theta, omega, omega_prev, res, dE, conv, need, &
                  selected, diag)
      IF (ALLOCATED(theta_sub)) DEALLOCATE (theta_sub)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_abba_mk

! **************************************************************************************************
!> \brief Block Davidson with thick restart for the NUM_EXC_EN lowest solutions of the full problem
!>        in the paired form of Olsen and Stratmann: with x^n = X^n - Y^n, y^n = X^n + Y^n, M = A+B and
!>        K = A-B, sum_jb M_ia,jb y_jb^n = Ω^n x_ia^n and sum_jb K_ia,jb x_jb^n = Ω^n y_ia^n. One
!>        orthonormal basis b_ia,m carries both vectors: Mr_mn = sum_ia b_ia,m (M b)_ia,n and
!>        Kr_mn = sum_ia b_ia,m (K b)_ia,n are solved by solve_reduced_paired for Ω^k, R_mk and L_mk with
!>        y_ia^k = sum_m b_ia,m R_mk, x_ia^k = sum_m b_ia,m L_mk and sum_ia x_ia^k y_ia^k = 1. The basis
!>        grows by the corrections t_ia,k = r_ia,k/(d_ia - Ω^k) of both residuals,
!>        rR_ia,k = sum_m (M b)_ia,m R_mk - Ω^k x_ia^k and rL_ia,k = sum_m (K b)_ia,m L_mk - Ω^k y_ia^k,
!>        d_ia chosen by PRECONDITIONER; the rR block enters before the rL block, which is nearly
!>        parallel to it for a small B
!> \param mv_env ...
!> \param mp2_env ...
!> \param unit_nr ...
!> \param exc_ens Ω^n in Hartree
!> \param fm_X X_ia^n on the process grid of mv_env, column n; created only for abba_ok
!> \param fm_Y Y_ia^n, as fm_X
!> \param abba_status abba_indefinite if A-B is not positive definite on the trial space
!> \param ab_margin smallest eigenvalue of A-B on the trial space in Hartree
!> \param fm_X_tda converged TDA vectors as initial guess, x = y = X
! **************************************************************************************************
   SUBROUTINE bse_davidson_abba_os(mv_env, mp2_env, unit_nr, exc_ens, fm_X, fm_Y, abba_status, &
                                   ab_margin, fm_X_tda)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: exc_ens
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_X, fm_Y
      INTEGER, INTENT(OUT)                               :: abba_status
      REAL(KIND=dp), INTENT(OUT)                         :: ab_margin
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_X_tda

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_davidson_abba_os'

      INTEGER :: block_size, guess_col, handle, iter, k, m, m_max, n_act, n_cand, n_dependent, &
         n_guess_cols, n_kernel, n_ov, n_pair, n_req, n_restart, n_seed, n_want, nrow_local, nt, &
         nt_guess, nt_L, nt_R, nt_sel
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected
      LOGICAL                                            :: all_conv, has_prev, stalled
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: conv, need
      REAL(KIND=dp)                                      :: dev, kr_min, t_iter
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, diag, omega, omega_prev, res, res_L, &
                                                            theta_sub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C_cand, C_work, Kr_new, L, Mr_new, Q, R
      TYPE(cp_blacs_env_type), POINTER                   :: diag_blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_b, fm_Kb, fm_Mb, fm_red_K, fm_red_M, &
                                                            fm_scratch, fm_work
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      abba_status = abba_ok
      ab_margin = 0.0_dp
      para_env => mv_env%para_env
      ! the trial vectors live on an npe x 1 grid, which is a bad one to diagonalize on
      NULLIFY (diag_blacs_env)
      CALL cp_blacs_env_create(diag_blacs_env, para_env)
      n_ov = mv_env%n_ov
      CALL davidson_sizes(mp2_env, n_ov, n_want, n_act, block_size)

      ALLOCATE (diag(n_ov))
      CALL bse_matvec_diagonal(mv_env, mp2_env%bse%preconditioner, diag)

      ! the current pairs and one block of corrections have to fit below the ceiling
      CALL subspace_ceiling(mp2_env, n_ov, n_act, 3, para_env, unit_nr, m_max, min_factor=4, work_factor=4)

      CALL bse_matvec_vector_struct(mv_env, m_max, fm_struct)
      CALL cp_fm_create(fm_b, fm_struct, name="fm_b_davidson", set_zero=.TRUE.)
      CALL cp_fm_create(fm_Mb, fm_struct, name="fm_Mb_davidson", set_zero=.TRUE.)
      CALL cp_fm_create(fm_Kb, fm_struct, name="fm_Kb_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      ! both residual blocks, or the restart basis
      CALL bse_matvec_vector_struct(mv_env, MAX(2*n_act, MIN(4*n_act, m_max)), fm_struct)
      CALL cp_fm_create(fm_work, fm_struct, name="fm_work_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL bse_matvec_vector_struct(mv_env, 2*block_size, fm_struct)
      CALL cp_fm_create(fm_scratch, fm_struct, name="fm_scratch_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=diag_blacs_env, &
                               nrow_global=m_max, ncol_global=m_max)
      CALL cp_fm_create(fm_red_M, fm_struct, name="fm_red_M_davidson", set_zero=.TRUE.)
      CALL cp_fm_create(fm_red_K, fm_struct, name="fm_red_K_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_get_info(fm_b, nrow_local=nrow_local)

      ALLOCATE (R(m_max, n_act), L(m_max, n_act), C_work(m_max, n_act))
      ALLOCATE (C_cand(m_max, 4*n_act), Q(m_max, 4*n_act))
      ALLOCATE (omega(n_act), omega_prev(n_act), res(n_act), res_L(n_act), dE(n_act), conv(n_act), &
                need(n_act), selected(n_act))
      C_cand(:, :) = 0.0_dp
      omega_prev(:) = 0.0_dp
      kr_min = 0.0_dp
      has_prev = .FALSE.
      stalled = .FALSE.
      all_conv = .FALSE.
      n_restart = 0
      n_kernel = 0
      n_dependent = 0

      ! guess: the TDA vectors if given (each a pair by itself), then the pairs x_k, y_k of the exact
      ! guess block or unit vectors on the lowest entries of the diagonal
      n_seed = 0
      IF (PRESENT(fm_X_tda)) THEN
         CALL cp_fm_get_info(fm_X_tda, ncol_global=n_seed)
         n_seed = MIN(n_seed, n_act)
         CALL cp_fm_to_fm(fm_X_tda, fm_b, n_seed)
      END IF
      n_pair = 0
      IF (mp2_env%bse%num_guess_transitions > 0) THEN
         CALL initial_guess_subblock(mv_env, mp2_env, diag, n_act, m_max, deg_thresh, fm_b, n_guess_cols, theta_sub, &
                                     unit_nr, n_seed, abba_status, ab_margin, n_pair)
      ELSE
         CALL initial_guess(diag, n_act, m_max, deg_thresh, fm_b, n_guess_cols, n_seed)
      END IF
      n_guess_cols = n_guess_cols + n_pair

      ! the guess enters the basis in blocks of the scratch width
      m = 0
      guess_col = 1
      DO WHILE (guess_col <= n_guess_cols)
         nt_guess = MIN(2*block_size, n_guess_cols - guess_col + 1)
         nt = nt_guess
         IF (guess_col > m + 1) THEN
            CALL cp_fm_to_fm(fm_b, fm_b, nt, guess_col, m + 1)
         END IF
         CALL extend_basis_os(mv_env, fm_b, fm_Mb, fm_Kb, fm_scratch, m, nt, n_kernel, n_dependent)
         guess_col = guess_col + nt_guess
         m = m + nt
      END DO

      IF (abba_status == abba_ok) THEN
         fm_b%local_data(1:nrow_local, m + 1:n_guess_cols) = 0.0_dp
         ! Mr_mn = sum_ia b_ia,m (M b)_ia,n and Kr_mn = sum_ia b_ia,m (K b)_ia,n on the guess basis
         CALL reduced_gram_blocks(fm_b, m, fm_Mb, 1, m, block_size, para_env, fm_red_M)
         CALL reduced_gram_blocks(fm_b, m, fm_Kb, 1, m, block_size, para_env, fm_red_K)

         CALL print_iteration_header('Block Davidson iterations for ABBA, Olsen-Stratmann paired subspace:', &
                                     unit_nr)

         DO iter = 1, mp2_env%bse%max_iter
            t_iter = m_walltime()

            ! Ω^k ascending with the coefficients R_mk of y^k and L_mk of x^k, sum_m R_mk L_mk = 1
            CALL solve_reduced_paired(fm_red_M, fm_red_K, m, n_act, para_env, omega, R, L, kr_min, &
                                      abba_status, diag_blacs_env)
            IF (abba_status /= abba_ok) EXIT

            ! thick restart: b, Mb and Kb are rotated alike; the current pairs come first, so that
            ! they survive the cap that leaves room for one block of corrections
            IF (m + 2*block_size > m_max .AND. m > 2*n_act .AND. m < n_ov) THEN
               C_cand(1:m, 1:n_act) = R(1:m, :)
               C_cand(1:m, n_act + 1:2*n_act) = L(1:m, :)
               n_cand = 2*n_act
               IF (has_prev) n_cand = 4*n_act
               CALL restart_basis(C_cand, m, n_cand, m_max - 2*block_size, Q, k)
               CALL rotate_in_place(fm_b, m, Q, k, fm_work)
               CALL rotate_in_place(fm_Mb, m, Q, k, fm_work)
               CALL rotate_in_place(fm_Kb, m, Q, k, fm_work)
               CALL reduced_rotate(fm_red_M, m, Q, k, para_env, diag_blacs_env)
               CALL reduced_rotate(fm_red_K, m, Q, k, para_env, diag_blacs_env)
               m = k
               n_restart = n_restart + 1
               CALL solve_reduced_paired(fm_red_M, fm_red_K, m, n_act, para_env, omega, R, L, kr_min, &
                                         abba_status, diag_blacs_env)
               IF (abba_status /= abba_ok) EXIT
            END IF

            ! rR_ia,k = sum_m (M b)_ia,m R_mk - Ω^k sum_m b_ia,m L_mk into the columns 1..n_act and
            ! rL_ia,k = sum_m (K b)_ia,m L_mk - Ω^k sum_m b_ia,m R_mk into the columns n_act+1..2 n_act
            DO k = 1, n_act
               C_work(1:m, k) = L(1:m, k)*omega(k)
            END DO
            CALL subspace_rotate(fm_Mb, m, R(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
            CALL subspace_rotate(fm_b, m, C_work(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
            DO k = 1, n_act
               C_work(1:m, k) = R(1:m, k)*omega(k)
            END DO
            CALL subspace_rotate(fm_Kb, m, L(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp, n_act + 1)
            CALL subspace_rotate(fm_b, m, C_work(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp, n_act + 1)
            CALL column_norms(fm_work, 1, n_act, para_env, res)
            CALL column_norms(fm_work, n_act + 1, n_act, para_env, res_L)
            res(:) = MAX(res(:), res_L(:))

            dE(:) = ABS(omega(:) - omega_prev(:))
            CALL convergence_flags(mp2_env, dE, res, conv)
            n_req = last_required_root(omega, n_want, n_act)
            IF (m == n_ov) conv(:) = .TRUE.
            CALL roots_in_need(conv, omega, res, n_req, need)
            all_conv = .NOT. ANY(need)
            CALL print_iteration(iter, m, COUNT(conv(1:n_want)), &
                                 MAX(MAXVAL(res(1:n_req)), MAXVAL(res, MASK=need)), &
                                 MAXVAL(dE(1:n_req)), omega(1), m_walltime() - t_iter, unit_nr)
            IF (all_conv) EXIT

            IF (iter == mp2_env%bse%max_iter) EXIT

            ! corrections t_ia,k = r_ia,k/(d_ia - Ω^k) of the lowest unconverged roots: the rR block
            ! against the basis, then the rL block of the same roots against the extended basis;
            ! both fit the room, which is odd only when the ceiling is the number of transitions
            CALL select_roots(mp2_env, need, res, MIN(block_size, (m_max - m + 1)/2), selected, nt_sel)
            nt_R = nt_sel
            CALL davidson_corrections(fm_work, selected(1:nt_R), omega, diag, fm_b, m + 1)
            CALL extend_basis_os(mv_env, fm_b, fm_Mb, fm_Kb, fm_scratch, m, nt_R, n_kernel, n_dependent)
            nt_L = MIN(nt_sel, m_max - m - nt_R)
            CALL davidson_corrections(fm_work, selected(1:nt_L), omega, diag, fm_b, m + nt_R + 1, n_act)
            CALL extend_basis_os(mv_env, fm_b, fm_Mb, fm_Kb, fm_scratch, m + nt_R, nt_L, n_kernel, n_dependent)
            nt = nt_R + nt_L
            IF (nt == 0 .AND. stalled) THEN
               CPABORT("BSE Davidson: correction vectors vanished before convergence")
            END IF
            stalled = nt == 0

            IF (nt > 0) THEN
               ! new columns of the symmetric Mr and Kr
               ALLOCATE (Mr_new(m + nt, nt), Kr_new(m + nt, nt))
               CALL subspace_gram(fm_b, 1, m + nt, fm_Mb, m + 1, nt, para_env, Mr_new)
               CALL subspace_gram(fm_b, 1, m + nt, fm_Kb, m + 1, nt, para_env, Kr_new)
               CALL reduced_extend(fm_red_M, Mr_new, m, nt)
               CALL reduced_extend(fm_red_K, Kr_new, m, nt)
               DEALLOCATE (Mr_new, Kr_new)
            END IF

            ! the pairs of this iteration are the second block of restart candidates
            C_cand(:, 2*n_act + 1:4*n_act) = 0.0_dp
            C_cand(1:m, 2*n_act + 1:3*n_act) = R(1:m, :)
            C_cand(1:m, 3*n_act + 1:4*n_act) = L(1:m, :)
            has_prev = .TRUE.
            omega_prev(:) = omega(:)
            m = m + nt
         END DO
      END IF

      IF (abba_status == abba_ok) THEN
         IF (.NOT. all_conv) CALL abort_unconverged(omega, res, conv, need, n_req, unit_nr)
         IF (mp2_env%bse%bse_debug_print) CALL print_tracked_roots(omega, res, n_want, n_act, unit_nr)

         IF (mp2_env%bse%bse_debug_print) THEN
            dev = orthonormality_deviation(fm_b, fm_b, m, block_size, para_env)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A10,T13,A,T63,ES18.6)') 'BSE|DEBUG|', &
                  'Max deviation of the basis from orthonormality', dev
            END IF
         END IF

         ! the smallest eigenvalue of A-B on the span of b is that of the last Kr
         ab_margin = kr_min

         IF (ALLOCATED(theta_sub)) CALL check_guess_bound(theta_sub, omega, n_req, unit_nr)
         CALL print_summary(iter, n_kernel, n_restart, n_want, n_req, n_act, n_ov, unit_nr, ab_margin, n_dependent)

         ! X_ia^n = sum_m b_ia,m (R_mn + L_mn)/2, Y_ia^n = sum_m b_ia,m (R_mn - L_mn)/2
         ALLOCATE (exc_ens(n_want))
         exc_ens(:) = omega(1:n_want)
         CALL bse_matvec_vector_struct(mv_env, n_want, fm_struct)
         CALL cp_fm_create(fm_X, fm_struct, name="fm_X_davidson", set_zero=.TRUE.)
         CALL cp_fm_create(fm_Y, fm_struct, name="fm_Y_davidson", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct)
         C_work(1:m, 1:n_want) = 0.5_dp*(R(1:m, 1:n_want) + L(1:m, 1:n_want))
         CALL subspace_rotate(fm_b, m, C_work(1:m, 1:n_want), n_want, fm_X, 1.0_dp, 0.0_dp)
         C_work(1:m, 1:n_want) = 0.5_dp*(R(1:m, 1:n_want) - L(1:m, 1:n_want))
         CALL subspace_rotate(fm_b, m, C_work(1:m, 1:n_want), n_want, fm_Y, 1.0_dp, 0.0_dp)
      END IF

      CALL cp_fm_release(fm_b)
      CALL cp_fm_release(fm_Mb)
      CALL cp_fm_release(fm_Kb)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_scratch)
      CALL cp_fm_release(fm_red_M)
      CALL cp_fm_release(fm_red_K)
      CALL cp_blacs_env_release(diag_blacs_env)
      DEALLOCATE (R, L, C_work, C_cand, Q, omega, omega_prev, res, res_L, dE, conv, need, selected, diag)
      IF (ALLOCATED(theta_sub)) DEALLOCATE (theta_sub)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_abba_os

! **************************************************************************************************
!> \brief Makes the nt columns T that follow column m of V orthonormal in the (A-B) inner product and
!>        stores their products: T_ia,k <- T_ia,k - sum_m V_ia,m sum_jb W_jb,m T_jb,k, Euclidean
!>        orthonormalisation of T, W_ia,k = sum_jb (A-B)_ia,jb T_jb,k, T <- T R^-1 and W <- W R^-1 with
!>        sum_p R_pk R_pl = sum_ia T_ia,k W_ia,l, and (M W)_ia,k = sum_jb (A+B)_ia,jb W_jb,k
!> \param mv_env ...
!> \param fm_V ...
!> \param fm_W ...
!> \param fm_MW ...
!> \param fm_scratch scratch of at least nt columns
!> \param m ...
!> \param nt on exit the number of columns kept
!> \param n_kernel counter of the kernel applications
!> \param n_dependent counter of the correction vectors dropped as dependent
!> \param abba_status ...
!> \param ab_margin smallest eigenvalue of T^T (A-B) T if that matrix is not positive definite
! **************************************************************************************************
   SUBROUTINE extend_basis_mk(mv_env, fm_V, fm_W, fm_MW, fm_scratch, m, nt, n_kernel, n_dependent, abba_status, &
                              ab_margin)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V, fm_W, fm_MW, fm_scratch
      INTEGER, INTENT(IN)                                :: m
      INTEGER, INTENT(INOUT)                             :: nt, n_kernel, n_dependent, abba_status
      REAL(KIND=dp), INTENT(INOUT)                       :: ab_margin

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'extend_basis_mk'

      INTEGER                                            :: handle, n_dropped
      LOGICAL                                            :: ok

      CALL timeset(routineN, handle)

      ! T_ia,k <- T_ia,k - sum_m V_ia,m sum_jb W_jb,m T_jb,k, then Euclidean orthonormalisation
      CALL project_out(fm_V, m + 1, nt, fm_V, m, mv_env%para_env, fm_W)
      CALL drop_small_columns(fm_V, m + 1, nt, mv_env%para_env)
      CALL cholqr2(fm_V, m + 1, nt, mv_env%para_env, ok, n_dropped)
      IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")
      n_dependent = n_dependent + n_dropped

      IF (nt > 0) THEN
         ! W_ia,k = sum_jb (A-B)_ia,jb T_jb,k: A T into W, B T into the scratch
         CALL bse_matvec_apply(mv_env, fm_V, m + 1, nt, fm_W, fm_scratch, 1)
         CALL columns_axpy(-1.0_dp, fm_scratch, 1, fm_W, m + 1, nt)
         n_kernel = n_kernel + nt

         ! T <- T R^-1 and W <- W R^-1 with sum_p R_pk R_pl = sum_ia T_ia,k W_ia,l
         CALL cholqr2_metric(fm_V, fm_W, m + 1, nt, mv_env%para_env, abba_status, ab_margin)

         IF (abba_status == abba_ok) THEN
            IF (extension_deviation(fm_V, fm_W, m, nt, mv_env%para_env) > ortho_tol) THEN
               CPABORT("BSE Davidson: the basis lost orthonormality in the (A-B) inner product")
            END IF
            ! (M W)_ia,k = sum_jb (A+B)_ia,jb W_jb,k: A W into MW, B W into the scratch
            CALL bse_matvec_apply(mv_env, fm_W, m + 1, nt, fm_MW, fm_scratch, 1)
            CALL columns_axpy(1.0_dp, fm_scratch, 1, fm_MW, m + 1, nt)
            n_kernel = n_kernel + nt
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE extend_basis_mk

! **************************************************************************************************
!> \brief Two Cholesky QR passes in the inner product of a symmetric K with the products KT at hand:
!>        G_kl = sum_ia T_ia,k (KT)_ia,l = sum_p R_pk R_pl, T <- T R^-1, KT <- KT R^-1. A failed
!>        factorization of a G with a negative eigenvalue reports K as indefinite, any other failure
!>        is repeated once with a diagonal shift and a third pass. The Gram is the only communication:
!>        R is factorised on one rank and broadcast, the triangular solve runs on the local rows
!> \param fm_T ...
!> \param fm_KT ...
!> \param first_col ...
!> \param nt ...
!> \param para_env ...
!> \param abba_status ...
!> \param eig_min smallest eigenvalue of G, set if it is negative
! **************************************************************************************************
   SUBROUTINE cholqr2_metric(fm_T, fm_KT, first_col, nt, para_env, abba_status, eig_min)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T, fm_KT
      INTEGER, INTENT(IN)                                :: first_col, nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(INOUT)                             :: abba_status
      REAL(KIND=dp), INTENT(INOUT)                       :: eig_min

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cholqr2_metric'

      INTEGER                                            :: handle, info, ipass, k, n_ov, n_pass, &
                                                            nrow_local
      REAL(KIND=dp)                                      :: shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G, R

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local, nrow_global=n_ov)
      ALLOCATE (G(nt, nt), R(nt, nt), eigval(nt))

      n_pass = 2
      ipass = 0
      DO WHILE (ipass < n_pass)
         ipass = ipass + 1
         CALL subspace_gram(fm_T, first_col, nt, fm_KT, first_col, nt, para_env, R)
         G(:, :) = 0.5_dp*(R(:, :) + TRANSPOSE(R(:, :)))
         R(:, :) = G(:, :)
         info = 0
         IF (para_env%is_source()) CALL DPOTRF('U', nt, R, nt, info)
         CALL para_env%bcast(info)
         IF (info /= 0) THEN
            IF (n_pass == 3) THEN
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson: orthonormalisation in the (A-B) inner product broke down")
            END IF
            R(:, :) = G(:, :)
            eigval(:) = 0.0_dp
            IF (para_env%is_source()) CALL diamat_all(R, eigval)
            CALL para_env%bcast(eigval)
            IF (eigval(1) < 0.0_dp) THEN
               abba_status = abba_indefinite
               eig_min = eigval(1)
               EXIT
            END IF
            n_pass = 3
            shift = 11.0_dp*(REAL(n_ov, dp)*REAL(nt, dp) + REAL(nt, dp)*REAL(nt + 1, dp))* &
                    EPSILON(1.0_dp)*NORM2(G)
            R(:, :) = G(:, :)
            DO k = 1, nt
               R(k, k) = R(k, k) + shift
            END DO
            IF (para_env%is_source()) CALL DPOTRF('U', nt, R, nt, info)
            CALL para_env%bcast(info)
            IF (info /= 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson: orthonormalisation in the (A-B) inner product broke down")
            END IF
         END IF
         CALL para_env%bcast(R)
         IF (nrow_local > 0) THEN
            CALL DTRSM('R', 'U', 'N', 'N', nrow_local, nt, 1.0_dp, R, nt, &
                       fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1))
            CALL DTRSM('R', 'U', 'N', 'N', nrow_local, nt, 1.0_dp, R, nt, &
                       fm_KT%local_data(:, first_col:first_col + nt - 1), SIZE(fm_KT%local_data, 1))
         END IF
      END DO

      DEALLOCATE (G, R, eigval)
      CALL timestop(handle)

   END SUBROUTINE cholqr2_metric

! **************************************************************************************************
!> \brief Paired reduced problem on the orthonormal basis b_ia,m: with Mr_mn = sum_ia b_ia,m (M b)_ia,n
!>        and Kr_mn = sum_ia b_ia,m (K b)_ia,n, sum_n Mr_mn R_nk = Ω_k L_mk and sum_n Kr_mn L_nk = Ω_k R_mk.
!>        Solved as in initial_guess_subblock: Kr = U diag(κ) U^T, S = U diag(√κ) U^T,
!>        sum_nop S_mn Mr_no S_op T_pk = λ_k T_mk, Ω_k = √λ_k, R_mk = sum_n S_mn T_nk/√Ω_k,
!>        L_mk = sum_n Mr_mn R_nk/Ω_k. κ_1 <= 0 reports K as indefinite on the span of b.
!>        Every m x m step runs on the grid of blacs_env: the two eigensolves through
!>        choose_eigv_solver and the four products through parallel_gemm; only the pair
!>        coefficients come back replicated
!> \param fm_red_M Mr in its live block, on the grid of blacs_env
!> \param fm_red_K Kr in its live block, on the grid of blacs_env
!> \param m ...
!> \param n_act number of pairs formed
!> \param para_env ...
!> \param omega Ω_k in Hartree
!> \param R coefficients of y_k in column k, zero below row m
!> \param L coefficients of x_k in column k, zero below row m
!> \param kr_min κ_1, the smallest eigenvalue of Kr
!> \param abba_status abba_indefinite when κ_1 <= 0; omega, R and L are then untouched
!> \param blacs_env grid the reduced problem is solved on
! **************************************************************************************************
   SUBROUTINE solve_reduced_paired(fm_red_M, fm_red_K, m, n_act, para_env, omega, R, L, kr_min, &
                                   abba_status, blacs_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red_M, fm_red_K
      INTEGER, INTENT(IN)                                :: m, n_act
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: omega
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: R, L
      REAL(KIND=dp), INTENT(OUT)                         :: kr_min
      INTEGER, INTENT(INOUT)                             :: abba_status
      TYPE(cp_blacs_env_type), INTENT(IN), TARGET        :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_reduced_paired'

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: kappa, lambda
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_H, fm_L, fm_Mr, fm_R, fm_S, fm_T, &
                                                            fm_U, fm_work

      CALL timeset(routineN, handle)
      CPASSERT(n_act <= m)

      ALLOCATE (kappa(m), lambda(m))
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=m, ncol_global=m)
      CALL cp_fm_create(fm_U, fm_struct, name="bse_paired_U")
      CALL cp_fm_create(fm_S, fm_struct, name="bse_paired_S")
      CALL cp_fm_create(fm_Mr, fm_struct, name="bse_paired_Mr")
      CALL cp_fm_create(fm_H, fm_struct, name="bse_paired_H")
      CALL cp_fm_create(fm_T, fm_struct, name="bse_paired_T")
      CALL cp_fm_create(fm_work, fm_struct, name="bse_paired_work")
      CALL cp_fm_struct_release(fm_struct)

      ! Kr = U diag(κ) U^T, on the symmetrised live block
      CALL cp_fm_to_fm_submat(fm_red_K, fm_work, m, m, 1, 1, 1, 1)
      CALL cp_fm_transpose(fm_work, fm_U)
      CALL cp_fm_scale_and_add(0.5_dp, fm_work, 0.5_dp, fm_U)
      kappa(:) = 0.0_dp
      CALL choose_eigv_solver(fm_work, fm_U, kappa)
      kr_min = kappa(1)

      IF (kappa(1) <= 0.0_dp) THEN
         abba_status = abba_indefinite
      ELSE
         ! S = U diag(√κ) U^T
         CALL cp_fm_to_fm(fm_U, fm_work)
         CALL cp_fm_column_scale(fm_work, SQRT(kappa))
         CALL parallel_gemm('N', 'T', m, m, m, 1.0_dp, fm_work, fm_U, 0.0_dp, fm_S)

         ! H_mp = sum_no S_mn Mr_no S_op, then T holds its eigenvectors
         CALL cp_fm_to_fm_submat(fm_red_M, fm_Mr, m, m, 1, 1, 1, 1)
         CALL cp_fm_transpose(fm_Mr, fm_T)
         CALL cp_fm_scale_and_add(0.5_dp, fm_Mr, 0.5_dp, fm_T)
         CALL parallel_gemm('N', 'N', m, m, m, 1.0_dp, fm_Mr, fm_S, 0.0_dp, fm_work)
         CALL parallel_gemm('N', 'N', m, m, m, 1.0_dp, fm_S, fm_work, 0.0_dp, fm_H)
         lambda(:) = 0.0_dp
         CALL choose_eigv_solver(fm_H, fm_T, lambda)
         IF (lambda(1) <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "BSE Davidson (OLSEN_STRATMANN): negative squared excitation energy. Matrix "// &
                          "(A+B) is not positive definite on the trial space, or the basis lost its "// &
                          "orthonormality.")
         END IF
         omega(1:n_act) = SQRT(lambda(1:n_act))

         ! R_mk = sum_n S_mn T_nk/√Ω_k, L_mk = sum_n Mr_mn R_nk/Ω_k
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=m, ncol_global=n_act)
         CALL cp_fm_create(fm_R, fm_struct, name="bse_paired_R")
         CALL cp_fm_create(fm_L, fm_struct, name="bse_paired_L")
         CALL cp_fm_struct_release(fm_struct)
         CALL parallel_gemm('N', 'N', m, n_act, m, 1.0_dp, fm_S, fm_T, 0.0_dp, fm_R)
         CALL cp_fm_column_scale(fm_R, 1.0_dp/SQRT(omega(1:n_act)))
         CALL parallel_gemm('N', 'N', m, n_act, m, 1.0_dp, fm_Mr, fm_R, 0.0_dp, fm_L)
         CALL cp_fm_column_scale(fm_L, 1.0_dp/omega(1:n_act))

         R(:, :) = 0.0_dp
         L(:, :) = 0.0_dp
         CALL cp_fm_get_submatrix(fm_R, R(1:m, 1:n_act))
         CALL cp_fm_get_submatrix(fm_L, L(1:m, 1:n_act))
         CALL cp_fm_release(fm_R)
         CALL cp_fm_release(fm_L)
      END IF

      CALL cp_fm_release(fm_U)
      CALL cp_fm_release(fm_S)
      CALL cp_fm_release(fm_Mr)
      CALL cp_fm_release(fm_H)
      CALL cp_fm_release(fm_T)
      CALL cp_fm_release(fm_work)
      DEALLOCATE (kappa, lambda)
      CALL timestop(handle)

   END SUBROUTINE solve_reduced_paired

! **************************************************************************************************
!> \brief Makes the nt columns T that follow column m of b orthonormal to b and to each other,
!>        T_ia,k <- T_ia,k - sum_m b_ia,m sum_jb b_jb,m T_jb,k and Cholesky QR, and stores the products
!>        (M b)_ia,k = sum_jb (A+B)_ia,jb b_jb,k, (K b)_ia,k = sum_jb (A-B)_ia,jb b_jb,k of the kept
!>        columns from one application of A and B
!> \param mv_env ...
!> \param fm_b ...
!> \param fm_Mb ...
!> \param fm_Kb ...
!> \param fm_scratch scratch of at least nt columns
!> \param m ...
!> \param nt on exit the number of columns kept
!> \param n_kernel counter of the kernel applications
!> \param n_dependent counter of the correction vectors dropped as dependent
! **************************************************************************************************
   SUBROUTINE extend_basis_os(mv_env, fm_b, fm_Mb, fm_Kb, fm_scratch, m, nt, n_kernel, n_dependent)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_b, fm_Mb, fm_Kb, fm_scratch
      INTEGER, INTENT(IN)                                :: m
      INTEGER, INTENT(INOUT)                             :: nt, n_kernel, n_dependent

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'extend_basis_os'

      INTEGER                                            :: handle, n_dropped
      LOGICAL                                            :: ok

      CALL timeset(routineN, handle)

      CALL project_out(fm_b, m + 1, nt, fm_b, m, mv_env%para_env)
      CALL drop_small_columns(fm_b, m + 1, nt, mv_env%para_env)
      CALL cholqr2(fm_b, m + 1, nt, mv_env%para_env, ok, n_dropped)
      IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")
      n_dependent = n_dependent + n_dropped

      IF (nt > 0) THEN
         IF (extension_deviation(fm_b, fm_b, m, nt, mv_env%para_env) > ortho_tol) THEN
            CPABORT("BSE Davidson: the basis lost orthonormality")
         END IF
         ! A b into Mb and B b into the scratch, then K b = A b - B b and M b = A b + B b
         CALL bse_matvec_apply(mv_env, fm_b, m + 1, nt, fm_Mb, fm_scratch, 1)
         CALL cp_fm_to_fm(fm_Mb, fm_Kb, nt, m + 1, m + 1)
         CALL columns_axpy(-1.0_dp, fm_scratch, 1, fm_Kb, m + 1, nt)
         CALL columns_axpy(1.0_dp, fm_scratch, 1, fm_Mb, m + 1, nt)
         n_kernel = n_kernel + nt
      END IF

      CALL timestop(handle)

   END SUBROUTINE extend_basis_os

! **************************************************************************************************
!> \brief Convergence flag of every tracked root from CONVERGENCE_CRITERION
!> \param mp2_env ...
!> \param dE energy changes in Hartree
!> \param res residual norms
!> \param conv ...
! **************************************************************************************************
   SUBROUTINE convergence_flags(mp2_env, dE, res, conv)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: dE, res
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: conv

      INTEGER                                            :: k

      DO k = 1, SIZE(conv)
         SELECT CASE (mp2_env%bse%convergence_criterion)
         CASE (bse_iter_en_cond)
            conv(k) = dE(k) < mp2_env%bse%eps_energy
         CASE (bse_iter_res_cond)
            conv(k) = res(k) < mp2_env%bse%eps_res
         CASE (bse_iter_or_cond)
            conv(k) = dE(k) < mp2_env%bse%eps_energy .OR. res(k) < mp2_env%bse%eps_res
         CASE (bse_iter_and_cond)
            conv(k) = dE(k) < mp2_env%bse%eps_energy .AND. res(k) < mp2_env%bse%eps_res
         END SELECT
      END DO

   END SUBROUTINE convergence_flags

! **************************************************************************************************
!> \brief Last root that has to converge: a multiplet that is cut by NUM_EXC_EN converges as a whole
!> \param energies ...
!> \param n_want ...
!> \param n_act ...
!> \return ...
! **************************************************************************************************
   FUNCTION last_required_root(energies, n_want, n_act) RESULT(n_req)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: energies
      INTEGER, INTENT(IN)                                :: n_want, n_act
      INTEGER                                            :: n_req

      n_req = n_want
      DO WHILE (n_req < n_act)
         IF (energies(n_req + 1) - energies(n_req) >= deg_thresh) EXIT
         n_req = n_req + 1
      END DO

   END FUNCTION last_required_root

! **************************************************************************************************
!> \brief Roots that still need correction vectors: the unconverged ones up to n_req, and a tracked
!>        root k above them as long as E_k - |r_k| < E_n_req. An eigenvalue of a symmetric matrix
!>        lies within |r_k| of E_k, so such a root can still belong to the requested states
!> \param conv ...
!> \param energies ...
!> \param res ...
!> \param n_req ...
!> \param need ...
! **************************************************************************************************
   SUBROUTINE roots_in_need(conv, energies, res, n_req, need)

      LOGICAL, DIMENSION(:), INTENT(IN)                  :: conv
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: energies, res
      INTEGER, INTENT(IN)                                :: n_req
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: need

      INTEGER                                            :: k

      DO k = 1, SIZE(need)
         need(k) = .NOT. conv(k)
         IF (k > n_req) need(k) = need(k) .AND. energies(k) - res(k) < energies(n_req)
      END DO

   END SUBROUTINE roots_in_need

! **************************************************************************************************
!> \brief The lowest roots in need, which receive a correction vector
!> \param mp2_env ...
!> \param need ...
!> \param res ...
!> \param n_max largest number of roots
!> \param selected ...
!> \param nt number of selected roots
! **************************************************************************************************
   SUBROUTINE select_roots(mp2_env, need, res, n_max, selected, nt)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: need
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: res
      INTEGER, INTENT(IN)                                :: n_max
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: selected
      INTEGER, INTENT(OUT)                               :: nt

      INTEGER                                            :: k

      selected(:) = 0
      nt = 0
      DO k = 1, SIZE(need)
         IF (.NOT. need(k)) CYCLE
         ! a root below EPS_RES only waits for its energy criterion
         IF (mp2_env%bse%convergence_criterion == bse_iter_and_cond .AND. &
             res(k) < mp2_env%bse%eps_res) CYCLE
         IF (nt == n_max) EXIT
         nt = nt + 1
         selected(nt) = k
      END DO

   END SUBROUTINE select_roots

! **************************************************************************************************
!> \brief Title and column heads of the iteration table
!> \param title ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_iteration_header(title, unit_nr)

      CHARACTER(LEN=*), INTENT(IN)                       :: title
      INTEGER, INTENT(IN)                                :: unit_nr

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', title
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
            '|r| in a.u., energies in eV, t the wall time of the iteration'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A5,T14,A8,T24,A6,T33,A12,T47,A12,T61,A12,T74,A7)') 'BSE|', &
            'Iter.', 'Z-space', 'Conv.', 'Max |r|', 'Max dE', 'Lowest E', 't (s)'
      END IF

   END SUBROUTINE print_iteration_header

! **************************************************************************************************
!> \brief One line of the iteration table
!> \param iter ...
!> \param m subspace dimension
!> \param n_conv ...
!> \param max_res ...
!> \param max_dE in Hartree
!> \param lowest_E in Hartree
!> \param t_iter wall time of the iteration in seconds
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_iteration(iter, m, n_conv, max_res, max_dE, lowest_E, t_iter, unit_nr)

      INTEGER, INTENT(IN)                                :: iter, m, n_conv
      REAL(KIND=dp), INTENT(IN)                          :: max_res, max_dE, lowest_E, t_iter
      INTEGER, INTENT(IN)                                :: unit_nr

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,I5,T14,I8,T24,I6,T33,ES12.4,T47,ES12.4,T61,F12.6,T74,F7.2)') &
            'BSE|', iter, m, n_conv, max_res, max_dE*evolt, lowest_E*evolt, t_iter
      END IF

   END SUBROUTINE print_iteration

! **************************************************************************************************
!> \brief Table of the requested roots and of the tracked roots in need, then the abort for MAX_ITER
!> \param energies in Hartree
!> \param res ...
!> \param conv ...
!> \param need ...
!> \param n_req ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE abort_unconverged(energies, res, conv, need, n_req, unit_nr)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: energies, res
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: conv, need
      INTEGER, INTENT(IN)                                :: n_req, unit_nr

      INTEGER                                            :: k

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A12,T30,A11,T50,A10,T67,A14)') 'BSE|', &
            'Excitation n', 'Energy (eV)', '|r| (a.u.)', 'Converged'
         DO k = 1, SIZE(need)
            IF (k > n_req .AND. .NOT. need(k)) CYCLE
            WRITE (unit_nr, '(T2,A4,T7,I12,T27,F14.6,T46,ES14.4,T67,L14)') 'BSE|', &
               k, energies(k)*evolt, res(k), conv(k)
         END DO
      END IF
      CPABORT("BSE Davidson: MAX_ITER reached before convergence")

   END SUBROUTINE abort_unconverged

! **************************************************************************************************
!> \brief Debug table of all tracked roots of a converged solve
!> \param energies in Hartree
!> \param res ...
!> \param n_want ...
!> \param n_act ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_tracked_roots(energies, res, n_want, n_act, unit_nr)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: energies, res
      INTEGER, INTENT(IN)                                :: n_want, n_act, unit_nr

      INTEGER                                            :: k

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A10,T13,A13,T27,A9,T45,A18,T67,A14)') 'BSE|DEBUG|', &
            'Tracked state', 'Requested', 'Energy (eV)', '|r| (a.u.)'
         DO k = 1, n_act
            WRITE (unit_nr, '(T2,A10,T13,I13,T27,L9,T45,F18.10,T67,ES14.4)') 'BSE|DEBUG|', &
               k, k <= n_want, energies(k)*evolt, res(k)
         END DO
      END IF

   END SUBROUTINE print_tracked_roots

! **************************************************************************************************
!> \brief Counters of a converged solve, the multiplet note and the warning for a degenerate group
!>        that reaches the last buffer state
!> \param iter ...
!> \param n_kernel ...
!> \param n_restart ...
!> \param n_want ...
!> \param n_req ...
!> \param n_act ...
!> \param n_ov ...
!> \param unit_nr ...
!> \param ab_margin smallest eigenvalue of A-B on the trial space in Hartree
!> \param n_dependent correction vectors dropped as dependent, printed when not zero
! **************************************************************************************************
   SUBROUTINE print_summary(iter, n_kernel, n_restart, n_want, n_req, n_act, n_ov, unit_nr, ab_margin, &
                            n_dependent)

      INTEGER, INTENT(IN)                                :: iter, n_kernel, n_restart, n_want, &
                                                            n_req, n_act, n_ov, unit_nr
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: ab_margin
      INTEGER, INTENT(IN), OPTIONAL                      :: n_dependent

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Davidson iterations', iter
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Kernel applications', n_kernel
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Thick restarts', n_restart
         IF (PRESENT(n_dependent)) THEN
            IF (n_dependent > 0) THEN
               WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Correction vectors dropped as dependent', n_dependent
            END IF
         END IF
         IF (PRESENT(ab_margin)) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T65,F16.6)') 'BSE|', &
               'A-B margin on the trial space (eV)', ab_margin*evolt
         END IF
         IF (n_req > n_want) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
               'Note: the requested states end inside a degenerate multiplet.'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      IF (n_req > n_want .AND. n_req == n_act .AND. n_act < n_ov) THEN
         CALL cp_warn(__LOCATION__, &
                      "BSE Davidson: the degenerate group at NUM_EXC_EN reaches the last buffer state, "// &
                      "further members may lie above it. Increase BSE_ITERAT%NUM_BUFFER_STATES.")
      END IF

   END SUBROUTINE print_summary

! **************************************************************************************************
!> \brief Number of wanted and tracked states and the block size, clamped to the problem size
!> \param mp2_env ...
!> \param n_ov ...
!> \param n_want ...
!> \param n_act wanted plus buffer states
!> \param block_size ...
! **************************************************************************************************
   SUBROUTINE davidson_sizes(mp2_env, n_ov, n_want, n_act, block_size)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: n_ov
      INTEGER, INTENT(OUT)                               :: n_want, n_act, block_size

      IF (mp2_env%bse%num_exc_en < 1) CPABORT("BSE_ITERAT%NUM_EXC_EN must be at least 1")
      IF (mp2_env%bse%num_buffer_states < 0) CPABORT("BSE_ITERAT%NUM_BUFFER_STATES must not be negative")
      IF (mp2_env%bse%block_size < 1 .AND. mp2_env%bse%block_size /= -1) THEN
         CPABORT("BSE_ITERAT%BLOCK_SIZE must be at least 1, or -1 for the default")
      END IF
      IF (mp2_env%bse%max_iter < 1) CPABORT("BSE_ITERAT%MAX_ITER must be at least 1")

      n_want = mp2_env%bse%num_exc_en
      IF (n_want > n_ov) THEN
         CPWARN("BSE_ITERAT%NUM_EXC_EN exceeds the number of transitions and is reduced to it.")
         n_want = n_ov
      END IF
      n_act = MIN(n_want + mp2_env%bse%num_buffer_states, n_ov)
      block_size = mp2_env%bse%block_size
      ! 32 was within 5 % of the best width on the L8 deck at 200 states; wider passes cost restarts
      IF (block_size == -1) block_size = MIN(32, n_act)
      IF (block_size > n_act) THEN
         CPWARN("BSE_ITERAT%BLOCK_SIZE exceeds the number of tracked states and is reduced to it.")
         block_size = n_act
      END IF

   END SUBROUTINE davidson_sizes

! **************************************************************************************************
!> \brief Largest subspace dimension, from MAX_SUBSPACE_FACTOR or from the free memory per rank
!> \param mp2_env ...
!> \param n_ov ...
!> \param n_act ...
!> \param n_arrays number of stored N_ov x m_max arrays
!> \param para_env ...
!> \param unit_nr ...
!> \param m_max ...
!> \param min_factor smallest subspace in vectors per tracked state, 3 by default
!> \param work_factor work columns per tracked state in the memory estimate, 2 by default
! **************************************************************************************************
   SUBROUTINE subspace_ceiling(mp2_env, n_ov, n_act, n_arrays, para_env, unit_nr, m_max, min_factor, &
                               work_factor)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: n_ov, n_act, n_arrays
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, INTENT(OUT)                               :: m_max
      INTEGER, INTENT(IN), OPTIONAL                      :: min_factor, work_factor

      REAL(KIND=dp), PARAMETER                           :: mem_fraction = 0.5_dp

      CHARACTER(LEN=16)                                  :: fac_str, mem_str
      INTEGER                                            :: min_fac, work_fac
      LOGICAL                                            :: from_memory
      REAL(KIND=dp)                                      :: GB_per_column, mem_avail_GB, n_columns

      min_fac = 3
      IF (PRESENT(min_factor)) min_fac = min_factor
      work_fac = 2
      IF (PRESENT(work_factor)) work_fac = work_factor
      WRITE (fac_str, '(I16)') min_fac

      ! one column of length N_ov, spread over all ranks
      GB_per_column = 8.0E-9_dp*REAL(n_ov, dp)/REAL(para_env%num_pe, dp)
      from_memory = .FALSE.

      IF (mp2_env%bse%max_subspace_factor == -1) THEN
         from_memory = .TRUE.
         CALL mp_mem_avail_per_rank_GB(para_env, mem_avail_GB)
         IF (mem_avail_GB > 0.0_dp) THEN
            ! work array of work_fac*n_act columns besides the n_arrays subspace arrays
            n_columns = mem_fraction*mem_avail_GB/GB_per_column - REAL(work_fac*n_act, dp)
            n_columns = MIN(n_columns/REAL(n_arrays, dp), 20.0_dp*REAL(n_act, dp))
            IF (n_columns < REAL(min_fac*n_act, dp) .AND. min_fac*n_act <= n_ov) THEN
               WRITE (mem_str, '(F16.3)') GB_per_column*REAL((min_fac*n_arrays + work_fac)*n_act, dp)/mem_fraction
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson: not enough memory for the smallest subspace of "// &
                             TRIM(ADJUSTL(fac_str))//" vectors per state. Free memory needed per MPI rank (GB): "// &
                             TRIM(ADJUSTL(mem_str)))
            END IF
            m_max = INT(n_columns)
         ELSE
            ! free memory not detectable on this system
            m_max = min_fac*n_act
         END IF
      ELSE IF (mp2_env%bse%max_subspace_factor >= min_fac) THEN
         m_max = mp2_env%bse%max_subspace_factor*n_act
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "BSE_ITERAT%MAX_SUBSPACE_FACTOR must be at least "//TRIM(ADJUSTL(fac_str))// &
                       " (or -1) for the chosen solver")
      END IF
      m_max = MIN(MAX(m_max, min_fac*n_act), n_ov)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (from_memory) THEN
            ! the memory estimate is capped at 20*n_act, so say which of the two bound it
            IF (m_max >= 20*n_act) THEN
               WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
                  'Maximum subspace dimension (20 x buffered states)', m_max
            ELSE
               WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
                  'Maximum subspace dimension (from free memory)', m_max
            END IF
            WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', &
               'Fraction of the free memory made available', mem_fraction
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
               'Maximum subspace dimension (MAX_SUBSPACE_FACTOR)', m_max
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', &
            'Memory of the subspace arrays per MPI rank (GB)', &
            GB_per_column*REAL(n_arrays*m_max + MIN(work_fac*n_act, m_max), dp)
      END IF

   END SUBROUTINE subspace_ceiling

! **************************************************************************************************
!> \brief Unit vectors on the lowest entries of the diagonal, never cutting a degenerate group
!>        The sorted diagonal is replicated; each rank writes the entries that fall on its own rows
!> \param diag ...
!> \param n_guess ...
!> \param n_max largest number of guess vectors
!> \param deg_thresh_diag ...
!> \param fm_Z ...
!> \param n_used ...
!> \param n_given leading columns that the caller has filled; they replace the lowest unit vectors
! **************************************************************************************************
   SUBROUTINE initial_guess(diag, n_guess, n_max, deg_thresh_diag, fm_Z, n_used, n_given)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: diag
      INTEGER, INTENT(IN)                                :: n_guess, n_max
      REAL(KIND=dp), INTENT(IN)                          :: deg_thresh_diag
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(OUT)                               :: n_used
      INTEGER, INTENT(IN), OPTIONAL                      :: n_given

      INTEGER                                            :: iloc, k, k_first, n_ov, nrow_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_sorted

      n_ov = SIZE(diag)
      ALLOCATE (diag_sorted(n_ov), idx(n_ov))
      diag_sorted(:) = diag(:)
      CALL sort(diag_sorted, n_ov, idx)

      n_used = n_guess
      DO WHILE (n_used < n_max)
         IF (diag_sorted(n_used + 1) - diag_sorted(n_used) >= deg_thresh_diag) EXIT
         n_used = n_used + 1
      END DO

      k_first = 1
      IF (PRESENT(n_given)) k_first = n_given + 1

      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      DO k = k_first, n_used
         DO iloc = 1, nrow_local
            IF (row_indices(iloc) == idx(k)) fm_Z%local_data(iloc, k) = 1.0_dp
         END DO
      END DO

      DEALLOCATE (diag_sorted, idx)

   END SUBROUTINE initial_guess

! **************************************************************************************************
!> \brief Initial guess from the exact A (and B) on the n_sub transitions with the lowest diagonal:
!>        TDA: eigenvectors of A_sub; ABBA: x = K_sub^-1/2 t with K_sub^1/2 M_sub K_sub^1/2 t = Ω^2 t,
!>        K = A-B, M = A+B. The Ritz values θ_k are upper bounds of the true energies (Rayleigh-Ritz for
!>        TDA, Cauchy-type interlacing of Bai and Li for the linear response problem), kept for
!>        check_guess_bound. A degenerate group at the upper boundary is taken as a whole. The block and
!>        its eigenvectors are replicated; each rank writes the entries that fall on its own rows
!> \param mv_env ...
!> \param mp2_env ...
!> \param diag ...
!> \param n_guess ...
!> \param n_max largest number of guess vectors
!> \param deg_thresh_diag ...
!> \param fm_Z ...
!> \param n_used ...
!> \param theta_sub allocated here, the n_sub Ritz values in Hartree
!> \param unit_nr ...
!> \param n_given leading columns that the caller has filled; they replace the lowest guess vectors
!> \param abba_status ABBA problem when present; abba_indefinite if K_sub is not positive definite
!> \param ab_margin smallest eigenvalue of K_sub in Hartree
!> \param n_pair ABBA: partners y_k = K_sub x_k/θ_k of the guess vectors x_k not given by the caller,
!>        written to the n_pair columns after n_used; n_used + n_pair stays within n_max
! **************************************************************************************************
   SUBROUTINE initial_guess_subblock(mv_env, mp2_env, diag, n_guess, n_max, deg_thresh_diag, fm_Z, n_used, &
                                     theta_sub, unit_nr, n_given, abba_status, ab_margin, n_pair)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: diag
      INTEGER, INTENT(IN)                                :: n_guess, n_max
      REAL(KIND=dp), INTENT(IN)                          :: deg_thresh_diag
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(OUT)                               :: n_used
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: theta_sub
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, INTENT(IN), OPTIONAL                      :: n_given
      INTEGER, INTENT(OUT), OPTIONAL                     :: abba_status
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: ab_margin
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_pair

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initial_guess_subblock'

      INTEGER                                            :: handle, iloc, k, k_first, n_ov, n_sub, &
                                                            n_top, nrow_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx, pos
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      LOGICAL                                            :: ok
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_sorted, eig_K
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: A_sub, B_sub, K_mhalf, K_phalf, U, V, Y
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      IF (PRESENT(n_pair)) n_pair = 0
      para_env => mv_env%para_env
      n_ov = SIZE(diag)
      n_sub = MIN(n_ov, MAX(mp2_env%bse%num_guess_transitions, n_guess))
      ALLOCATE (diag_sorted(n_ov), idx(n_ov), pos(n_ov))
      diag_sorted(:) = diag(:)
      CALL sort(diag_sorted, n_ov, idx)

      ALLOCATE (A_sub(n_sub, n_sub), V(n_sub, n_sub), theta_sub(n_sub))
      ok = .TRUE.
      IF (PRESENT(abba_status)) THEN
         abba_status = abba_ok
         ALLOCATE (B_sub(n_sub, n_sub), U(n_sub, n_sub), K_phalf(n_sub, n_sub), K_mhalf(n_sub, n_sub), &
                   eig_K(n_sub))
         CALL bse_matvec_subblock(mv_env, idx(1:n_sub), A_sub, B_sub)
         ! K_sub^±1/2 from its eigen decomposition U diag(κ) U^T
         B_sub(:, :) = A_sub(:, :) - B_sub(:, :)
         CALL solve_replicated(B_sub, n_sub, para_env, eig_K, U)
         ab_margin = eig_K(1)
         IF (eig_K(1) <= 0.0_dp) THEN
            abba_status = abba_indefinite
            ok = .FALSE.
            n_used = 0
         ELSE
            DO k = 1, n_sub
               K_phalf(:, k) = U(:, k)*SQRT(eig_K(k))
               K_mhalf(:, k) = U(:, k)/SQRT(eig_K(k))
            END DO
            K_phalf(:, :) = MATMUL(K_phalf, TRANSPOSE(U))
            K_mhalf(:, :) = MATMUL(K_mhalf, TRANSPOSE(U))
            ! M_sub = 2 A_sub - K_sub, then H = K^1/2 M K^1/2 and x = K^-1/2 t
            A_sub(:, :) = 2.0_dp*A_sub(:, :) - B_sub(:, :)
            A_sub(:, :) = MATMUL(K_phalf, MATMUL(A_sub, K_phalf))
            CALL solve_replicated(A_sub, n_sub, para_env, theta_sub, U)
            IF (theta_sub(1) <= 0.0_dp) THEN
               CPABORT("BSE Davidson: the guess block gives a non-positive squared excitation energy")
            END IF
            theta_sub(:) = SQRT(theta_sub(:))
            V(:, :) = MATMUL(K_mhalf, U)
            IF (PRESENT(n_pair)) THEN
               ! partners y_k = K_sub x_k/θ_k, B_sub holds K_sub
               ALLOCATE (Y(n_sub, n_sub))
               Y(:, :) = MATMUL(B_sub, V)
               DO k = 1, n_sub
                  Y(:, k) = Y(:, k)/theta_sub(k)
               END DO
            END IF
         END IF
         DEALLOCATE (B_sub, U, K_phalf, K_mhalf, eig_K)
      ELSE
         CALL bse_matvec_subblock(mv_env, idx(1:n_sub), A_sub)
         CALL solve_replicated(A_sub, n_sub, para_env, theta_sub, V)
      END IF

      IF (ok) THEN
         k_first = 1
         IF (PRESENT(n_given)) k_first = n_given + 1

         n_used = n_guess
         n_top = MIN(n_max, n_sub)
         ! the partners of the columns k_first..n_used share the n_max columns with the guess
         IF (PRESENT(n_pair)) n_top = MIN(n_top, (n_max + k_first - 1)/2)
         DO WHILE (n_used < n_top)
            IF (theta_sub(n_used + 1) - theta_sub(n_used) >= deg_thresh_diag) EXIT
            n_used = n_used + 1
         END DO

         pos(:) = 0
         DO k = 1, n_sub
            pos(idx(k)) = k
         END DO
         CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
         DO k = k_first, n_used
            DO iloc = 1, nrow_local
               IF (pos(row_indices(iloc)) > 0) fm_Z%local_data(iloc, k) = V(pos(row_indices(iloc)), k)
            END DO
         END DO
         IF (PRESENT(n_pair)) THEN
            n_pair = MIN(n_used - k_first + 1, n_max - n_used)
            DO k = 1, n_pair
               DO iloc = 1, nrow_local
                  IF (pos(row_indices(iloc)) > 0) THEN
                     fm_Z%local_data(iloc, n_used + k) = Y(pos(row_indices(iloc)), k_first + k - 1)
                  END IF
               END DO
            END DO
         END IF

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Transitions in the exact guess block', n_sub
            IF (mp2_env%bse%bse_debug_print) THEN
               WRITE (unit_nr, '(T2,A10,T13,A,T67,F14.6)') 'BSE|DEBUG|', &
                  'Lowest Ritz value of the guess block (eV)', theta_sub(1)*evolt
            END IF
         END IF
      END IF

      DEALLOCATE (diag_sorted, idx, pos, A_sub, V)
      IF (ALLOCATED(Y)) DEALLOCATE (Y)

      CALL timestop(handle)

   END SUBROUTINE initial_guess_subblock

! **************************************************************************************************
!> \brief A converged energy above the Ritz value of the guess block with the same index is a
!>        provable miss: both are upper bounds of the true eigenvalue with that index
!> \param theta_sub ...
!> \param energies converged energies in Hartree
!> \param n_req ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE check_guess_bound(theta_sub, energies, n_req, unit_nr)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: theta_sub, energies
      INTEGER, INTENT(IN)                                :: n_req, unit_nr

      REAL(KIND=dp)                                      :: excess

      excess = MAXVAL(energies(1:n_req) - theta_sub(1:n_req))
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A,T65,F16.6)') 'BSE|', &
            'Largest excess over the guess-block bound (eV)', excess*evolt
      END IF
      IF (excess > deg_thresh) THEN
         CALL cp_warn(__LOCATION__, &
                      "A converged BSE state lies above the upper bound given by the guess block, "// &
                      "so a lower state was missed. Raise BSE_ITERAT%NUM_GUESS_TRANSITIONS or "// &
                      "BSE_ITERAT%NUM_BUFFER_STATES.")
      END IF

   END SUBROUTINE check_guess_bound

! **************************************************************************************************
!> \brief Eigenpairs of the symmetrised live block of the reduced matrix on the grid of blacs_env:
!>        theta ascending, and the first n_vec eigenvectors gathered into C(1:m, 1:n_vec)
!> \param fm_red ...
!> \param m ...
!> \param n_vec eigenvectors gathered; the drivers never read more than n_act
!> \param para_env ...
!> \param blacs_env ...
!> \param theta ...
!> \param C ...
! **************************************************************************************************
   SUBROUTINE solve_reduced(fm_red, m, n_vec, para_env, blacs_env, theta, C)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      INTEGER, INTENT(IN)                                :: m, n_vec
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env
      TYPE(cp_blacs_env_type), INTENT(IN), TARGET        :: blacs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: theta
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: C

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_reduced'

      INTEGER                                            :: handle, n_get
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      TYPE(cp_fm_type)                                   :: fm_a, fm_t

      CALL timeset(routineN, handle)

      CALL reduced_live_block(fm_red, m, para_env, blacs_env, "bse_reduced", fm_a)
      CALL cp_fm_create(fm_t, fm_a%matrix_struct, name="bse_reduced_vectors")
      ! a <- (a + a^T)/2: the eigensolver reads one triangle
      CALL cp_fm_transpose(fm_a, fm_t)
      CALL cp_fm_scale_and_add(0.5_dp, fm_a, 0.5_dp, fm_t)
      ALLOCATE (eigval(m))
      eigval(:) = 0.0_dp
      CALL choose_eigv_solver(fm_a, fm_t, eigval)

      theta(:) = 0.0_dp
      theta(1:m) = eigval(:)
      C(:, :) = 0.0_dp
      n_get = MIN(m, n_vec)
      CALL cp_fm_get_submatrix(fm_t, C(1:m, 1:n_get), 1, 1, m, n_get)

      DEALLOCATE (eigval)
      CALL cp_fm_release(fm_a)
      CALL cp_fm_release(fm_t)

      CALL timestop(handle)

   END SUBROUTINE solve_reduced

! **************************************************************************************************
!> \brief Eigenpairs of a small replicated symmetric matrix, symmetrised, solved on one rank and
!>        broadcast: the fixed-size blocks of the guess and of cholqr2, never the reduced matrix
!> \param G ...
!> \param m ...
!> \param para_env ...
!> \param theta ascending
!> \param C all m eigenvectors in the columns
! **************************************************************************************************
   SUBROUTINE solve_replicated(G, m, para_env, theta, C)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: G
      INTEGER, INTENT(IN)                                :: m
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: theta
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: C

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_replicated'

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G_sym

      CALL timeset(routineN, handle)

      ALLOCATE (G_sym(m, m), eigval(m))
      G_sym(:, :) = 0.5_dp*(G(1:m, 1:m) + TRANSPOSE(G(1:m, 1:m)))
      eigval(:) = 0.0_dp
      IF (para_env%is_source()) CALL diamat_all(G_sym, eigval)
      CALL para_env%bcast(G_sym)
      CALL para_env%bcast(eigval)

      theta(:) = 0.0_dp
      theta(1:m) = eigval(:)
      C(:, :) = 0.0_dp
      C(1:m, 1:m) = G_sym(:, :)
      DEALLOCATE (G_sym, eigval)

      CALL timestop(handle)

   END SUBROUTINE solve_replicated

! **************************************************************************************************
!> \brief Orthonormal coefficient vectors of the restart basis from the candidate columns in their
!>        order of priority (modified Gram-Schmidt in the replicated m-space), at most n_cap of them
!> \param cand candidate coefficient columns, zero-padded to m rows: the current Ritz vectors, then
!>        the previous ones
!> \param m ...
!> \param n_cand ...
!> \param n_cap ...
!> \param Q m x n_new restart basis
!> \param n_new ...
! **************************************************************************************************
   SUBROUTINE restart_basis(cand, m, n_cand, n_cap, Q, n_new)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cand
      INTEGER, INTENT(IN)                                :: m, n_cand, n_cap
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: Q
      INTEGER, INTENT(OUT)                               :: n_new

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'restart_basis'
      REAL(KIND=dp), PARAMETER                           :: norm_drop = 1.0E-8_dp

      INTEGER                                            :: handle, ipass, j, k
      REAL(KIND=dp)                                      :: norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: basis

      CALL timeset(routineN, handle)

      ! never more columns than candidates: the paired solver's cap is close to m_max
      ALLOCATE (basis(m, MIN(n_cap, n_cand)))
      n_new = 0
      DO k = 1, n_cand
         IF (n_new == n_cap) EXIT
         basis(:, n_new + 1) = cand(1:m, k)
         ! a previous vector is nearly parallel to its successor, one pass leaves the basis
         ! non-orthonormal at the 1e-6 level and stalls the residuals there
         DO ipass = 1, 2
            DO j = 1, n_new
               basis(:, n_new + 1) = basis(:, n_new + 1) - &
                                     DOT_PRODUCT(basis(:, j), basis(:, n_new + 1))*basis(:, j)
            END DO
            norm = NORM2(basis(:, n_new + 1))
            IF (norm < norm_drop) EXIT
            basis(:, n_new + 1) = basis(:, n_new + 1)/norm
         END DO
         IF (norm < norm_drop) CYCLE
         n_new = n_new + 1
      END DO
      Q(:, :) = 0.0_dp
      Q(1:m, 1:n_new) = basis(:, 1:n_new)
      DEALLOCATE (basis)

      CALL timestop(handle)

   END SUBROUTINE restart_basis

! **************************************************************************************************
!> \brief out_ia,k = alpha sum_m V_ia,m C_mk + beta out_ia,k for m <= nv, k <= nc on the local rows ia
!>        On the npe x 1 grid whole columns are local, so the rotation is one local DGEMM
!> \param fm_V ...
!> \param nv ...
!> \param C ...
!> \param nc ...
!> \param fm_out ...
!> \param alpha ...
!> \param beta ...
!> \param out_col first column of fm_out written, 1 by default
! **************************************************************************************************
   SUBROUTINE subspace_rotate(fm_V, nv, C, nc, fm_out, alpha, beta, out_col)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: nv
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(IN)                                      :: C
      INTEGER, INTENT(IN)                                :: nc
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_out
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta
      INTEGER, INTENT(IN), OPTIONAL                      :: out_col

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'subspace_rotate'

      INTEGER                                            :: handle, nrow_local, o_col

      CALL timeset(routineN, handle)

      o_col = 1
      IF (PRESENT(out_col)) o_col = out_col
      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local)
      IF (nrow_local > 0) THEN
         CALL DGEMM('N', 'N', nrow_local, nc, nv, alpha, fm_V%local_data, SIZE(fm_V%local_data, 1), &
                    C, SIZE(C, 1), beta, fm_out%local_data(:, o_col:o_col + nc - 1), SIZE(fm_out%local_data, 1))
      END IF

      CALL timestop(handle)

   END SUBROUTINE subspace_rotate

! **************************************************************************************************
!> \brief V_ia,k <- sum_p V_ia,p Q_pk for p <= m, k <= nq through the work matrix
!> \param fm_V ...
!> \param m ...
!> \param Q ...
!> \param nq ...
!> \param fm_work ...
! **************************************************************************************************
   SUBROUTINE rotate_in_place(fm_V, m, Q, nq, fm_work)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Q
      INTEGER, INTENT(IN)                                :: nq
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_work

      INTEGER                                            :: nrow_local

      CALL subspace_rotate(fm_V, m, Q(1:m, 1:nq), nq, fm_work, 1.0_dp, 0.0_dp)
      CALL cp_fm_to_fm(fm_work, fm_V, nq)
      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local)
      fm_V%local_data(1:nrow_local, nq + 1:m) = 0.0_dp

   END SUBROUTINE rotate_in_place

! **************************************************************************************************
!> \brief T_ia,t0+k-1 <- T_ia,t0+k-1 + alpha U_ia,u0+k-1 for k <= n
!>        On the npe x 1 grid whole columns are local, so the update needs no communication
!> \param alpha ...
!> \param fm_U ...
!> \param u0 first column read
!> \param fm_T ...
!> \param t0 first column written
!> \param n ...
! **************************************************************************************************
   SUBROUTINE columns_axpy(alpha, fm_U, u0, fm_T, t0, n)

      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_U
      INTEGER, INTENT(IN)                                :: u0
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: t0, n

      INTEGER                                            :: nrow_local

      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      fm_T%local_data(1:nrow_local, t0:t0 + n - 1) = fm_T%local_data(1:nrow_local, t0:t0 + n - 1) + &
                                                     alpha*fm_U%local_data(1:nrow_local, u0:u0 + n - 1)

   END SUBROUTINE columns_axpy

! **************************************************************************************************
!> \brief G_kl = sum_ia U_ia,u0+k-1 V_ia,v0+l-1 for k <= nu, l <= nv, local product and sum over the ranks
!>        The npe x 1 grid splits the rows over the ranks, which is what the one allreduce completes
!> \param fm_U ...
!> \param u0 ...
!> \param nu ...
!> \param fm_V ...
!> \param v0 ...
!> \param nv ...
!> \param para_env ...
!> \param G ...
! **************************************************************************************************
   SUBROUTINE subspace_gram(fm_U, u0, nu, fm_V, v0, nv, para_env, G)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_U
      INTEGER, INTENT(IN)                                :: u0, nu
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: v0, nv
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(OUT)                                     :: G

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'subspace_gram'

      INTEGER                                            :: handle, nrow_local

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_U, nrow_local=nrow_local)
      G(:, :) = 0.0_dp
      IF (nrow_local > 0) THEN
         CALL DGEMM('T', 'N', nu, nv, nrow_local, 1.0_dp, &
                    fm_U%local_data(:, u0:u0 + nu - 1), SIZE(fm_U%local_data, 1), &
                    fm_V%local_data(:, v0:v0 + nv - 1), SIZE(fm_V%local_data, 1), &
                    0.0_dp, G, SIZE(G, 1))
      END IF
      CALL para_env%sum(G)

      CALL timestop(handle)

   END SUBROUTINE subspace_gram

! **************************************************************************************************
!> \brief fm_red(1:nu, v0:v0+nv-1) <- sum_ia U_ia,k V_ia,v0+l-1 for k <= nu, l <= nv, written in
!>        column blocks of nb so the replicated temporary is nu x nb and never nu x nv
!> \param fm_U ...
!> \param nu ...
!> \param fm_V ...
!> \param v0 first column of fm_V read and first column of fm_red written
!> \param nv ...
!> \param nb width of one column block
!> \param para_env ...
!> \param fm_red the distributed reduced matrix, on the diagonalisation grid
! **************************************************************************************************
   SUBROUTINE reduced_gram_blocks(fm_U, nu, fm_V, v0, nv, nb, para_env, fm_red)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_U
      INTEGER, INTENT(IN)                                :: nu
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: v0, nv, nb
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reduced_gram_blocks'

      INTEGER                                            :: c0, handle, nc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk

      CALL timeset(routineN, handle)

      ALLOCATE (blk(nu, nb))
      DO c0 = v0, v0 + nv - 1, nb
         nc = MIN(nb, v0 + nv - c0)
         CALL subspace_gram(fm_U, 1, nu, fm_V, c0, nc, para_env, blk(:, 1:nc))
         CALL cp_fm_set_submatrix(fm_red, blk(:, 1:nc), 1, c0, nu, nc)
      END DO
      DEALLOCATE (blk)

      CALL timestop(handle)

   END SUBROUTINE reduced_gram_blocks

! **************************************************************************************************
!> \brief fm_red(1:m+nt, m+1:m+nt) <- blk and fm_red(m+1:m+nt, 1:m) <- blk(1:m, :)^T, the two
!>        blocks a basis extension by nt columns adds to the reduced matrix
!> \param fm_red ...
!> \param blk (m+nt) x nt, the new columns' overlaps with the whole extended basis
!> \param m ...
!> \param nt ...
! **************************************************************************************************
   SUBROUTINE reduced_extend(fm_red, blk, m, nt)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: blk
      INTEGER, INTENT(IN)                                :: m, nt

      CALL cp_fm_set_submatrix(fm_red, blk, 1, m + 1, m + nt, nt)
      IF (m > 0) CALL cp_fm_set_submatrix(fm_red, blk(1:m, :), m + 1, 1, nt, m, transpose=.TRUE.)

   END SUBROUTINE reduced_extend

! **************************************************************************************************
!> \brief A new m x m matrix on blacs_env holding the live block fm_red(1:m, 1:m)
!> \param fm_red ...
!> \param m ...
!> \param para_env ...
!> \param blacs_env ...
!> \param name ...
!> \param fm_live created here, released by the caller
! **************************************************************************************************
   SUBROUTINE reduced_live_block(fm_red, m, para_env, blacs_env, name, fm_live)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      INTEGER, INTENT(IN)                                :: m
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env
      TYPE(cp_blacs_env_type), INTENT(IN), TARGET        :: blacs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_live

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=m, ncol_global=m)
      CALL cp_fm_create(fm_live, fm_struct, name=name)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm_submat(fm_red, fm_live, m, m, 1, 1, 1, 1)

   END SUBROUTINE reduced_live_block

! **************************************************************************************************
!> \brief fm_red(1:k, 1:k) <- Q^T fm_red(1:m, 1:m) Q for the m x k restart basis Q, the rest of the
!>        matrix zeroed: two distributed products, no replicated m x m
!> \param fm_red ...
!> \param m ...
!> \param Q replicated, zero-padded to m rows
!> \param k ...
!> \param para_env ...
!> \param blacs_env ...
! **************************************************************************************************
   SUBROUTINE reduced_rotate(fm_red, m, Q, k, para_env, blacs_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Q
      INTEGER, INTENT(IN)                                :: k
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env
      TYPE(cp_blacs_env_type), INTENT(IN), TARGET        :: blacs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reduced_rotate'

      INTEGER                                            :: handle
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_GQ, fm_live, fm_Q, fm_QGQ

      CALL timeset(routineN, handle)

      CALL reduced_live_block(fm_red, m, para_env, blacs_env, "bse_reduced_live", fm_live)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=m, ncol_global=k)
      CALL cp_fm_create(fm_Q, fm_struct, name="bse_restart_basis")
      CALL cp_fm_create(fm_GQ, fm_struct, name="bse_reduced_GQ")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=k, ncol_global=k)
      CALL cp_fm_create(fm_QGQ, fm_struct, name="bse_reduced_QGQ")
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_set_submatrix(fm_Q, Q(1:m, 1:k))
      CALL parallel_gemm('N', 'N', m, k, m, 1.0_dp, fm_live, fm_Q, 0.0_dp, fm_GQ)
      CALL parallel_gemm('T', 'N', k, k, m, 1.0_dp, fm_Q, fm_GQ, 0.0_dp, fm_QGQ)
      CALL cp_fm_set_all(fm_red, 0.0_dp)
      CALL cp_fm_to_fm_submat(fm_QGQ, fm_red, k, k, 1, 1, 1, 1)

      CALL cp_fm_release(fm_live)
      CALL cp_fm_release(fm_Q)
      CALL cp_fm_release(fm_GQ)
      CALL cp_fm_release(fm_QGQ)

      CALL timestop(handle)

   END SUBROUTINE reduced_rotate

! **************************************************************************************************
!> \brief max_kl |sum_ia U_ia,k V_ia,l - δ_kl| over k, l <= m, in column blocks of nb
!> \param fm_U ...
!> \param fm_V ...
!> \param m ...
!> \param nb ...
!> \param para_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION orthonormality_deviation(fm_U, fm_V, m, nb, para_env) RESULT(dev)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_U, fm_V
      INTEGER, INTENT(IN)                                :: m, nb
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp)                                      :: dev

      INTEGER                                            :: c0, l, nc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk

      ALLOCATE (blk(m, nb))
      dev = 0.0_dp
      DO c0 = 1, m, nb
         nc = MIN(nb, m - c0 + 1)
         CALL subspace_gram(fm_U, 1, m, fm_V, c0, nc, para_env, blk(:, 1:nc))
         DO l = 1, nc
            blk(c0 + l - 1, l) = blk(c0 + l - 1, l) - 1.0_dp
         END DO
         dev = MAX(dev, MAXVAL(ABS(blk(:, 1:nc))))
      END DO
      DEALLOCATE (blk)

   END FUNCTION orthonormality_deviation

! **************************************************************************************************
!> \brief Euclidean norms of nt columns of a row-distributed matrix
!>        Each rank sums the squares of its own rows of the npe x 1 grid, then one allreduce
!> \param fm_T ...
!> \param first_col ...
!> \param nt ...
!> \param para_env ...
!> \param norms ...
! **************************************************************************************************
   SUBROUTINE column_norms(fm_T, first_col, nt, para_env, norms)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col, nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: norms

      INTEGER                                            :: k, nrow_local

      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      norms(:) = 0.0_dp
      DO k = 1, nt
         norms(k) = SUM(fm_T%local_data(1:nrow_local, first_col + k - 1)**2)
      END DO
      CALL para_env%sum(norms)
      norms(:) = SQRT(norms(:))

   END SUBROUTINE column_norms

! **************************************************************************************************
!> \brief Davidson corrections t_ia,k = r_ia,k/(d_ia - θ_k) of the selected roots k, written to fm_Z
!>        Each rank divides the entries of its own rows of the npe x 1 grid, no communication
!> \param fm_R residuals, column r_offset + k belongs to root k
!> \param selected ...
!> \param theta ...
!> \param diag ...
!> \param fm_Z ...
!> \param first_col ...
!> \param r_offset columns of fm_R before the residuals, 0 by default
! **************************************************************************************************
   SUBROUTINE davidson_corrections(fm_R, selected, theta, diag, fm_Z, first_col, r_offset)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_R
      INTEGER, DIMENSION(:), INTENT(IN)                  :: selected
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: theta, diag
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(IN)                                :: first_col
      INTEGER, INTENT(IN), OPTIONAL                      :: r_offset

      CHARACTER(LEN=*), PARAMETER :: routineN = 'davidson_corrections'
      REAL(KIND=dp), PARAMETER :: eref_scale = 0.99_dp, threshold = 16.0_dp*EPSILON(1.0_dp)

      INTEGER                                            :: handle, iloc, it, k, nrow_local, r_off
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp)                                      :: denom

! near-degeneracy guard as in tddfpt_compute_residual_vects

      CALL timeset(routineN, handle)

      r_off = 0
      IF (PRESENT(r_offset)) r_off = r_offset
      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      DO it = 1, SIZE(selected)
         k = selected(it)
         DO iloc = 1, nrow_local
            denom = diag(row_indices(iloc)) - theta(k)
            IF (ABS(denom) < threshold) denom = denom + (1.0_dp - eref_scale)*theta(k)
            fm_Z%local_data(iloc, first_col + it - 1) = fm_R%local_data(iloc, r_off + k)/denom
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE davidson_corrections

! **************************************************************************************************
!> \brief T_ia,k <- T_ia,k - sum_m V_ia,m sum_jb U_jb,m T_jb,k for m <= nv and nt columns k, twice;
!>        U = V for orthonormal columns V, else the dual basis with sum_ia U_ia,m V_ia,n = δ_mn
!>        One allreduce per pass inside the Gram, the subtraction stays on the rows of the npe x 1 grid
!> \param fm_T ...
!> \param first_col ...
!> \param nt ...
!> \param fm_V ...
!> \param nv ...
!> \param para_env ...
!> \param fm_dual U
! **************************************************************************************************
   SUBROUTINE project_out(fm_T, first_col, nt, fm_V, nv, para_env, fm_dual)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col, nt
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: nv
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_dual

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'project_out'

      INTEGER                                            :: handle, ipass, nrow_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coeff

      IF (nt == 0 .OR. nv == 0) RETURN
      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      ALLOCATE (coeff(nv, nt))
      DO ipass = 1, 2
         IF (PRESENT(fm_dual)) THEN
            CALL subspace_gram(fm_dual, 1, nv, fm_T, first_col, nt, para_env, coeff)
         ELSE
            CALL subspace_gram(fm_V, 1, nv, fm_T, first_col, nt, para_env, coeff)
         END IF
         IF (nrow_local > 0) THEN
            CALL DGEMM('N', 'N', nrow_local, nt, nv, -1.0_dp, &
                       fm_V%local_data, SIZE(fm_V%local_data, 1), coeff, nv, 1.0_dp, &
                       fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1))
         END IF
      END DO
      DEALLOCATE (coeff)
      CALL timestop(handle)

   END SUBROUTINE project_out

! **************************************************************************************************
!> \brief Normalises nt columns and removes those that vanished in the projection
!>        Only the norms communicate; the scaling and the zero-fill stay on the local rows
!> \param fm_T ...
!> \param first_col ...
!> \param nt on exit the number of columns kept
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE drop_small_columns(fm_T, first_col, nt, para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col
      INTEGER, INTENT(INOUT)                             :: nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      REAL(KIND=dp), PARAMETER                           :: norm_drop = 1.0E-10_dp

      INTEGER                                            :: k, n_kept, nrow_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: norms

      IF (nt == 0) RETURN
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      ALLOCATE (norms(nt))
      CALL column_norms(fm_T, first_col, nt, para_env, norms)
      n_kept = 0
      DO k = 1, nt
         IF (norms(k) < norm_drop) CYCLE
         n_kept = n_kept + 1
         fm_T%local_data(1:nrow_local, first_col + n_kept - 1) = &
            fm_T%local_data(1:nrow_local, first_col + k - 1)/norms(k)
      END DO
      fm_T%local_data(1:nrow_local, first_col + n_kept:first_col + nt - 1) = 0.0_dp
      nt = n_kept
      DEALLOCATE (norms)

   END SUBROUTINE drop_small_columns

! **************************************************************************************************
!> \brief Orthonormalises nt columns by two Cholesky QR passes on the replicated Gram matrix. A
!>        failed factorization means a numerically dependent block: the directions of the Gram matrix
!>        below the drop tolerance are removed (canonical orthogonalisation), nt shrinks, and the
!>        passes start again. The Gram is the only communication: R is factorised on one rank and
!>        broadcast, the triangular solve runs on the rows of the npe x 1 grid
!> \param fm_T ...
!> \param first_col ...
!> \param nt ...
!> \param para_env ...
!> \param ok ...
!> \param n_dropped columns removed as dependent
! **************************************************************************************************
   SUBROUTINE cholqr2(fm_T, first_col, nt, para_env, ok, n_dropped)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col
      INTEGER, INTENT(INOUT)                             :: nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(OUT)                               :: ok
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_dropped

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cholqr2'
      REAL(KIND=dp), PARAMETER                           :: eig_drop = 1.0E-10_dp

      INTEGER                                            :: handle, info, ipass, k, n_drop_calls, &
                                                            n_keep, nrow_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eig
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G, R, T_new, U

      ok = .TRUE.
      IF (PRESENT(n_dropped)) n_dropped = 0
      IF (nt == 0) RETURN
      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      ALLOCATE (G(nt, nt), R(nt, nt))

      n_drop_calls = 0
      ipass = 0
      DO WHILE (ipass < 2)
         ipass = ipass + 1
         CALL subspace_gram(fm_T, first_col, nt, fm_T, first_col, nt, para_env, G)
         R(:, :) = G(:, :)
         info = 0
         IF (para_env%is_source()) CALL DPOTRF('U', nt, R, nt, info)
         CALL para_env%bcast(info)
         IF (info /= 0) THEN
            n_drop_calls = n_drop_calls + 1
            IF (n_drop_calls > 2) THEN
               ok = .FALSE.
               EXIT
            END IF
            ! T <- T U Λ^-1/2 over the eigenpairs of G above the drop tolerance
            ALLOCATE (eig(nt), U(nt, nt))
            CALL solve_replicated(G, nt, para_env, eig, U)
            n_keep = COUNT(eig > eig_drop*eig(nt))
            IF (n_keep == 0) THEN
               ok = .FALSE.
               DEALLOCATE (eig, U)
               EXIT
            END IF
            DO k = 1, n_keep
               U(:, nt - n_keep + k) = U(:, nt - n_keep + k)/SQRT(eig(nt - n_keep + k))
            END DO
            IF (nrow_local > 0) THEN
               ALLOCATE (T_new(nrow_local, n_keep))
               CALL DGEMM('N', 'N', nrow_local, n_keep, nt, 1.0_dp, &
                          fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1), &
                          U(:, nt - n_keep + 1:nt), nt, 0.0_dp, T_new, nrow_local)
               fm_T%local_data(1:nrow_local, first_col:first_col + n_keep - 1) = T_new(:, :)
               fm_T%local_data(1:nrow_local, first_col + n_keep:first_col + nt - 1) = 0.0_dp
               DEALLOCATE (T_new)
            END IF
            IF (PRESENT(n_dropped)) n_dropped = n_dropped + nt - n_keep
            nt = n_keep
            DEALLOCATE (eig, U, G, R)
            ALLOCATE (G(nt, nt), R(nt, nt))
            ipass = 0
            CYCLE
         END IF
         CALL para_env%bcast(R)
         IF (nrow_local > 0) THEN
            CALL DTRSM('R', 'U', 'N', 'N', nrow_local, nt, 1.0_dp, R, nt, &
                       fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1))
         END IF
      END DO

      DEALLOCATE (G, R)
      CALL timestop(handle)

   END SUBROUTINE cholqr2

! **************************************************************************************************
!> \brief Largest deviation of the nt new columns from orthonormality against the whole basis,
!>        sum_ia V_ia,m D_ia,n - δ_mn for m = 1..nv+nt, n = nv+1..nv+nt; D = V in the Euclidean case
!> \param fm_V ...
!> \param fm_dual ...
!> \param nv ...
!> \param nt ...
!> \param para_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION extension_deviation(fm_V, fm_dual, nv, nt, para_env) RESULT(dev)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V, fm_dual
      INTEGER, INTENT(IN)                                :: nv, nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp)                                      :: dev

      INTEGER                                            :: k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: S

      dev = 0.0_dp
      IF (nt == 0) RETURN
      ALLOCATE (S(nv + nt, nt))
      CALL subspace_gram(fm_V, 1, nv + nt, fm_dual, nv + 1, nt, para_env, S)
      DO k = 1, nt
         S(nv + k, k) = S(nv + k, k) - 1.0_dp
      END DO
      dev = MAXVAL(ABS(S))
      DEALLOCATE (S)

   END FUNCTION extension_deviation

! **************************************************************************************************
!> \brief Debug check of the Davidson result against the full diagonalization of the explicit A:
!>        energies, and per multiplet the smallest singular value of X_ref^T X
!> \param fm_A_explicit ...
!> \param exc_ens ...
!> \param fm_X ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE bse_davidson_refcheck_tda(fm_A_explicit, exc_ens, fm_X, unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A_explicit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_ens
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_davidson_refcheck_tda'

      INTEGER                                            :: diag_info, handle, n_ov, n_ref, n_want
      REAL(KIND=dp)                                      :: min_overlap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ref_ens
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: X_dav, X_ref
      TYPE(cp_fm_type)                                   :: fm_eigvec, fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_A_explicit, nrow_global=n_ov)
      n_want = SIZE(exc_ens)

      ! the eigensolver destroys its input
      CALL cp_fm_create(fm_work, fm_A_explicit%matrix_struct)
      CALL cp_fm_to_fm(fm_A_explicit, fm_work)
      CALL cp_fm_create(fm_eigvec, fm_A_explicit%matrix_struct)
      ALLOCATE (ref_ens(n_ov))
      CALL choose_eigv_solver(fm_work, fm_eigvec, ref_ens, diag_info)
      IF (diag_info /= 0) CPABORT("Reference diagonalization of A failed in the BSE Davidson check")
      CALL cp_fm_release(fm_work)

      n_ref = last_required_root(ref_ens, n_want, n_ov)
      ALLOCATE (X_ref(n_ov, n_ref), X_dav(n_ov, n_want))
      CALL cp_fm_get_submatrix(fm_eigvec, X_ref, 1, 1, n_ov, n_ref)
      CALL cp_fm_get_submatrix(fm_X, X_dav)
      CALL cp_fm_release(fm_eigvec)

      CALL min_multiplet_overlap(ref_ens, n_want, n_ref, X_ref, X_dav, min_overlap)
      CALL print_refcheck(MAXVAL(ABS(exc_ens(:) - ref_ens(1:n_want))), min_overlap, unit_nr)

      DEALLOCATE (ref_ens, X_ref, X_dav)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_refcheck_tda

! **************************************************************************************************
!> \brief Debug check of an ABBA Davidson result against the Hermitian reduction of the explicit A
!>        and B, C = (A-B)^1/2 (A+B) (A-B)^1/2 with C T^n = (Ω^n)^2 T^n,
!>        (X+Y)^n = (Ω^n)^-1/2 (A-B)^1/2 T^n, (X-Y)^n = (Ω^n)^1/2 (A-B)^-1/2 T^n: energies, and per
!>        multiplet the smallest singular value of X_ref^T X - Y_ref^T Y
!> \param fm_A_explicit ...
!> \param fm_B_explicit ...
!> \param exc_ens ...
!> \param fm_X ...
!> \param fm_Y ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE bse_davidson_refcheck_abba(fm_A_explicit, fm_B_explicit, exc_ens, fm_X, fm_Y, mp2_env, &
                                         unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A_explicit, fm_B_explicit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_ens
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X, fm_Y
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_davidson_refcheck_abba'

      INTEGER                                            :: diag_info, handle, k, n_ov, n_ref, n_want
      REAL(KIND=dp)                                      :: min_overlap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ref_ens
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: diff_ref, sum_ref, X_dav, X_ref, Y_dav, &
                                                            Y_ref
      TYPE(cp_fm_type)                                   :: fm_C, fm_inv_sqrt_A_minus_B, &
                                                            fm_sqrt_A_minus_B, fm_T, fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_A_explicit, nrow_global=n_ov)
      n_want = SIZE(exc_ens)

      CALL create_hermitian_form_of_ABBA(fm_A_explicit, fm_B_explicit, fm_C, fm_sqrt_A_minus_B, &
                                         fm_inv_sqrt_A_minus_B, unit_nr, mp2_env, 0.0_dp)
      CALL cp_fm_create(fm_T, fm_C%matrix_struct)
      ALLOCATE (ref_ens(n_ov))
      CALL choose_eigv_solver(fm_C, fm_T, ref_ens, diag_info)
      IF (diag_info /= 0) CPABORT("Reference diagonalization of C failed in the BSE Davidson check")
      IF (ref_ens(1) <= 0.0_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "Reference matrix C has a non-positive eigenvalue in the BSE Davidson check")
      END IF
      ref_ens(:) = SQRT(ref_ens(:))

      n_ref = last_required_root(ref_ens, n_want, n_ov)
      ALLOCATE (sum_ref(n_ov, n_ref), diff_ref(n_ov, n_ref))
      CALL cp_fm_create(fm_work, fm_C%matrix_struct)
      CALL parallel_gemm("N", "N", n_ov, n_ref, n_ov, 1.0_dp, fm_sqrt_A_minus_B, fm_T, 0.0_dp, fm_work)
      CALL cp_fm_get_submatrix(fm_work, sum_ref, 1, 1, n_ov, n_ref)
      CALL parallel_gemm("N", "N", n_ov, n_ref, n_ov, 1.0_dp, fm_inv_sqrt_A_minus_B, fm_T, 0.0_dp, fm_work)
      CALL cp_fm_get_submatrix(fm_work, diff_ref, 1, 1, n_ov, n_ref)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_T)
      CALL cp_fm_release(fm_C)
      CALL cp_fm_release(fm_sqrt_A_minus_B)
      CALL cp_fm_release(fm_inv_sqrt_A_minus_B)

      ALLOCATE (X_ref(n_ov, n_ref), Y_ref(n_ov, n_ref), X_dav(n_ov, n_want), Y_dav(n_ov, n_want))
      DO k = 1, n_ref
         sum_ref(:, k) = sum_ref(:, k)/SQRT(ref_ens(k))
         diff_ref(:, k) = diff_ref(:, k)*SQRT(ref_ens(k))
      END DO
      X_ref(:, :) = 0.5_dp*(sum_ref(:, :) + diff_ref(:, :))
      Y_ref(:, :) = 0.5_dp*(sum_ref(:, :) - diff_ref(:, :))
      CALL cp_fm_get_submatrix(fm_X, X_dav)
      CALL cp_fm_get_submatrix(fm_Y, Y_dav)

      CALL min_multiplet_overlap(ref_ens, n_want, n_ref, X_ref, X_dav, min_overlap, Y_ref, Y_dav)
      CALL print_refcheck(MAXVAL(ABS(exc_ens(:) - ref_ens(1:n_want))), min_overlap, unit_nr)

      DEALLOCATE (ref_ens, sum_ref, diff_ref, X_ref, Y_ref, X_dav, Y_dav)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_refcheck_abba

! **************************************************************************************************
!> \brief Smallest singular value over the multiplets of O_kl = sum_ia (X_ref,ia^k X_ia^l - Y_ref,ia^k Y_ia^l),
!>        k over a multiplet of the reference and l over the Davidson states inside it
!> \param ref_ens ...
!> \param n_want ...
!> \param n_ref reference states up to the end of the multiplet that contains state n_want
!> \param X_ref ...
!> \param X_dav ...
!> \param min_overlap ...
!> \param Y_ref ...
!> \param Y_dav ...
! **************************************************************************************************
   SUBROUTINE min_multiplet_overlap(ref_ens, n_want, n_ref, X_ref, X_dav, min_overlap, Y_ref, Y_dav)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ref_ens
      INTEGER, INTENT(IN)                                :: n_want, n_ref
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: X_ref, X_dav
      REAL(KIND=dp), INTENT(OUT)                         :: min_overlap
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: Y_ref, Y_dav

      INTEGER                                            :: info, k1, k2, k2_dav, lwork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_vals, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap
      REAL(KIND=dp), DIMENSION(1, 1)                     :: dummy

      min_overlap = 1.0_dp
      k1 = 1
      DO WHILE (k1 <= n_want)
         k2 = k1
         DO WHILE (k2 < n_ref)
            IF (ref_ens(k2 + 1) - ref_ens(k2) >= deg_thresh) EXIT
            k2 = k2 + 1
         END DO
         k2_dav = MIN(k2, n_want)
         ALLOCATE (overlap(k2 - k1 + 1, k2_dav - k1 + 1), sing_vals(k2_dav - k1 + 1))
         overlap(:, :) = MATMUL(TRANSPOSE(X_ref(:, k1:k2)), X_dav(:, k1:k2_dav))
         IF (PRESENT(Y_ref)) THEN
            overlap(:, :) = overlap(:, :) - MATMUL(TRANSPOSE(Y_ref(:, k1:k2)), Y_dav(:, k1:k2_dav))
         END IF
         lwork = 5*(k2 - k1 + 1) + 10
         ALLOCATE (work(lwork))
         CALL DGESVD('N', 'N', k2 - k1 + 1, k2_dav - k1 + 1, overlap, k2 - k1 + 1, sing_vals, &
                     dummy, 1, dummy, 1, work, lwork, info)
         IF (info /= 0) CPABORT("SVD failed in the BSE Davidson check")
         min_overlap = MIN(min_overlap, MINVAL(sing_vals))
         DEALLOCATE (overlap, sing_vals, work)
         k1 = k2 + 1
      END DO

   END SUBROUTINE min_multiplet_overlap

! **************************************************************************************************
!> \brief The two lines of a reference check
!> \param dev_E largest energy deviation in Hartree
!> \param min_overlap ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_refcheck(dev_E, min_overlap, unit_nr)

      REAL(KIND=dp), INTENT(IN)                          :: dev_E, min_overlap
      INTEGER, INTENT(IN)                                :: unit_nr

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A10,T13,A,T63,ES18.6)') 'BSE|DEBUG|', &
            'Max energy deviation vs full diagonalization (eV)', dev_E*evolt
         WRITE (unit_nr, '(T2,A10,T13,A,T63,F18.12)') 'BSE|DEBUG|', &
            'Min subspace overlap vs full diagonalization', min_overlap
      END IF

   END SUBROUTINE print_refcheck

END MODULE bse_davidson

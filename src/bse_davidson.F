!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block Davidson solvers for the lowest excitations of the Bethe-Salpeter equation on top
!>        of the matrix-free application of A and B: the TDA problem sum_jb A_ia,jb X_jb = Ω X_ia
!>        (bse_davidson_tda) and the full problem, either as sum_jb [(A+B)(A-B)]_ia,jb x_jb = Ω^2 x_ia
!>        (bse_davidson_abba_mk) or as the pair sum_jb (A+B)_ia,jb y_jb = Ω x_ia, sum_jb (A-B)_ia,jb x_jb = Ω y_ia
!>        (bse_davidson_abba_os). Every solver projects its problem onto an orthonormal basis, takes the
!>        lowest Ritz pairs, extends the basis by the preconditioned residuals of the roots still in need
!>        and restarts thick at the subspace ceiling; the helpers below are shared by the three drivers
!> \par History
!>      09.2026 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_davidson

   USE bibliography,                    ONLY: Bai2012,&
                                              Davidson1975,&
                                              Fukaya2014,&
                                              Olsen1988,&
                                              Stathopoulos1998,&
                                              Stratmann1998,&
                                              Vecharynski2017,&
                                              cite_reference
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
                                              bse_iter_res_cond,&
                                              bse_memcheck_abort,&
                                              bse_memcheck_clamp,&
                                              bse_memcheck_off,&
                                              bse_memcheck_warn
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

   ! the three drivers, for the memory estimate
   INTEGER, PARAMETER, PRIVATE :: driver_tda = 1, driver_mk = 2, driver_os = 3

   PUBLIC :: bse_davidson_tda, bse_davidson_abba_mk, bse_davidson_abba_os, &
             bse_davidson_refcheck

CONTAINS

! **************************************************************************************************
!> \brief Block Davidson with thick restart for the NUM_EXC_EN lowest solutions of
!>        sum_jb A_ia,jb X_jb^n = Ω^n X_ia^n. With the orthonormal basis Z_ia,m the reduced matrix is
!>        Ar_mn = sum_ia Z_ia,m (A Z)_ia,n with eigenpairs (θ_k, c_mk), X_ia^k = sum_m Z_ia,m c_mk,
!>        and the basis grows by the corrections t_ia,k = r_ia,k/(d_ia - θ_k) of the residuals
!>        r_ia,k = sum_m (A Z)_ia,m c_mk - θ_k X_ia^k, d_ia being the diagonal chosen by PRECONDITIONER.
!>        Arrays: Ar_mn (the reduced A) in fm_red_A on the process grid, c_mk in coef_ritz, θ_k c_mk in
!>        coef_ritz_theta
!> \param mv_env slabs, prefactors and transition energies of the matrix-free A
!> \param mp2_env mp2_env%bse carries the BSE_ITERAT settings and BSE_DEBUG_PRINT
!> \param unit_nr output unit, positive on the writing rank only
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coef_restart, coef_restart_cand, &
                                                            coef_ritz, coef_ritz_theta, red_block
      TYPE(cp_blacs_env_type), POINTER                   :: diag_blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_AZ, fm_red_A, fm_work, fm_Z
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      CALL cite_reference(Davidson1975)
      CALL cite_reference(Stathopoulos1998)
      CALL cite_reference(Fukaya2014)

      para_env => mv_env%para_env
      ! every local_data use below takes whole rows per rank: the npe x 1 grid of bse_matvec_create
      CPASSERT(mv_env%blacs_env%num_pe(2) == 1)
      ! the trial vectors live on an npe x 1 grid, which is a bad one to diagonalize on
      NULLIFY (diag_blacs_env)
      CALL cp_blacs_env_create(diag_blacs_env, para_env)
      n_ov = mv_env%n_ov
      CALL davidson_sizes(mp2_env, n_ov, unit_nr, n_want, n_act, block_size)

      ALLOCATE (diag(n_ov))
      CALL bse_matvec_diagonal(mv_env, mp2_env%bse%preconditioner, diag)

      CALL subspace_ceiling(mp2_env, mv_env, driver_tda, n_act, block_size, unit_nr, m_max)

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

      ! c_mk and θ_k c_mk, the restart candidates (current and previous Ritz vectors) and the restart basis Q_pk
      ALLOCATE (coef_ritz(m_max, n_act), coef_restart_cand(m_max, 2*n_act), coef_ritz_theta(m_max, n_act))
      ALLOCATE (coef_restart(m_max, 2*n_act))
      ALLOCATE (theta(m_max), theta_prev(n_act), res(n_act), dE(n_act), conv(n_act), need(n_act), &
                selected(n_act))
      coef_restart_cand(:, :) = 0.0_dp
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
      ! both guess routines return orthonormal columns by construction; one Gram matrix checks it
      IF (extension_deviation(fm_Z, fm_Z, 0, m, para_env) > ortho_tol) THEN
         CPABORT("BSE Davidson: the initial guess is not orthonormal")
      END IF
      ! (A Z)_ia,m of the guess and Ar_mn = sum_ia Z_ia,m (A Z)_ia,n
      CALL bse_matvec_apply(mv_env, fm_Z, 1, m, fm_AZ)
      n_kernel = m
      CALL reduced_gram_blocks(fm_Z, m, fm_AZ, 1, m, block_size, para_env, fm_red_A)

      CALL print_iteration_header('Block Davidson iterations within the TDA:', unit_nr)

      all_conv = .FALSE.
      DO iter = 1, mp2_env%bse%max_iter
         t_iter = m_walltime()

         ! sum_n Ar_mn c_nk = θ_k c_mk, θ ascending
         CALL solve_reduced(fm_red_A, m, n_act, para_env, diag_blacs_env, theta, coef_ritz)

         ! thick restart onto the current and the previous Ritz vectors, no kernel application; with the
         ! restart basis Q_pk in coef_restart, Z_ia,k <- sum_p Z_ia,p Q_pk, A Z alike, Ar_kl <- sum_pq Q_pk Ar_pq Q_ql
         IF (m + block_size > m_max .AND. m > n_act .AND. m < n_ov) THEN
            coef_restart_cand(1:m, 1:n_act) = coef_ritz(1:m, 1:n_act)
            n_cand = n_act
            IF (has_prev) n_cand = 2*n_act
            CALL restart_basis(coef_restart_cand, m, n_cand, 2*n_act, coef_restart, k)
            CALL rotate_in_place(fm_Z, m, coef_restart, k, fm_work)
            CALL rotate_in_place(fm_AZ, m, coef_restart, k, fm_work)
            CALL reduced_rotate(fm_red_A, m, coef_restart, k, para_env, diag_blacs_env)
            m = k
            n_restart = n_restart + 1
            CALL solve_reduced(fm_red_A, m, n_act, para_env, diag_blacs_env, theta, coef_ritz)
         END IF

         ! r_ia,k = sum_m (A Z)_ia,m c_mk - θ_k sum_m Z_ia,m c_mk for all tracked roots k
         DO k = 1, n_act
            coef_ritz_theta(1:m, k) = coef_ritz(1:m, k)*theta(k)
         END DO
         CALL subspace_rotate(fm_AZ, m, coef_ritz(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
         CALL subspace_rotate(fm_Z, m, coef_ritz_theta(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
         CALL column_norms(fm_work, 1, n_act, para_env, res)

         dE(:) = ABS(theta(1:n_act) - theta_prev(:))
         CALL convergence_flags(mp2_env, dE, res, conv)
         n_req = last_required_root(theta, n_want, n_act)
         ! a subspace that spans all transitions is exact
         IF (m == n_ov) conv(:) = .TRUE.
         CALL roots_in_need(conv, theta, res, n_req, need)
         all_conv = .NOT. ANY(need)
         ! the masked MAXVAL is -HUGE when no root is in need; the unmasked term covers that case
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
         ! no new vector: with an unchanged basis the next iteration finds dE = 0, which lets the roots that
         ! only wait for the energy criterion converge; a second empty pass means the residuals cannot improve
         IF (nt == 0 .AND. stalled) THEN
            CPABORT("BSE Davidson: correction vectors vanished before convergence")
         END IF
         stalled = nt == 0
         CALL cholqr2(fm_Z, m + 1, nt, para_env, ok, n_dropped)
         ! cholqr2 gives up after two rounds of dropping dependent directions, or when none is left
         IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")
         n_dependent = n_dependent + n_dropped
         ! post-condition of the extension, one Gram matrix; the residual norms below rely on it
         IF (extension_deviation(fm_Z, fm_Z, m, nt, para_env) > ortho_tol) THEN
            CPABORT("BSE Davidson: the basis lost orthonormality")
         END IF

         IF (nt > 0) THEN
            CALL bse_matvec_apply(mv_env, fm_Z, m + 1, nt, fm_AZ)
            n_kernel = n_kernel + nt

            ! new columns Ar_mn = sum_ia Z_ia,m (A Z)_ia,n of the symmetric reduced matrix
            ALLOCATE (red_block(m + nt, nt))
            CALL subspace_gram(fm_Z, 1, m + nt, fm_AZ, m + 1, nt, para_env, red_block)
            CALL reduced_extend(fm_red_A, red_block, m, nt)
            DEALLOCATE (red_block)
         END IF

         ! the Ritz vectors of this iteration are the second block of restart candidates
         coef_restart_cand(:, n_act + 1:2*n_act) = 0.0_dp
         coef_restart_cand(1:m, n_act + 1:2*n_act) = coef_ritz(1:m, 1:n_act)
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
      CALL subspace_rotate(fm_Z, m, coef_ritz(1:m, 1:n_want), n_want, fm_X, 1.0_dp, 0.0_dp)

      CALL cp_fm_release(fm_Z)
      CALL cp_fm_release(fm_AZ)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_red_A)
      CALL cp_blacs_env_release(diag_blacs_env)
      DEALLOCATE (coef_ritz, coef_restart_cand, coef_ritz_theta, coef_restart, theta, theta_prev, res, dE, conv, &
                  need, selected, diag)
      IF (ALLOCATED(theta_sub)) DEALLOCATE (theta_sub)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_tda

! **************************************************************************************************
!> \brief Block Davidson with thick restart for the NUM_EXC_EN lowest solutions of the full problem,
!>        written as sum_jb [(A+B)(A-B)]_ia,jb x_jb^n = (Ω^n)^2 x_ia^n with x^n = X^n - Y^n, which is
!>        symmetric in the inner product <u,v> = sum_ia,jb u_ia (A-B)_ia,jb v_jb. The basis V_ia,m is
!>        orthonormal in that inner product, W_ia,m = sum_jb (A-B)_ia,jb V_jb,m, and the reduced matrix
!>        Hr_mn = sum_ia,jb W_ia,m (A+B)_ia,jb W_jb,n has the eigenpairs ((Ω^k)^2, c_mk). Then
!>        x_ia^k = (Ω^k)^1/2 sum_m V_ia,m c_mk and y_ia^k = (X^k+Y^k)_ia = (Ω^k)^-1/2 sum_m W_ia,m c_mk
!>        obey sum_ia x_ia^k y_ia^k = 1, the residual is r_ia,k = sum_jb (A+B)_ia,jb y_jb^k - Ω^k x_ia^k
!>        and the correction t_ia,k = r_ia,k/(d_ia^2 - (Ω^k)^2), d_ia chosen by PRECONDITIONER.
!>        Arrays: Hr_mn in fm_red_M on the process grid, c_mk in coef_ritz, (Ω^k)^1/2 c_mk in coef_x,
!>        (Ω^k)^-1/2 c_mk in coef_y, Ω^k (Ω^k)^1/2 c_mk in coef_x_omega
!> \param mv_env slabs, prefactors and transition energies of the matrix-free A and B
!> \param mp2_env mp2_env%bse carries the BSE_ITERAT settings and BSE_DEBUG_PRINT
!> \param unit_nr output unit, positive on the writing rank only
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, diag, diag_sq, omega, omega_prev, &
                                                            omega_sq, res, theta, theta_sub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coef_restart, coef_restart_cand, &
                                                            coef_ritz, coef_x, coef_x_omega, &
                                                            coef_y, red_block
      TYPE(cp_blacs_env_type), POINTER                   :: diag_blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_metric, fm_metric_vec, fm_MW, &
                                                            fm_red_M, fm_scratch, fm_V, fm_W, &
                                                            fm_work
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      CALL cite_reference(Davidson1975)
      CALL cite_reference(Vecharynski2017)
      CALL cite_reference(Bai2012)
      CALL cite_reference(Stathopoulos1998)
      CALL cite_reference(Fukaya2014)

      abba_status = abba_ok
      ab_margin = 0.0_dp
      para_env => mv_env%para_env
      ! every local_data use below takes whole rows per rank: the npe x 1 grid of bse_matvec_create
      CPASSERT(mv_env%blacs_env%num_pe(2) == 1)
      ! the trial vectors live on an npe x 1 grid, which is a bad one to diagonalize on
      NULLIFY (diag_blacs_env)
      CALL cp_blacs_env_create(diag_blacs_env, para_env)
      n_ov = mv_env%n_ov
      CALL davidson_sizes(mp2_env, n_ov, unit_nr, n_want, n_act, block_size)

      ALLOCATE (diag(n_ov))
      CALL bse_matvec_diagonal(mv_env, mp2_env%bse%preconditioner, diag)

      CALL subspace_ceiling(mp2_env, mv_env, driver_mk, n_act, block_size, unit_nr, m_max)

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

      ! c_mk, the restart candidates and basis, then (Ω^k)^1/2 c_mk, (Ω^k)^-1/2 c_mk and Ω^k (Ω^k)^1/2 c_mk
      ALLOCATE (coef_ritz(m_max, n_act), coef_restart_cand(m_max, 2*n_act), coef_restart(m_max, 2*n_act))
      ALLOCATE (coef_x(m_max, n_act), coef_y(m_max, n_act), coef_x_omega(m_max, n_act))
      ALLOCATE (theta(m_max), omega(n_act), omega_prev(n_act), res(n_act), dE(n_act), conv(n_act), &
                need(n_act), selected(n_act))
      coef_restart_cand(:, :) = 0.0_dp
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
      ALLOCATE (diag_sq(n_ov))
      diag_sq(:) = diag(:)**2

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
         ! Hr_mn = sum_ia W_ia,m (M W)_ia,n on the guess basis
         CALL reduced_gram_blocks(fm_W, m, fm_MW, 1, m, block_size, para_env, fm_red_M)

         CALL print_iteration_header('Block Davidson iterations for ABBA, (A+B)(A-B) x = E^2 x:', unit_nr)

         DO iter = 1, mp2_env%bse%max_iter
            t_iter = m_walltime()

            ! sum_n Hr_mn c_nk = θ_k c_mk with θ_k = (Ω^k)^2, ascending
            CALL solve_reduced(fm_red_M, m, n_act, para_env, diag_blacs_env, theta, coef_ritz)

            ! thick restart: V, W and MW are rotated alike by Q_pk in coef_restart, which keeps
            ! sum_ia V_ia,m W_ia,n = δ_mn, and Hr_kl <- sum_pq Q_pk Hr_pq Q_ql
            IF (m + block_size > m_max .AND. m > n_act .AND. m < n_ov) THEN
               coef_restart_cand(1:m, 1:n_act) = coef_ritz(1:m, 1:n_act)
               n_cand = n_act
               IF (has_prev) n_cand = 2*n_act
               CALL restart_basis(coef_restart_cand, m, n_cand, 2*n_act, coef_restart, k)
               CALL rotate_in_place(fm_V, m, coef_restart, k, fm_work)
               CALL rotate_in_place(fm_W, m, coef_restart, k, fm_work)
               CALL rotate_in_place(fm_MW, m, coef_restart, k, fm_work)
               CALL reduced_rotate(fm_red_M, m, coef_restart, k, para_env, diag_blacs_env)
               m = k
               n_restart = n_restart + 1
               CALL solve_reduced(fm_red_M, m, n_act, para_env, diag_blacs_env, theta, coef_ritz)
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
               coef_y(1:m, k) = coef_ritz(1:m, k)/SQRT(omega(k))
               coef_x_omega(1:m, k) = coef_ritz(1:m, k)*omega(k)*SQRT(omega(k))
            END DO
            CALL subspace_rotate(fm_MW, m, coef_y(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
            CALL subspace_rotate(fm_V, m, coef_x_omega(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
            CALL column_norms(fm_work, 1, n_act, para_env, res)

            dE(:) = ABS(omega(:) - omega_prev(:))
            CALL convergence_flags(mp2_env, dE, res, conv)
            n_req = last_required_root(omega, n_want, n_act)
            ! a subspace that spans all transitions is exact
            IF (m == n_ov) conv(:) = .TRUE.
            CALL roots_in_need(conv, omega, res, n_req, need)
            all_conv = .NOT. ANY(need)
            ! the masked MAXVAL is -HUGE when no root is in need; the unmasked term covers that case
            CALL print_iteration(iter, m, COUNT(conv(1:n_want)), &
                                 MAX(MAXVAL(res(1:n_req)), MAXVAL(res, MASK=need)), &
                                 MAXVAL(dE(1:n_req)), omega(1), m_walltime() - t_iter, unit_nr)
            IF (all_conv) EXIT

            IF (iter == mp2_env%bse%max_iter) EXIT

            ! corrections t_ia,k = r_ia,k/(d_ia^2 - (Ω^k)^2) for the lowest unconverged roots
            CALL select_roots(mp2_env, need, res, MIN(block_size, m_max - m), selected, nt)
            CALL davidson_corrections(fm_work, selected(1:nt), theta, diag_sq, fm_V, m + 1)
            CALL extend_basis_mk(mv_env, fm_V, fm_W, fm_MW, fm_scratch, m, nt, n_kernel, n_dependent, abba_status, &
                                 ab_margin)
            IF (abba_status /= abba_ok) EXIT
            ! no new vector: with an unchanged basis the next iteration finds dE = 0, which lets the roots that
            ! only wait for the energy criterion converge; a second empty pass means the residuals cannot improve
            IF (nt == 0 .AND. stalled) THEN
               CPABORT("BSE Davidson: correction vectors vanished before convergence")
            END IF
            stalled = nt == 0

            IF (nt > 0) THEN
               ! new columns Hr_mn = sum_ia W_ia,m (M W)_ia,n of the symmetric reduced matrix
               ALLOCATE (red_block(m + nt, nt))
               CALL subspace_gram(fm_W, 1, m + nt, fm_MW, m + 1, nt, para_env, red_block)
               CALL reduced_extend(fm_red_M, red_block, m, nt)
               DEALLOCATE (red_block)
            END IF

            ! the Ritz vectors of this iteration are the second block of restart candidates
            coef_restart_cand(:, n_act + 1:2*n_act) = 0.0_dp
            coef_restart_cand(1:m, n_act + 1:2*n_act) = coef_ritz(1:m, 1:n_act)
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
         ! inverse of the largest eigenvalue of the metric sum_ia V_ia,m V_ia,n
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=diag_blacs_env, &
                                  nrow_global=m, ncol_global=m)
         CALL cp_fm_create(fm_metric, fm_struct, name="bse_basis_metric")
         CALL cp_fm_create(fm_metric_vec, fm_struct, name="bse_basis_metric_vectors")
         CALL cp_fm_struct_release(fm_struct)
         CALL reduced_gram_blocks(fm_V, m, fm_V, 1, m, block_size, para_env, fm_metric)
         ALLOCATE (omega_sq(m))
         omega_sq(:) = 0.0_dp
         CALL choose_eigv_solver(fm_metric, fm_metric_vec, omega_sq)
         ab_margin = 1.0_dp/omega_sq(m)
         DEALLOCATE (omega_sq)
         CALL cp_fm_release(fm_metric)
         CALL cp_fm_release(fm_metric_vec)

         IF (ALLOCATED(theta_sub)) CALL check_guess_bound(theta_sub, omega, n_req, unit_nr)
         CALL print_summary(iter, n_kernel, n_restart, n_want, n_req, n_act, n_ov, unit_nr, ab_margin, n_dependent)

         ! x_ia^n = (Ω^n)^1/2 sum_m V_ia,m c_mn, y_ia^n = (Ω^n)^-1/2 sum_m W_ia,m c_mn,
         ! X_ia^n = (y_ia^n + x_ia^n)/2, Y_ia^n = (y_ia^n - x_ia^n)/2
         ALLOCATE (exc_ens(n_want))
         exc_ens(:) = omega(1:n_want)
         DO k = 1, n_want
            coef_x(1:m, k) = coef_ritz(1:m, k)*SQRT(omega(k))
            coef_y(1:m, k) = coef_ritz(1:m, k)/SQRT(omega(k))
         END DO
         CALL bse_matvec_vector_struct(mv_env, n_want, fm_struct)
         CALL cp_fm_create(fm_X, fm_struct, name="fm_X_davidson", set_zero=.TRUE.)
         CALL cp_fm_create(fm_Y, fm_struct, name="fm_Y_davidson", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct)
         CALL subspace_rotate(fm_W, m, coef_y(1:m, 1:n_want), n_want, fm_X, 0.5_dp, 0.0_dp)
         CALL subspace_rotate(fm_V, m, coef_x(1:m, 1:n_want), n_want, fm_X, 0.5_dp, 1.0_dp)
         CALL subspace_rotate(fm_W, m, coef_y(1:m, 1:n_want), n_want, fm_Y, 0.5_dp, 0.0_dp)
         CALL subspace_rotate(fm_V, m, coef_x(1:m, 1:n_want), n_want, fm_Y, -0.5_dp, 1.0_dp)
      END IF

      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_W)
      CALL cp_fm_release(fm_MW)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_scratch)
      CALL cp_fm_release(fm_red_M)
      CALL cp_blacs_env_release(diag_blacs_env)
      DEALLOCATE (coef_ritz, coef_restart_cand, coef_x, coef_y, coef_x_omega, coef_restart, theta, omega, omega_prev, &
                  res, dE, conv, need, selected, diag, diag_sq)
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
!>        parallel to it for a small B. Arrays: Mr_mn and Kr_mn in fm_red_M and fm_red_K on the process
!>        grid, R_mk in coef_right, L_mk in coef_left
!> \param mv_env slabs, prefactors and transition energies of the matrix-free A and B
!> \param mp2_env mp2_env%bse carries the BSE_ITERAT settings and BSE_DEBUG_PRINT
!> \param unit_nr output unit, positive on the writing rank only
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coef_left, coef_restart, &
                                                            coef_restart_cand, coef_right, &
                                                            coef_work, red_block_K, red_block_M
      TYPE(cp_blacs_env_type), POINTER                   :: diag_blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_b, fm_Kb, fm_Mb, fm_red_K, fm_red_M, &
                                                            fm_scratch, fm_work
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      CALL cite_reference(Davidson1975)
      CALL cite_reference(Olsen1988)
      CALL cite_reference(Stratmann1998)
      CALL cite_reference(Bai2012)
      CALL cite_reference(Stathopoulos1998)
      CALL cite_reference(Fukaya2014)

      abba_status = abba_ok
      ab_margin = 0.0_dp
      para_env => mv_env%para_env
      ! every local_data use below takes whole rows per rank: the npe x 1 grid of bse_matvec_create
      CPASSERT(mv_env%blacs_env%num_pe(2) == 1)
      ! the trial vectors live on an npe x 1 grid, which is a bad one to diagonalize on
      NULLIFY (diag_blacs_env)
      CALL cp_blacs_env_create(diag_blacs_env, para_env)
      n_ov = mv_env%n_ov
      CALL davidson_sizes(mp2_env, n_ov, unit_nr, n_want, n_act, block_size)

      ALLOCATE (diag(n_ov))
      CALL bse_matvec_diagonal(mv_env, mp2_env%bse%preconditioner, diag)

      ! the current pairs and one block of corrections have to fit below the ceiling
      CALL subspace_ceiling(mp2_env, mv_env, driver_os, n_act, block_size, unit_nr, m_max, min_factor=4)

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

      ! R_mk, L_mk and a work block, then the restart candidates (current and previous pairs) and the basis Q_pk
      ALLOCATE (coef_right(m_max, n_act), coef_left(m_max, n_act), coef_work(m_max, n_act))
      ALLOCATE (coef_restart_cand(m_max, 4*n_act), coef_restart(m_max, 4*n_act))
      ALLOCATE (omega(n_act), omega_prev(n_act), res(n_act), res_L(n_act), dE(n_act), conv(n_act), &
                need(n_act), selected(n_act))
      coef_restart_cand(:, :) = 0.0_dp
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
            CALL solve_reduced_paired(fm_red_M, fm_red_K, m, n_act, para_env, omega, coef_right, coef_left, kr_min, &
                                      abba_status, diag_blacs_env)
            ! the smallest eigenvalue of A-B on the span of b is that of the last Kr
            ab_margin = kr_min
            IF (abba_status /= abba_ok) EXIT

            ! thick restart: b, Mb and Kb are rotated alike; the current pairs come first, so that
            ! they survive the cap that leaves room for one block of corrections; with Q_pk in
            ! coef_restart, b_ia,k <- sum_p b_ia,p Q_pk and Mr_kl <- sum_pq Q_pk Mr_pq Q_ql, Kr alike
            IF (m + 2*block_size > m_max .AND. m > 2*n_act .AND. m < n_ov) THEN
               coef_restart_cand(1:m, 1:n_act) = coef_right(1:m, :)
               coef_restart_cand(1:m, n_act + 1:2*n_act) = coef_left(1:m, :)
               n_cand = 2*n_act
               IF (has_prev) n_cand = 4*n_act
               CALL restart_basis(coef_restart_cand, m, n_cand, m_max - 2*block_size, coef_restart, k)
               CALL rotate_in_place(fm_b, m, coef_restart, k, fm_work)
               CALL rotate_in_place(fm_Mb, m, coef_restart, k, fm_work)
               CALL rotate_in_place(fm_Kb, m, coef_restart, k, fm_work)
               CALL reduced_rotate(fm_red_M, m, coef_restart, k, para_env, diag_blacs_env)
               CALL reduced_rotate(fm_red_K, m, coef_restart, k, para_env, diag_blacs_env)
               m = k
               n_restart = n_restart + 1
               CALL solve_reduced_paired(fm_red_M, fm_red_K, m, n_act, para_env, omega, coef_right, coef_left, kr_min, &
                                         abba_status, diag_blacs_env)
               ab_margin = kr_min
               IF (abba_status /= abba_ok) EXIT
            END IF

            ! rR_ia,k = sum_m (M b)_ia,m R_mk - Ω^k sum_m b_ia,m L_mk into the columns 1..n_act and
            ! rL_ia,k = sum_m (K b)_ia,m L_mk - Ω^k sum_m b_ia,m R_mk into the columns n_act+1..2 n_act
            DO k = 1, n_act
               coef_work(1:m, k) = coef_left(1:m, k)*omega(k)
            END DO
            CALL subspace_rotate(fm_Mb, m, coef_right(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
            CALL subspace_rotate(fm_b, m, coef_work(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
            DO k = 1, n_act
               coef_work(1:m, k) = coef_right(1:m, k)*omega(k)
            END DO
            CALL subspace_rotate(fm_Kb, m, coef_left(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp, n_act + 1)
            CALL subspace_rotate(fm_b, m, coef_work(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp, n_act + 1)
            CALL column_norms(fm_work, 1, n_act, para_env, res)
            CALL column_norms(fm_work, n_act + 1, n_act, para_env, res_L)
            res(:) = MAX(res(:), res_L(:))

            dE(:) = ABS(omega(:) - omega_prev(:))
            CALL convergence_flags(mp2_env, dE, res, conv)
            n_req = last_required_root(omega, n_want, n_act)
            ! a subspace that spans all transitions is exact
            IF (m == n_ov) conv(:) = .TRUE.
            CALL roots_in_need(conv, omega, res, n_req, need)
            all_conv = .NOT. ANY(need)
            ! the masked MAXVAL is -HUGE when no root is in need; the unmasked term covers that case
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
            ! no new vector: with an unchanged basis the next iteration finds dE = 0, which lets the roots that
            ! only wait for the energy criterion converge; a second empty pass means the residuals cannot improve
            IF (nt == 0 .AND. stalled) THEN
               CPABORT("BSE Davidson: correction vectors vanished before convergence")
            END IF
            stalled = nt == 0

            IF (nt > 0) THEN
               ! new columns n of Mr_mn = sum_ia b_ia,m (M b)_ia,n and Kr_mn = sum_ia b_ia,m (K b)_ia,n
               ALLOCATE (red_block_M(m + nt, nt), red_block_K(m + nt, nt))
               CALL subspace_gram(fm_b, 1, m + nt, fm_Mb, m + 1, nt, para_env, red_block_M)
               CALL subspace_gram(fm_b, 1, m + nt, fm_Kb, m + 1, nt, para_env, red_block_K)
               CALL reduced_extend(fm_red_M, red_block_M, m, nt)
               CALL reduced_extend(fm_red_K, red_block_K, m, nt)
               DEALLOCATE (red_block_M, red_block_K)
            END IF

            ! the pairs of this iteration are the second block of restart candidates
            coef_restart_cand(:, 2*n_act + 1:4*n_act) = 0.0_dp
            coef_restart_cand(1:m, 2*n_act + 1:3*n_act) = coef_right(1:m, :)
            coef_restart_cand(1:m, 3*n_act + 1:4*n_act) = coef_left(1:m, :)
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

         IF (ALLOCATED(theta_sub)) CALL check_guess_bound(theta_sub, omega, n_req, unit_nr)
         CALL print_summary(iter, n_kernel, n_restart, n_want, n_req, n_act, n_ov, unit_nr, ab_margin, n_dependent)

         ! X_ia^n = sum_m b_ia,m (R_mn + L_mn)/2, Y_ia^n = sum_m b_ia,m (R_mn - L_mn)/2
         ALLOCATE (exc_ens(n_want))
         exc_ens(:) = omega(1:n_want)
         CALL bse_matvec_vector_struct(mv_env, n_want, fm_struct)
         CALL cp_fm_create(fm_X, fm_struct, name="fm_X_davidson", set_zero=.TRUE.)
         CALL cp_fm_create(fm_Y, fm_struct, name="fm_Y_davidson", set_zero=.TRUE.)
         CALL cp_fm_struct_release(fm_struct)
         coef_work(1:m, 1:n_want) = 0.5_dp*(coef_right(1:m, 1:n_want) + coef_left(1:m, 1:n_want))
         CALL subspace_rotate(fm_b, m, coef_work(1:m, 1:n_want), n_want, fm_X, 1.0_dp, 0.0_dp)
         coef_work(1:m, 1:n_want) = 0.5_dp*(coef_right(1:m, 1:n_want) - coef_left(1:m, 1:n_want))
         CALL subspace_rotate(fm_b, m, coef_work(1:m, 1:n_want), n_want, fm_Y, 1.0_dp, 0.0_dp)
      END IF

      CALL cp_fm_release(fm_b)
      CALL cp_fm_release(fm_Mb)
      CALL cp_fm_release(fm_Kb)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_scratch)
      CALL cp_fm_release(fm_red_M)
      CALL cp_fm_release(fm_red_K)
      CALL cp_blacs_env_release(diag_blacs_env)
      DEALLOCATE (coef_right, coef_left, coef_work, coef_restart_cand, coef_restart, omega, omega_prev, res, res_L, &
                  dE, conv, need, selected, diag)
      IF (ALLOCATED(theta_sub)) DEALLOCATE (theta_sub)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_abba_os

! **************************************************************************************************
!> \brief Makes the nt columns T that follow column m of V orthonormal in the (A-B) inner product and
!>        stores their products: T_ia,k <- T_ia,k - sum_m V_ia,m sum_jb W_jb,m T_jb,k, Euclidean
!>        orthonormalisation of T, W_ia,k = sum_jb (A-B)_ia,jb T_jb,k, T_ia,k <- sum_l T_ia,l (R^-1)_lk and
!>        W alike with sum_p R_pk R_pl = sum_ia T_ia,k W_ia,l, and (M W)_ia,k = sum_jb (A+B)_ia,jb W_jb,k
!> \param mv_env the matrix-free A and B
!> \param fm_V basis; columns 1..m are kept, m+1..m+nt hold T on entry and the kept columns on exit
!> \param fm_W (A-B) V, extended in the same columns
!> \param fm_MW (A+B) W, extended in the same columns
!> \param fm_scratch scratch of at least nt columns
!> \param m columns of the basis before the extension
!> \param nt on exit the number of columns kept
!> \param n_kernel counter of the kernel applications
!> \param n_dependent counter of the correction vectors dropped as dependent
!> \param abba_status abba_indefinite if the Gram matrix of T in the (A-B) inner product is not positive definite
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
      ! cholqr2 gives up after two rounds of dropping dependent directions, or when none is left
      IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")
      n_dependent = n_dependent + n_dropped

      IF (nt > 0) THEN
         ! W_ia,k = sum_jb (A-B)_ia,jb T_jb,k: A T into W, B T into the scratch
         CALL bse_matvec_apply(mv_env, fm_V, m + 1, nt, fm_W, fm_scratch, 1)
         CALL columns_axpy(-1.0_dp, fm_scratch, 1, fm_W, m + 1, nt)
         n_kernel = n_kernel + nt

         ! T_ia,k <- sum_l T_ia,l (R^-1)_lk and W alike, with sum_p R_pk R_pl = sum_ia T_ia,k W_ia,l
         CALL cholqr2_metric(fm_V, fm_W, m + 1, nt, mv_env%para_env, abba_status, ab_margin)

         IF (abba_status == abba_ok) THEN
            ! post-condition of the extension, one Gram matrix; the residual norms below rely on it
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
!>        G_kl = sum_ia T_ia,k (KT)_ia,l = sum_p R_pk R_pl, T_ia,k <- sum_l T_ia,l (R^-1)_lk, KT alike. A failed
!>        factorization of a G with a negative eigenvalue reports K as indefinite, any other failure
!>        is repeated once with a diagonal shift and a third pass. The Gram is the only communication:
!>        R is factorised on one rank and broadcast, the triangular solve runs on the local rows
!> \param fm_T the block to orthonormalise, columns first_col .. first_col+nt-1
!> \param fm_KT K T in the same columns, transformed alike
!> \param first_col first column of the block
!> \param nt number of columns in the block
!> \param para_env communicator of the row-distributed matrices
!> \param abba_status abba_indefinite if gram has a negative eigenvalue, else untouched
!> \param eig_min smallest eigenvalue of gram, set if it is negative
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: chol_factor, gram, gram_raw

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local, nrow_global=n_ov)
      ALLOCATE (gram_raw(nt, nt), gram(nt, nt), chol_factor(nt, nt), eigval(nt))

      n_pass = 2
      ipass = 0
      DO WHILE (ipass < n_pass)
         ipass = ipass + 1
         ! G_kl = (sum_ia T_ia,k (KT)_ia,l + sum_ia T_ia,l (KT)_ia,k)/2 = sum_p R_pk R_pl
         CALL subspace_gram(fm_T, first_col, nt, fm_KT, first_col, nt, para_env, gram_raw)
         gram(:, :) = 0.5_dp*(gram_raw(:, :) + TRANSPOSE(gram_raw(:, :)))
         chol_factor(:, :) = gram(:, :)
         info = 0
         IF (para_env%is_source()) CALL DPOTRF('U', nt, chol_factor, nt, info)
         CALL para_env%bcast(info)
         IF (info /= 0) THEN
            IF (n_pass == 3) THEN
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson: orthonormalisation in the (A-B) inner product broke down")
            END IF
            chol_factor(:, :) = gram(:, :)
            eigval(:) = 0.0_dp
            IF (para_env%is_source()) CALL diamat_all(chol_factor, eigval)
            CALL para_env%bcast(eigval)
            IF (eigval(1) < 0.0_dp) THEN
               abba_status = abba_indefinite
               eig_min = eigval(1)
               EXIT
            END IF
            ! a positive G that failed numerically: factorise G_kl + shift δ_kl and pass once more, the
            ! shift of shifted Cholesky QR (Fukaya et al. 2020), 11 (N_ov nt + nt (nt+1)) EPSILON |G|_F
            n_pass = 3
            shift = 11.0_dp*(REAL(n_ov, dp)*REAL(nt, dp) + REAL(nt, dp)*REAL(nt + 1, dp))* &
                    EPSILON(1.0_dp)*NORM2(gram)
            chol_factor(:, :) = gram(:, :)
            DO k = 1, nt
               chol_factor(k, k) = chol_factor(k, k) + shift
            END DO
            IF (para_env%is_source()) CALL DPOTRF('U', nt, chol_factor, nt, info)
            CALL para_env%bcast(info)
            IF (info /= 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson: orthonormalisation in the (A-B) inner product broke down")
            END IF
         END IF
         ! T_ia,k <- sum_l T_ia,l (R^-1)_lk and KT alike, on the local rows
         CALL para_env%bcast(chol_factor)
         IF (nrow_local > 0) THEN
            CALL DTRSM('R', 'U', 'N', 'N', nrow_local, nt, 1.0_dp, chol_factor, nt, &
                       fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1))
            CALL DTRSM('R', 'U', 'N', 'N', nrow_local, nt, 1.0_dp, chol_factor, nt, &
                       fm_KT%local_data(:, first_col:first_col + nt - 1), SIZE(fm_KT%local_data, 1))
         END IF
      END DO

      DEALLOCATE (gram_raw, gram, chol_factor, eigval)
      CALL timestop(handle)

   END SUBROUTINE cholqr2_metric

! **************************************************************************************************
!> \brief Paired reduced problem on the orthonormal basis b_ia,m,
!>        Mr_mn = sum_ia b_ia,m (M b)_ia,n and Kr_mn = sum_ia b_ia,m (K b)_ia,n,
!>        sum_n Mr_mn R_nk = Ω_k L_mk and sum_n Kr_mn L_nk = Ω_k R_mk,
!>        solved as in initial_guess_subblock through the symmetric form:
!>        Kr_mn = sum_p U_mp κ_p U_np and S_mn = sum_p U_mp √κ_p U_np, the square root of Kr,
!>        sum_nop S_mn Mr_no S_op T_pk = λ_k T_mk with Ω_k = √λ_k,
!>        R_mk = sum_n S_mn T_nk/√Ω_k and L_mk = sum_n Mr_mn R_nk/Ω_k.
!>        κ_1 <= 0 reports K as indefinite on the span of b.
!>        Every m x m step runs on the grid of blacs_env: the two eigensolves through
!>        choose_eigv_solver and the four products through parallel_gemm; only the pair
!>        coefficients come back replicated, R_mk in coef_right and L_mk in coef_left
!> \param fm_red_M Mr in its live block, on the grid of blacs_env
!> \param fm_red_K Kr in its live block, on the grid of blacs_env
!> \param m dimension of the live block
!> \param n_act number of pairs formed
!> \param para_env communicator of the grid
!> \param omega Ω_k in Hartree
!> \param coef_right coefficients of y_k in column k, zero below row m
!> \param coef_left coefficients of x_k in column k, zero below row m
!> \param kr_min κ_1, the smallest eigenvalue of Kr
!> \param abba_status abba_indefinite when κ_1 <= 0; omega, coef_right and coef_left are then untouched
!> \param blacs_env grid the reduced problem is solved on
! **************************************************************************************************
   SUBROUTINE solve_reduced_paired(fm_red_M, fm_red_K, m, n_act, para_env, omega, coef_right, coef_left, kr_min, &
                                   abba_status, blacs_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red_M, fm_red_K
      INTEGER, INTENT(IN)                                :: m, n_act
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: omega
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: coef_right, coef_left
      REAL(KIND=dp), INTENT(OUT)                         :: kr_min
      INTEGER, INTENT(INOUT)                             :: abba_status
      TYPE(cp_blacs_env_type), INTENT(IN), TARGET        :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_reduced_paired'

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: kappa, lambda
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_H, fm_H_eigvec, fm_K_eigvec, &
                                                            fm_K_sqrt, fm_L, fm_Mr, fm_R, fm_work

      CALL timeset(routineN, handle)
      ! n_act pairs are extracted from an m-dimensional basis
      CPASSERT(n_act <= m)

      ! Kr_mn = sum_p U_mp κ_p U_np, S_mn = sum_p U_mp √κ_p U_np, H_mp = sum_no S_mn Mr_no S_op, T its eigenvectors:
      ! κ, U in fm_K_eigvec, S in fm_K_sqrt, the symmetrised Mr, H, T in fm_H_eigvec and a work matrix,
      ! all m x m on the grid of blacs_env
      ALLOCATE (kappa(m), lambda(m))
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=m, ncol_global=m)
      CALL cp_fm_create(fm_K_eigvec, fm_struct, name="bse_paired_U")
      CALL cp_fm_create(fm_K_sqrt, fm_struct, name="bse_paired_S")
      CALL cp_fm_create(fm_Mr, fm_struct, name="bse_paired_Mr")
      CALL cp_fm_create(fm_H, fm_struct, name="bse_paired_H")
      CALL cp_fm_create(fm_H_eigvec, fm_struct, name="bse_paired_T")
      CALL cp_fm_create(fm_work, fm_struct, name="bse_paired_work")
      CALL cp_fm_struct_release(fm_struct)

      ! Kr_mn = sum_p U_mp κ_p U_np on the symmetrised live block, (Kr_mn + Kr_nm)/2
      CALL cp_fm_to_fm_submat(fm_red_K, fm_work, m, m, 1, 1, 1, 1)
      CALL cp_fm_transpose(fm_work, fm_K_eigvec)
      CALL cp_fm_scale_and_add(0.5_dp, fm_work, 0.5_dp, fm_K_eigvec)
      kappa(:) = 0.0_dp
      CALL choose_eigv_solver(fm_work, fm_K_eigvec, kappa)
      kr_min = kappa(1)

      IF (kappa(1) <= 0.0_dp) THEN
         abba_status = abba_indefinite
      ELSE
         ! S_mn = sum_p U_mp √κ_p U_np
         CALL cp_fm_to_fm(fm_K_eigvec, fm_work)
         CALL cp_fm_column_scale(fm_work, SQRT(kappa))
         CALL parallel_gemm('N', 'T', m, m, m, 1.0_dp, fm_work, fm_K_eigvec, 0.0_dp, fm_K_sqrt)

         ! H_mp = sum_no S_mn Mr_no S_op with Mr_no <- (Mr_no + Mr_on)/2, then T holds its eigenvectors
         CALL cp_fm_to_fm_submat(fm_red_M, fm_Mr, m, m, 1, 1, 1, 1)
         CALL cp_fm_transpose(fm_Mr, fm_H_eigvec)
         CALL cp_fm_scale_and_add(0.5_dp, fm_Mr, 0.5_dp, fm_H_eigvec)
         CALL parallel_gemm('N', 'N', m, m, m, 1.0_dp, fm_Mr, fm_K_sqrt, 0.0_dp, fm_work)
         CALL parallel_gemm('N', 'N', m, m, m, 1.0_dp, fm_K_sqrt, fm_work, 0.0_dp, fm_H)
         lambda(:) = 0.0_dp
         CALL choose_eigv_solver(fm_H, fm_H_eigvec, lambda)
         ! H is congruent to Mr = sum_ia b (M b), so λ_1 <= 0 means (A+B) indefinite on the span of b, or
         ! round-off in a nearly dependent basis, which extension_deviation excludes at every extension
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
         CALL parallel_gemm('N', 'N', m, n_act, m, 1.0_dp, fm_K_sqrt, fm_H_eigvec, 0.0_dp, fm_R)
         CALL cp_fm_column_scale(fm_R, 1.0_dp/SQRT(omega(1:n_act)))
         CALL parallel_gemm('N', 'N', m, n_act, m, 1.0_dp, fm_Mr, fm_R, 0.0_dp, fm_L)
         CALL cp_fm_column_scale(fm_L, 1.0_dp/omega(1:n_act))

         coef_right(:, :) = 0.0_dp
         coef_left(:, :) = 0.0_dp
         CALL cp_fm_get_submatrix(fm_R, coef_right(1:m, 1:n_act))
         CALL cp_fm_get_submatrix(fm_L, coef_left(1:m, 1:n_act))
         CALL cp_fm_release(fm_R)
         CALL cp_fm_release(fm_L)
      END IF

      CALL cp_fm_release(fm_K_eigvec)
      CALL cp_fm_release(fm_K_sqrt)
      CALL cp_fm_release(fm_Mr)
      CALL cp_fm_release(fm_H)
      CALL cp_fm_release(fm_H_eigvec)
      CALL cp_fm_release(fm_work)
      DEALLOCATE (kappa, lambda)
      CALL timestop(handle)

   END SUBROUTINE solve_reduced_paired

! **************************************************************************************************
!> \brief Makes the nt columns T that follow column m of b orthonormal to b and to each other,
!>        T_ia,k <- T_ia,k - sum_m b_ia,m sum_jb b_jb,m T_jb,k, then T_ia,k <- sum_l T_ia,l (R^-1)_lk with
!>        sum_p R_pk R_pl = sum_ia T_ia,k T_ia,l (cholqr2), and stores the products
!>        (M b)_ia,k = sum_jb (A+B)_ia,jb b_jb,k, (K b)_ia,k = sum_jb (A-B)_ia,jb b_jb,k of the kept
!>        columns from one application of A and B
!> \param mv_env the matrix-free A and B
!> \param fm_b basis; columns 1..m are kept, m+1..m+nt hold T on entry and the kept columns on exit
!> \param fm_Mb (A+B) b, extended in the same columns
!> \param fm_Kb (A-B) b, extended in the same columns
!> \param fm_scratch scratch of at least nt columns
!> \param m columns of the basis before the extension
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
      ! cholqr2 gives up after two rounds of dropping dependent directions, or when none is left
      IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")
      n_dependent = n_dependent + n_dropped

      IF (nt > 0) THEN
         ! post-condition of the extension, one Gram matrix; the residual norms below rely on it
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
!> \brief Convergence flag of every tracked root k from CONVERGENCE_CRITERION: EN dE_k < EPS_ENERGY,
!>        RES |r_k| < EPS_RES, OR either of the two, AND both
!> \param mp2_env CONVERGENCE_CRITERION, EPS_ENERGY and EPS_RES
!> \param dE energy changes in Hartree
!> \param res residual norms
!> \param conv one flag per tracked root
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
!> \param energies tracked energies in Hartree, ascending
!> \param n_want requested states, NUM_EXC_EN
!> \param n_act tracked states, the upper limit
!> \return the last root that has to converge, n_want <= n_req <= n_act
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
!>        lies within |r_k| of E_k, so such a root can still belong to the requested states; the ABBA
!>        solvers, whose residuals are not those of a symmetric problem, use the same rule as a heuristic
!> \param conv convergence flags of the tracked roots
!> \param energies tracked energies in Hartree, ascending
!> \param res residual norms
!> \param n_req last required root
!> \param need one flag per tracked root, set here
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
!> \param mp2_env CONVERGENCE_CRITERION and EPS_RES
!> \param need roots in need, from roots_in_need
!> \param res residual norms
!> \param n_max largest number of roots
!> \param selected indices of the selected roots in ascending order, zero beyond nt
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
!> \param title line printed above the column heads
!> \param unit_nr output unit, positive on the writing rank only
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
!> \param iter iteration count
!> \param m subspace dimension
!> \param n_conv converged roots among the requested ones
!> \param max_res largest residual norm over the required roots and the roots in need
!> \param max_dE in Hartree
!> \param lowest_E in Hartree
!> \param t_iter wall time of the iteration in seconds
!> \param unit_nr output unit, positive on the writing rank only
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
!> \param res residual norms
!> \param conv convergence flags
!> \param need roots in need
!> \param n_req last required root
!> \param unit_nr output unit, positive on the writing rank only
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
      ! unconverged roots would enter postprocess_bse as converged states; the table above is the diagnostic
      CPABORT("BSE Davidson: MAX_ITER reached before convergence")

   END SUBROUTINE abort_unconverged

! **************************************************************************************************
!> \brief Debug table of all tracked roots of a converged solve
!> \param energies in Hartree
!> \param res residual norms
!> \param n_want requested states
!> \param n_act tracked states
!> \param unit_nr output unit, positive on the writing rank only
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
!> \param iter iterations run
!> \param n_kernel kernel applications
!> \param n_restart thick restarts
!> \param n_want requested states
!> \param n_req last required root
!> \param n_act tracked states
!> \param n_ov number of transitions
!> \param unit_nr output unit, positive on the writing rank only
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

      IF (n_req > n_want .AND. n_req == n_act .AND. n_act < n_ov .AND. unit_nr > 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "BSE Davidson: the degenerate group at NUM_EXC_EN reaches the last buffer state, "// &
                      "further members may lie above it. Increase BSE_ITERAT%NUM_BUFFER_STATES.")
      END IF

   END SUBROUTINE print_summary

! **************************************************************************************************
!> \brief Number of wanted and tracked states and the block size, clamped to the problem size
!> \param mp2_env NUM_EXC_EN, NUM_BUFFER_STATES, BLOCK_SIZE and MAX_ITER
!> \param n_ov number of transitions
!> \param unit_nr output unit, positive on the writing rank only
!> \param n_want NUM_EXC_EN, at most n_ov
!> \param n_act wanted plus buffer states
!> \param block_size BLOCK_SIZE, or MIN(32, n_act) for -1; at most n_act
! **************************************************************************************************
   SUBROUTINE davidson_sizes(mp2_env, n_ov, unit_nr, n_want, n_act, block_size)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: n_ov, unit_nr
      INTEGER, INTENT(OUT)                               :: n_want, n_act, block_size

      IF (mp2_env%bse%num_exc_en < 1) CPABORT("BSE_ITERAT%NUM_EXC_EN must be at least 1")
      IF (mp2_env%bse%num_buffer_states < 0) CPABORT("BSE_ITERAT%NUM_BUFFER_STATES must not be negative")
      IF (mp2_env%bse%block_size < 1 .AND. mp2_env%bse%block_size /= -1) THEN
         CPABORT("BSE_ITERAT%BLOCK_SIZE must be at least 1, or -1 for the default")
      END IF
      IF (mp2_env%bse%max_iter < 1) CPABORT("BSE_ITERAT%MAX_ITER must be at least 1")

      n_want = mp2_env%bse%num_exc_en
      IF (n_want > n_ov) THEN
         IF (unit_nr > 0) CALL cp_warn(__LOCATION__, &
                                       "BSE_ITERAT%NUM_EXC_EN exceeds the number of transitions and is reduced to it.")
         n_want = n_ov
      END IF
      n_act = MIN(n_want + mp2_env%bse%num_buffer_states, n_ov)
      block_size = mp2_env%bse%block_size
      ! 32 was within 5 % of the best width on the L8 deck at 200 states; wider passes cost restarts
      IF (block_size == -1) block_size = MIN(32, n_act)
      IF (block_size > n_act) THEN
         IF (unit_nr > 0) CALL cp_warn(__LOCATION__, &
                                       "BSE_ITERAT%BLOCK_SIZE exceeds the number of tracked states and is reduced to it.")
         block_size = n_act
      END IF

   END SUBROUTINE davidson_sizes

! **************************************************************************************************
!> \brief Largest subspace dimension m_max: MAX_SUBSPACE_FACTOR x n_act, or 20 x n_act for -1, with
!>        the solver's arrays at m_max checked against the memory budget per rank as MEMORY_CHECK says;
!>        CLAMP lowers m_max to the largest that fits, never below min_factor x n_act
!> \param mp2_env MAX_SUBSPACE_FACTOR, MEMORY_CHECK and MEMORY_BUDGET_GB
!> \param mv_env n_ov, the RI share per rank, block_cols and mem_fraction
!> \param driver driver_tda, driver_mk or driver_os
!> \param n_act tracked states
!> \param block_size correction vectors per iteration
!> \param unit_nr output unit, positive on the writing rank only
!> \param m_max the ceiling, min_factor x n_act <= m_max <= n_ov
!> \param min_factor smallest subspace in vectors per tracked state, 3 by default
! **************************************************************************************************
   SUBROUTINE subspace_ceiling(mp2_env, mv_env, driver, n_act, block_size, unit_nr, m_max, min_factor)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      INTEGER, INTENT(IN)                                :: driver, n_act, block_size, unit_nr
      INTEGER, INTENT(OUT)                               :: m_max
      INTEGER, INTENT(IN), OPTIONAL                      :: min_factor

      CHARACTER(LEN=16)                                  :: bud_str, fac_str, mem_str
      INTEGER                                            :: m_floor, m_hi, m_lo, m_mid, m_start, &
                                                            min_fac, n_ri_max, nb
      LOGICAL                                            :: from_free, over, skipped
      REAL(KIND=dp)                                      :: budget_GB, dist_GB, mem_avail_GB, repl_GB

      min_fac = 3
      IF (PRESENT(min_factor)) min_fac = min_factor
      WRITE (fac_str, '(I16)') min_fac

      ! the start value: the keyword, or for -1 the generous 20 x n_act that the memory check below trims to what fits
      IF (mp2_env%bse%max_subspace_factor == -1) THEN
         m_start = 20*n_act
      ELSE IF (mp2_env%bse%max_subspace_factor >= min_fac) THEN
         m_start = mp2_env%bse%max_subspace_factor*n_act
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "BSE_ITERAT%MAX_SUBSPACE_FACTOR must be at least "//TRIM(ADJUSTL(fac_str))// &
                       " (or -1) for the chosen solver")
      END IF
      m_floor = MIN(min_fac*n_act, mv_env%n_ov)
      m_start = MIN(MAX(m_start, m_floor), mv_env%n_ov)
      m_max = m_start

      ! the widest kernel call of the driver caps the pass width; the largest RI share of any rank
      ! keeps the estimate, and so m_max, identical on all ranks
      SELECT CASE (driver)
      CASE (driver_tda)
         nb = n_act
      CASE (driver_mk)
         nb = block_size
      CASE DEFAULT
         nb = 2*block_size
      END SELECT
      IF (mv_env%block_cols > 0) nb = MIN(nb, mv_env%block_cols)
      n_ri_max = mv_env%n_ri_loc
      CALL mv_env%para_env%max(n_ri_max)

      ! the budget: as given, or the share mem_fraction of the free memory per rank
      from_free = mp2_env%bse%memory_budget_gb < 0.0_dp
      skipped = .FALSE.
      budget_GB = mp2_env%bse%memory_budget_gb
      IF (from_free) THEN
         CALL mp_mem_avail_per_rank_GB(mv_env%para_env, mem_avail_GB)
         budget_GB = mv_env%mem_fraction*mem_avail_GB
         skipped = mem_avail_GB <= 0.0_dp
         IF (skipped .AND. mp2_env%bse%memory_check /= bse_memcheck_off .AND. unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, "BSE Davidson: free memory not detectable, memory check skipped")
         END IF
      END IF

      CALL davidson_footprint(driver, mv_env, m_max, n_act, block_size, nb, n_ri_max, dist_GB, repl_GB)
      over = .NOT. skipped .AND. mp2_env%bse%memory_check /= bse_memcheck_off .AND. &
             dist_GB + repl_GB > budget_GB

      IF (over .AND. mp2_env%bse%memory_check == bse_memcheck_clamp) THEN
         CALL davidson_footprint(driver, mv_env, m_floor, n_act, block_size, nb, n_ri_max, dist_GB, repl_GB)
         IF (dist_GB + repl_GB > budget_GB) THEN
            WRITE (mem_str, '(F16.3)') dist_GB + repl_GB
            WRITE (bud_str, '(F16.3)') budget_GB
            CALL cp_abort(__LOCATION__, &
                          "BSE Davidson: not enough memory for the smallest subspace of "//TRIM(ADJUSTL(fac_str))// &
                          " vectors per state: "//TRIM(ADJUSTL(mem_str))//" GB per MPI rank, the budget is "// &
                          TRIM(ADJUSTL(bud_str))//" GB. Raise MEMORY_BUDGET_GB or use more MPI ranks.")
         END IF
         ! the footprint grows with m: bisect for the largest m within the budget
         m_lo = m_floor
         m_hi = m_start
         DO WHILE (m_hi - m_lo > 1)
            m_mid = (m_lo + m_hi)/2
            CALL davidson_footprint(driver, mv_env, m_mid, n_act, block_size, nb, n_ri_max, dist_GB, repl_GB)
            IF (dist_GB + repl_GB <= budget_GB) THEN
               m_lo = m_mid
            ELSE
               m_hi = m_mid
            END IF
         END DO
         m_max = m_lo
         CALL davidson_footprint(driver, mv_env, m_max, n_act, block_size, nb, n_ri_max, dist_GB, repl_GB)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (mp2_env%bse%max_subspace_factor == -1) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Maximum subspace dimension (20 x buffered states)', &
               m_start
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Maximum subspace dimension (MAX_SUBSPACE_FACTOR)', &
               m_start
         END IF
         IF (m_max < m_start) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Maximum subspace dimension after the memory check', &
               m_max
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory per MPI rank at the subspace ceiling (GB)'
         WRITE (unit_nr, '(T2,A4,T9,A,T67,F14.3)') 'BSE|', 'Distributed arrays', dist_GB
         WRITE (unit_nr, '(T2,A4,T9,A,T67,F14.3)') 'BSE|', 'Replicated arrays', repl_GB
         WRITE (unit_nr, '(T2,A4,T9,A,T67,F14.3)') 'BSE|', 'Total', dist_GB + repl_GB
         IF (.NOT. skipped) THEN
            IF (from_free) THEN
               WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Fraction of the free memory made available', &
                  mv_env%mem_fraction
            END IF
            WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Memory budget per MPI rank (GB)', budget_GB
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (mp2_env%bse%memory_check == bse_memcheck_off) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: OFF'
         ELSE IF (skipped) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: skipped, free memory not detectable'
         ELSE
            SELECT CASE (mp2_env%bse%memory_check)
            CASE (bse_memcheck_warn)
               IF (over) THEN
                  WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: WARN, the total exceeds the budget'
               ELSE
                  WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: within the budget, no WARN raised'
               END IF
            CASE (bse_memcheck_clamp)
               IF (over) THEN
                  WRITE (unit_nr, '(T2,A4,T7,A,I0,A,I0)') 'BSE|', &
                     'Memory check: CLAMP applied, subspace ceiling reduced from ', m_start, ' to ', m_max
               ELSE
                  WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: within the budget, no CLAMP applied'
               END IF
            CASE (bse_memcheck_abort)
               IF (over) THEN
                  WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: ABORT, the total exceeds the budget'
               ELSE
                  WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Memory check: within the budget, no ABORT raised'
               END IF
            END SELECT
         END IF
      END IF

      IF (over) THEN
         WRITE (mem_str, '(F16.3)') dist_GB + repl_GB
         WRITE (bud_str, '(F16.3)') budget_GB
         SELECT CASE (mp2_env%bse%memory_check)
         CASE (bse_memcheck_warn)
            IF (unit_nr > 0) CALL cp_warn(__LOCATION__, &
                                          "BSE Davidson: the solver's arrays at the subspace ceiling need "// &
                                          TRIM(ADJUSTL(mem_str))//" GB per MPI rank, above the budget of "// &
                                          TRIM(ADJUSTL(bud_str))//" GB.")
         CASE (bse_memcheck_clamp)
            IF (unit_nr > 0) CALL cp_warn(__LOCATION__, &
                                          "BSE Davidson: subspace ceiling reduced to fit the memory budget; raise "// &
                                          "MEMORY_BUDGET_GB or lower MAX_SUBSPACE_FACTOR to silence this.")
         CASE (bse_memcheck_abort)
            CALL cp_abort(__LOCATION__, &
                          "BSE Davidson: the solver's arrays at the subspace ceiling need "//TRIM(ADJUSTL(mem_str))// &
                          " GB per MPI rank, the budget is "//TRIM(ADJUSTL(bud_str))//" GB. Raise MEMORY_BUDGET_GB, "// &
                          "lower MAX_SUBSPACE_FACTOR, or set MEMORY_CHECK CLAMP.")
         END SELECT
      END IF

   END SUBROUTINE subspace_ceiling

! **************************************************************************************************
!> \brief Memory per MPI rank of the solver's arrays at subspace dimension m, split into the family
!>        spread over the ranks and the family every rank holds whole. The counts are those of the
!>        ALLOCATE and cp_fm_create statements of the drivers at a thick restart, their peak:
!>        distributed 8 B (N_ov (n_basis m + n_work + n_scratch) + (n_red + n_solve) m^2)/P,
!>        replicated 8 B (c_coef m n_act + m block_size + n_mv N_ov nb + n_ri nb + n_virt virt^2)
!> \param driver driver_tda, driver_mk or driver_os
!> \param mv_env n_ov, virt and the communicator
!> \param m subspace dimension
!> \param n_act tracked states
!> \param block_size correction vectors per iteration
!> \param nb kernel pass width
!> \param n_ri RI functions of a rank, the largest over the ranks
!> \param dist_GB distributed family, falls as 1/P
!> \param repl_GB replicated family, flat in P
! **************************************************************************************************
   SUBROUTINE davidson_footprint(driver, mv_env, m, n_act, block_size, nb, n_ri, dist_GB, repl_GB)

      INTEGER, INTENT(IN)                                :: driver
      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      INTEGER, INTENT(IN)                                :: m, n_act, block_size, nb, n_ri
      REAL(KIND=dp), INTENT(OUT)                         :: dist_GB, repl_GB

      INTEGER                                            :: c_coef, n_basis, n_mv, n_red, n_scratch, &
                                                            n_solve, n_virt, n_work
      REAL(KIND=dp)                                      :: m_real, n_ov_real

      ! distributed 8 B (N_ov (n_basis m + n_work + n_scratch) + (n_red + n_solve) m^2)/P: the trial vector
      ! arrays, the work and scratch columns, the reduced and the per-solve m x m matrices;
      ! replicated 8 B (c_coef m n_act + m block_size + n_mv N_ov nb + n_ri nb + n_virt virt^2): the m x n_act
      ! coefficient blocks with the restart work array, the pass buffers of the kernel, virt x virt for B
      SELECT CASE (driver)
      CASE (driver_tda)
         n_basis = 2
         n_work = MIN(2*n_act, m)
         n_scratch = 0
         n_red = 1
         n_solve = 2
         c_coef = 8
         n_mv = 3
         n_virt = 0
      CASE (driver_mk)
         n_basis = 3
         n_work = MIN(2*n_act, m)
         n_scratch = block_size
         n_red = 1
         n_solve = 2
         c_coef = 10
         n_mv = 4
         n_virt = 1
      CASE (driver_os)
         n_basis = 3
         n_work = MAX(2*n_act, MIN(4*n_act, m))
         n_scratch = 2*block_size
         n_red = 2
         n_solve = 6
         c_coef = 15
         n_mv = 4
         n_virt = 1
      CASE DEFAULT
         ! unreachable for the three driver constants; the DEFAULT keeps every counter defined
         CPABORT("BSE Davidson: unknown driver in the memory estimate")
      END SELECT

      ! in REAL: m^2 overflows INTEGER(4) above 46340
      m_real = REAL(m, dp)
      n_ov_real = REAL(mv_env%n_ov, dp)
      dist_GB = 8.0E-9_dp*(n_ov_real*REAL(n_basis*m + n_work + n_scratch, dp) + REAL(n_red + n_solve, dp)*m_real*m_real)/ &
                REAL(mv_env%para_env%num_pe, dp)
      repl_GB = 8.0E-9_dp*(REAL(c_coef, dp)*m_real*REAL(n_act, dp) + m_real*REAL(block_size, dp) + &
                           REAL(n_mv, dp)*n_ov_real*REAL(nb, dp) + REAL(n_ri, dp)*REAL(nb, dp) + &
                           REAL(n_virt, dp)*REAL(mv_env%virt, dp)**2)

   END SUBROUTINE davidson_footprint

! **************************************************************************************************
!> \brief Unit vectors Z_ia,k = δ_ia,diag_order(k) on the entries diag_order(1..n_used) with the lowest diagonal, never
!>        cutting a degenerate group. The sorted diagonal is replicated; each rank writes the entries
!>        that fall on its own rows
!> \param diag the diagonal chosen by PRECONDITIONER, replicated
!> \param n_guess guess vectors wanted
!> \param n_max largest number of guess vectors
!> \param deg_thresh_diag entries of the sorted diagonal closer than this form one degenerate group
!> \param fm_Z receives the unit vectors in its columns 1..n_used
!> \param n_used guess vectors written, n_guess <= n_used <= n_max
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
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: diag_order
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_sorted

      n_ov = SIZE(diag)
      ALLOCATE (diag_sorted(n_ov), diag_order(n_ov))
      diag_sorted(:) = diag(:)
      CALL sort(diag_sorted, n_ov, diag_order)

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
            IF (row_indices(iloc) == diag_order(k)) fm_Z%local_data(iloc, k) = 1.0_dp
         END DO
      END DO

      DEALLOCATE (diag_sorted, diag_order)

   END SUBROUTINE initial_guess

! **************************************************************************************************
!> \brief Initial guess from the exact A (and B) on the n_sub transitions with the lowest diagonal:
!>        TDA: sum_q A_sub,pq V_qk = θ_k V_pk; ABBA: x_pk = sum_q (K_sub^-1/2)_pq t_qk with
!>        sum_qrs (K_sub^1/2)_pq M_sub,qr (K_sub^1/2)_rs t_sk = θ_k^2 t_pk, K = A-B, M = A+B. The Ritz
!>        values θ_k are upper bounds of the true energies (Rayleigh-Ritz for TDA, Cauchy-type
!>        interlacing of Bai and Li for the linear response problem), kept for check_guess_bound. A
!>        degenerate group at the upper boundary is taken as a whole. The block and its eigenvectors
!>        are replicated; each rank writes the entries that fall on its own rows
!> \param mv_env the matrix-free A and B, for the exact block
!> \param mp2_env NUM_GUESS_TRANSITIONS and BSE_DEBUG_PRINT
!> \param diag the diagonal chosen by PRECONDITIONER, replicated; its lowest entries select the block
!> \param n_guess guess vectors wanted
!> \param n_max largest number of guess vectors
!> \param deg_thresh_diag Ritz values closer than this form one degenerate group
!> \param fm_Z receives the guess vectors in its columns 1..n_used
!> \param n_used guess vectors written, n_guess <= n_used <= n_max; 0 for abba_indefinite
!> \param theta_sub allocated here, the n_sub Ritz values in Hartree
!> \param unit_nr output unit, positive on the writing rank only
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
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: diag_order, rank_in_block
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      LOGICAL                                            :: ok
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_sorted, eig_K
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: A_sub, B_sub, eigvec_H, eigvec_K, &
                                                            guess_x, guess_y, H_sub, K_mhalf, &
                                                            K_phalf, K_sub, M_sub
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      IF (PRESENT(n_pair)) n_pair = 0
      para_env => mv_env%para_env
      n_ov = SIZE(diag)
      n_sub = MIN(n_ov, MAX(mp2_env%bse%num_guess_transitions, n_guess))
      ALLOCATE (diag_sorted(n_ov), diag_order(n_ov), rank_in_block(n_ov))
      diag_sorted(:) = diag(:)
      CALL sort(diag_sorted, n_ov, diag_order)

      ALLOCATE (A_sub(n_sub, n_sub), guess_x(n_sub, n_sub), theta_sub(n_sub))
      ok = .TRUE.
      IF (PRESENT(abba_status)) THEN
         abba_status = abba_ok
         ALLOCATE (B_sub(n_sub, n_sub), K_sub(n_sub, n_sub), M_sub(n_sub, n_sub), H_sub(n_sub, n_sub), &
                   eigvec_K(n_sub, n_sub), eigvec_H(n_sub, n_sub), K_phalf(n_sub, n_sub), &
                   K_mhalf(n_sub, n_sub), eig_K(n_sub))
         CALL bse_matvec_subblock(mv_env, diag_order(1:n_sub), A_sub, B_sub)
         ! (K_sub^±1/2)_pq = sum_r U_pr κ_r^±1/2 U_qr from K_sub,pq = sum_r U_pr κ_r U_qr
         K_sub(:, :) = A_sub(:, :) - B_sub(:, :)
         CALL solve_replicated(K_sub, n_sub, para_env, eig_K, eigvec_K)
         ab_margin = eig_K(1)
         IF (eig_K(1) <= 0.0_dp) THEN
            abba_status = abba_indefinite
            ok = .FALSE.
            n_used = 0
         ELSE
            DO k = 1, n_sub
               K_phalf(:, k) = eigvec_K(:, k)*SQRT(eig_K(k))
               K_mhalf(:, k) = eigvec_K(:, k)/SQRT(eig_K(k))
            END DO
            K_phalf(:, :) = MATMUL(K_phalf, TRANSPOSE(eigvec_K))
            K_mhalf(:, :) = MATMUL(K_mhalf, TRANSPOSE(eigvec_K))
            ! M_sub = 2 A_sub - K_sub, H_ps = sum_qr (K^1/2)_pq M_qr (K^1/2)_rs, x_pk = sum_q (K^-1/2)_pq t_qk
            M_sub(:, :) = 2.0_dp*A_sub(:, :) - K_sub(:, :)
            H_sub(:, :) = MATMUL(K_phalf, MATMUL(M_sub, K_phalf))
            CALL solve_replicated(H_sub, n_sub, para_env, theta_sub, eigvec_H)
            ! the same failure as in solve_reduced_paired: with K_sub positive definite, M_sub is indefinite
            IF (theta_sub(1) <= 0.0_dp) THEN
               CPABORT("BSE Davidson: the guess block gives a non-positive squared excitation energy")
            END IF
            theta_sub(:) = SQRT(theta_sub(:))
            guess_x(:, :) = MATMUL(K_mhalf, eigvec_H)
            IF (PRESENT(n_pair)) THEN
               ! partners guess_y_k = K_sub guess_x_k/θ_k
               ALLOCATE (guess_y(n_sub, n_sub))
               guess_y(:, :) = MATMUL(K_sub, guess_x)
               DO k = 1, n_sub
                  guess_y(:, k) = guess_y(:, k)/theta_sub(k)
               END DO
            END IF
         END IF
         DEALLOCATE (B_sub, K_sub, M_sub, H_sub, eigvec_K, eigvec_H, K_phalf, K_mhalf, eig_K)
      ELSE
         CALL bse_matvec_subblock(mv_env, diag_order(1:n_sub), A_sub)
         CALL solve_replicated(A_sub, n_sub, para_env, theta_sub, guess_x)
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

         rank_in_block(:) = 0
         DO k = 1, n_sub
            rank_in_block(diag_order(k)) = k
         END DO
         CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
         DO k = k_first, n_used
            DO iloc = 1, nrow_local
               IF (rank_in_block(row_indices(iloc)) > 0) fm_Z%local_data(iloc, k) = guess_x(rank_in_block(row_indices(iloc)), k)
            END DO
         END DO
         IF (PRESENT(n_pair)) THEN
            n_pair = MIN(n_used - k_first + 1, n_max - n_used)
            DO k = 1, n_pair
               DO iloc = 1, nrow_local
                  IF (rank_in_block(row_indices(iloc)) > 0) THEN
                     fm_Z%local_data(iloc, n_used + k) = guess_y(rank_in_block(row_indices(iloc)), k_first + k - 1)
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

      DEALLOCATE (diag_sorted, diag_order, rank_in_block, A_sub, guess_x)
      IF (ALLOCATED(guess_y)) DEALLOCATE (guess_y)

      CALL timestop(handle)

   END SUBROUTINE initial_guess_subblock

! **************************************************************************************************
!> \brief A converged energy above the Ritz value of the guess block with the same index is a
!>        provable miss: both are upper bounds of the true eigenvalue with that index
!> \param theta_sub Ritz values of the guess block in Hartree, ascending
!> \param energies converged energies in Hartree
!> \param n_req roots compared, the required ones
!> \param unit_nr output unit, positive on the writing rank only
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
      IF (excess > deg_thresh .AND. unit_nr > 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "A converged BSE state lies above the upper bound given by the guess block, "// &
                      "so a lower state was missed. Raise BSE_ITERAT%NUM_GUESS_TRANSITIONS or "// &
                      "BSE_ITERAT%NUM_BUFFER_STATES.")
      END IF

   END SUBROUTINE check_guess_bound

! **************************************************************************************************
!> \brief sum_q (red_pq + red_qp)/2 c_qk = θ_k c_pk for p, q <= m, red_pq the live block of fm_red, solved
!>        on the grid of blacs_env: θ_k in theta ascending, c_pk of the first n_vec k gathered into
!>        coef_ritz
!> \param fm_red the reduced matrix; its live block is the first m rows and columns
!> \param m dimension of the live block
!> \param n_vec eigenvectors gathered; the drivers never read more than n_act
!> \param para_env communicator of the grid
!> \param blacs_env grid the reduced problem is solved on
!> \param theta θ_k in the first m entries, zero beyond
!> \param coef_ritz c_pk in the first m rows and n_vec columns, zero beyond
! **************************************************************************************************
   SUBROUTINE solve_reduced(fm_red, m, n_vec, para_env, blacs_env, theta, coef_ritz)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      INTEGER, INTENT(IN)                                :: m, n_vec
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env
      TYPE(cp_blacs_env_type), INTENT(IN), TARGET        :: blacs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: theta
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: coef_ritz

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_reduced'

      INTEGER                                            :: handle, n_get
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      TYPE(cp_fm_type)                                   :: fm_reduced_eigvec, fm_reduced_sym

      CALL timeset(routineN, handle)

      CALL reduced_live_block(fm_red, m, para_env, blacs_env, "bse_reduced", fm_reduced_sym)
      CALL cp_fm_create(fm_reduced_eigvec, fm_reduced_sym%matrix_struct, name="bse_reduced_vectors")
      ! a_pq <- (red_pq + red_qp)/2: the eigensolver reads one triangle
      CALL cp_fm_transpose(fm_reduced_sym, fm_reduced_eigvec)
      CALL cp_fm_scale_and_add(0.5_dp, fm_reduced_sym, 0.5_dp, fm_reduced_eigvec)
      ALLOCATE (eigval(m))
      eigval(:) = 0.0_dp
      CALL choose_eigv_solver(fm_reduced_sym, fm_reduced_eigvec, eigval)

      theta(:) = 0.0_dp
      theta(1:m) = eigval(:)
      coef_ritz(:, :) = 0.0_dp
      n_get = MIN(m, n_vec)
      CALL cp_fm_get_submatrix(fm_reduced_eigvec, coef_ritz(1:m, 1:n_get), 1, 1, m, n_get)

      DEALLOCATE (eigval)
      CALL cp_fm_release(fm_reduced_sym)
      CALL cp_fm_release(fm_reduced_eigvec)

      CALL timestop(handle)

   END SUBROUTINE solve_reduced

! **************************************************************************************************
!> \brief sum_q (mat_pq + mat_qp)/2 coef_qk = θ_k coef_pk for p, q <= m, a small replicated block solved
!>        on one rank and broadcast: the fixed-size blocks of the guess and of cholqr2, never the
!>        reduced matrix
!> \param mat the block; its first m rows and columns are read
!> \param m dimension of the block
!> \param para_env communicator; the solve runs on its source rank
!> \param theta ascending
!> \param coef all m eigenvectors in the columns
! **************************************************************************************************
   SUBROUTINE solve_replicated(mat, m, para_env, theta, coef)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: mat
      INTEGER, INTENT(IN)                                :: m
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: theta
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: coef

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_replicated'

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_sym

      CALL timeset(routineN, handle)

      ALLOCATE (mat_sym(m, m), eigval(m))
      mat_sym(:, :) = 0.5_dp*(mat(1:m, 1:m) + TRANSPOSE(mat(1:m, 1:m)))
      eigval(:) = 0.0_dp
      IF (para_env%is_source()) CALL diamat_all(mat_sym, eigval)
      CALL para_env%bcast(mat_sym)
      CALL para_env%bcast(eigval)

      theta(:) = 0.0_dp
      theta(1:m) = eigval(:)
      coef(:, :) = 0.0_dp
      coef(1:m, 1:m) = mat_sym(:, :)
      DEALLOCATE (mat_sym, eigval)

      CALL timestop(handle)

   END SUBROUTINE solve_replicated

! **************************************************************************************************
!> \brief Orthonormal coefficient vectors of the restart basis from the candidate columns in their
!>        order of priority, at most n_cap of them: per candidate k, u_p = cand_pk -
!>        sum_j (sum_q coef_restart_qj cand_qk) coef_restart_pj over the vectors j kept so far, twice,
!>        then coef_restart_p,new = u_p/|u| unless |u| is below coef_norm_drop
!> \param cand candidate coefficient columns, zero-padded to m rows: the current Ritz vectors, then
!>        the previous ones
!> \param m rows of the coefficient vectors, the current basis size
!> \param n_cand candidate columns
!> \param n_cap largest number of vectors kept
!> \param coef_restart m x n_new restart basis
!> \param n_new vectors kept, at most MIN(n_cap, n_cand)
! **************************************************************************************************
   SUBROUTINE restart_basis(cand, m, n_cand, n_cap, coef_restart, n_new)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cand
      INTEGER, INTENT(IN)                                :: m, n_cand, n_cap
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: coef_restart
      INTEGER, INTENT(OUT)                               :: n_new

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'restart_basis'
      REAL(KIND=dp), PARAMETER                           :: coef_norm_drop = 1.0E-8_dp

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
            IF (norm < coef_norm_drop) EXIT
            basis(:, n_new + 1) = basis(:, n_new + 1)/norm
         END DO
         ! relative: a candidate that keeps less than this of its unit norm repeats the vectors kept
         IF (norm < coef_norm_drop) CYCLE
         n_new = n_new + 1
      END DO
      coef_restart(:, :) = 0.0_dp
      coef_restart(1:m, 1:n_new) = basis(:, 1:n_new)
      DEALLOCATE (basis)

      CALL timestop(handle)

   END SUBROUTINE restart_basis

! **************************************************************************************************
!> \brief out_ia,k = alpha sum_m V_ia,m coef_mk + beta out_ia,k for m <= nv, k <= nc on the local rows ia
!>        On the npe x 1 grid whole columns are local, so the rotation is one local DGEMM
!> \param fm_V basis, columns 1..nv read
!> \param nv columns of fm_V read
!> \param coef nv x nc coefficients, replicated
!> \param nc columns written
!> \param fm_out receives the rotated columns
!> \param alpha scale of the product
!> \param beta scale of the previous content of fm_out
!> \param out_col first column of fm_out written, 1 by default
! **************************************************************************************************
   SUBROUTINE subspace_rotate(fm_V, nv, coef, nc, fm_out, alpha, beta, out_col)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: nv
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(IN)                                      :: coef
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
                    coef, SIZE(coef, 1), beta, fm_out%local_data(:, o_col:o_col + nc - 1), SIZE(fm_out%local_data, 1))
      END IF

      CALL timestop(handle)

   END SUBROUTINE subspace_rotate

! **************************************************************************************************
!> \brief V_ia,k <- sum_p V_ia,p coef_pk for p <= m, k <= nq through the work matrix
!> \param fm_V basis; on exit columns 1..nq hold the rotation and nq+1..m are zero
!> \param m columns of fm_V read
!> \param coef restart basis, m x nq
!> \param nq columns of the rotated basis
!> \param fm_work at least nq columns
! **************************************************************************************************
   SUBROUTINE rotate_in_place(fm_V, m, coef, nq, fm_work)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: coef
      INTEGER, INTENT(IN)                                :: nq
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_work

      INTEGER                                            :: nrow_local

      CALL subspace_rotate(fm_V, m, coef(1:m, 1:nq), nq, fm_work, 1.0_dp, 0.0_dp)
      CALL cp_fm_to_fm(fm_work, fm_V, nq)
      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local)
      fm_V%local_data(1:nrow_local, nq + 1:m) = 0.0_dp

   END SUBROUTINE rotate_in_place

! **************************************************************************************************
!> \brief T_ia,t0+k-1 <- T_ia,t0+k-1 + alpha U_ia,u0+k-1 for k <= n
!>        On the npe x 1 grid whole columns are local, so the update needs no communication
!> \param alpha scale of the added columns
!> \param fm_U columns added
!> \param u0 first column read
!> \param fm_T columns updated
!> \param t0 first column written
!> \param n number of columns
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
!> \param fm_U left factor
!> \param u0 first column of fm_U
!> \param nu columns of fm_U
!> \param fm_V right factor
!> \param v0 first column of fm_V
!> \param nv columns of fm_V
!> \param para_env communicator of the row-distributed matrices
!> \param gram nu x nv, replicated
! **************************************************************************************************
   SUBROUTINE subspace_gram(fm_U, u0, nu, fm_V, v0, nv, para_env, gram)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_U
      INTEGER, INTENT(IN)                                :: u0, nu
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: v0, nv
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(OUT)                                     :: gram

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'subspace_gram'

      INTEGER                                            :: handle, nrow_local

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_U, nrow_local=nrow_local)
      gram(:, :) = 0.0_dp
      IF (nrow_local > 0) THEN
         CALL DGEMM('T', 'N', nu, nv, nrow_local, 1.0_dp, &
                    fm_U%local_data(:, u0:u0 + nu - 1), SIZE(fm_U%local_data, 1), &
                    fm_V%local_data(:, v0:v0 + nv - 1), SIZE(fm_V%local_data, 1), &
                    0.0_dp, gram, SIZE(gram, 1))
      END IF
      CALL para_env%sum(gram)

      CALL timestop(handle)

   END SUBROUTINE subspace_gram

! **************************************************************************************************
!> \brief Overlaps of the whole basis U with a block of nv columns of V, written into the distributed
!>        reduced matrix: fm_red_k,v0+l-1 <- sum_ia U_ia,k V_ia,v0+l-1 for k <= nu, l <= nv, in column
!>        blocks of nb so that the replicated temporary is nu x nb and never nu x nv
!> \param fm_U the basis, columns 1..nu
!> \param nu columns of fm_U, rows of fm_red written
!> \param fm_V the block of columns
!> \param v0 first column of fm_V read and first column of fm_red written
!> \param nv columns of fm_V read
!> \param nb width of one column block
!> \param para_env communicator of the row-distributed matrices
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gram_block

      CALL timeset(routineN, handle)

      ALLOCATE (gram_block(nu, nb))
      DO c0 = v0, v0 + nv - 1, nb
         nc = MIN(nb, v0 + nv - c0)
         CALL subspace_gram(fm_U, 1, nu, fm_V, c0, nc, para_env, gram_block(:, 1:nc))
         CALL cp_fm_set_submatrix(fm_red, gram_block(:, 1:nc), 1, c0, nu, nc)
      END DO
      DEALLOCATE (gram_block)

      CALL timestop(handle)

   END SUBROUTINE reduced_gram_blocks

! **************************************************************************************************
!> \brief fm_red_p,m+l <- red_block_pl for p <= m+nt and fm_red_m+l,p <- red_block_pl for p <= m,
!>        l <= nt: the two blocks a basis extension by nt columns adds to the symmetric reduced matrix
!> \param fm_red the reduced matrix, at least m+nt rows and columns
!> \param red_block (m+nt) x nt, the new columns' overlaps with the whole extended basis
!> \param m live dimension before the extension
!> \param nt columns added
! **************************************************************************************************
   SUBROUTINE reduced_extend(fm_red, red_block, m, nt)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: red_block
      INTEGER, INTENT(IN)                                :: m, nt

      CALL cp_fm_set_submatrix(fm_red, red_block, 1, m + 1, m + nt, nt)
      IF (m > 0) CALL cp_fm_set_submatrix(fm_red, red_block(1:m, :), m + 1, 1, nt, m, transpose=.TRUE.)

   END SUBROUTINE reduced_extend

! **************************************************************************************************
!> \brief A new m x m matrix on blacs_env holding the live block fm_red(1:m, 1:m)
!> \param fm_red the reduced matrix
!> \param m dimension of the live block
!> \param para_env communicator of the grid
!> \param blacs_env grid of the new matrix
!> \param name name of the new matrix
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
!> \brief fm_red_ij <- sum_pq coef_restart_pi fm_red_pq coef_restart_qj for p, q <= m and i, j <= k,
!>        the rest of the matrix zeroed: two distributed products, no replicated m x m
!> \param fm_red the reduced matrix, rotated in place
!> \param m live dimension before the restart
!> \param coef_restart replicated, zero-padded to m rows
!> \param k vectors of the restart basis, the live dimension after
!> \param para_env communicator of the grid
!> \param blacs_env grid of fm_red and of the products
! **************************************************************************************************
   SUBROUTINE reduced_rotate(fm_red, m, coef_restart, k, para_env, blacs_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_red
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: coef_restart
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

      CALL cp_fm_set_submatrix(fm_Q, coef_restart(1:m, 1:k))
      ! GQ_pj = sum_q fm_red_pq coef_restart_qj, then QGQ_ij = sum_p coef_restart_pi GQ_pj
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
!> \param fm_U left basis
!> \param fm_V right basis, fm_U for the Euclidean case
!> \param m columns compared
!> \param nb width of one column block
!> \param para_env communicator of the row-distributed matrices
!> \return the largest deviation
! **************************************************************************************************
   FUNCTION orthonormality_deviation(fm_U, fm_V, m, nb, para_env) RESULT(dev)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_U, fm_V
      INTEGER, INTENT(IN)                                :: m, nb
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp)                                      :: dev

      INTEGER                                            :: c0, l, nc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gram_block

      ALLOCATE (gram_block(m, nb))
      dev = 0.0_dp
      DO c0 = 1, m, nb
         nc = MIN(nb, m - c0 + 1)
         CALL subspace_gram(fm_U, 1, m, fm_V, c0, nc, para_env, gram_block(:, 1:nc))
         DO l = 1, nc
            gram_block(c0 + l - 1, l) = gram_block(c0 + l - 1, l) - 1.0_dp
         END DO
         dev = MAX(dev, MAXVAL(ABS(gram_block(:, 1:nc))))
      END DO
      DEALLOCATE (gram_block)

   END FUNCTION orthonormality_deviation

! **************************************************************************************************
!> \brief Euclidean norms of nt columns of a row-distributed matrix
!>        Each rank sums the squares of its own rows of the npe x 1 grid, then one allreduce
!> \param fm_T row-distributed matrix
!> \param first_col first column measured
!> \param nt number of columns
!> \param para_env communicator of the rows
!> \param norms nt Euclidean norms
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
!> \param selected indices of the roots corrected
!> \param theta Ritz values, indexed by root
!> \param diag d_ia, replicated
!> \param fm_Z receives the corrections
!> \param first_col first column of fm_Z written
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

      CALL timeset(routineN, handle)

      r_off = 0
      IF (PRESENT(r_offset)) r_off = r_offset
      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      DO it = 1, SIZE(selected)
         k = selected(it)
         DO iloc = 1, nrow_local
            denom = diag(row_indices(iloc)) - theta(k)
            ! near-degeneracy guard as in tddfpt_compute_residual_vects: threshold detects a numerically
            ! zero d_ia - θ_k, which the shift (1 - eref_scale) θ_k replaces by the denominator at 0.99 θ_k
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
!> \param fm_T columns projected
!> \param first_col first column of fm_T projected
!> \param nt number of columns projected
!> \param fm_V the basis projected out
!> \param nv columns of fm_V
!> \param para_env communicator of the row-distributed matrices
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
!> \param fm_T row-distributed matrix; the kept columns move to the front of the block
!> \param first_col first column of the block
!> \param nt on exit the number of columns kept
!> \param para_env communicator of the rows
! **************************************************************************************************
   SUBROUTINE drop_small_columns(fm_T, first_col, nt, para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col
      INTEGER, INTENT(INOUT)                             :: nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      REAL(KIND=dp), PARAMETER                           :: col_norm_drop = 1.0E-10_dp

      INTEGER                                            :: k, n_kept, nrow_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: norms

      IF (nt == 0) RETURN
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      ALLOCATE (norms(nt))
      CALL column_norms(fm_T, first_col, nt, para_env, norms)
      n_kept = 0
      DO k = 1, nt
         ! absolute: the corrections r/(d - θ) arrive unnormalised, a column below this is numerically zero
         IF (norms(k) < col_norm_drop) CYCLE
         n_kept = n_kept + 1
         fm_T%local_data(1:nrow_local, first_col + n_kept - 1) = &
            fm_T%local_data(1:nrow_local, first_col + k - 1)/norms(k)
      END DO
      fm_T%local_data(1:nrow_local, first_col + n_kept:first_col + nt - 1) = 0.0_dp
      nt = n_kept
      DEALLOCATE (norms)

   END SUBROUTINE drop_small_columns

! **************************************************************************************************
!> \brief Orthonormalises nt columns by two Cholesky QR passes on the replicated Gram matrix,
!>        T_ia,k <- sum_l T_ia,l (R^-1)_lk with sum_p R_pk R_pl = G_kl = sum_ia T_ia,k T_ia,l. A
!>        failed factorization means a numerically dependent block: the directions of the Gram matrix
!>        below the drop tolerance are removed (canonical orthogonalisation), nt shrinks, and the
!>        passes start again. The Gram is the only communication: R is factorised on one rank and
!>        broadcast, the triangular solve runs on the rows of the npe x 1 grid
!> \param fm_T the block to orthonormalise, columns first_col .. first_col+nt-1
!> \param first_col first column of the block
!> \param nt columns of the block; on exit the number kept
!> \param para_env communicator of the row-distributed matrix
!> \param ok .FALSE. when two drop rounds still leave no factorisable Gram matrix, or no direction survives
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
      REAL(KIND=dp), PARAMETER                           :: eig_drop_rel = 1.0E-10_dp

      INTEGER                                            :: handle, info, ipass, k, n_drop_calls, &
                                                            n_keep, nrow_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eig
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: chol_factor, gram, gram_eigvec, T_reduced

      ok = .TRUE.
      IF (PRESENT(n_dropped)) n_dropped = 0
      IF (nt == 0) RETURN
      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      ALLOCATE (gram(nt, nt), chol_factor(nt, nt))

      n_drop_calls = 0
      ipass = 0
      DO WHILE (ipass < 2)
         ipass = ipass + 1
         ! G_kl = sum_ia T_ia,k T_ia,l = sum_p R_pk R_pl
         CALL subspace_gram(fm_T, first_col, nt, fm_T, first_col, nt, para_env, gram)
         chol_factor(:, :) = gram(:, :)
         info = 0
         IF (para_env%is_source()) CALL DPOTRF('U', nt, chol_factor, nt, info)
         CALL para_env%bcast(info)
         IF (info /= 0) THEN
            n_drop_calls = n_drop_calls + 1
            IF (n_drop_calls > 2) THEN
               ok = .FALSE.
               EXIT
            END IF
            ! T_ia,k <- sum_l T_ia,l U_lk Λ_k^-1/2 over the eigenpairs of G above the drop tolerance,
            ! eig_drop_rel of the largest eigenvalue: directions below it are numerically dependent
            ALLOCATE (eig(nt), gram_eigvec(nt, nt))
            CALL solve_replicated(gram, nt, para_env, eig, gram_eigvec)
            n_keep = COUNT(eig > eig_drop_rel*eig(nt))
            IF (n_keep == 0) THEN
               ok = .FALSE.
               DEALLOCATE (eig, gram_eigvec)
               EXIT
            END IF
            DO k = 1, n_keep
               gram_eigvec(:, nt - n_keep + k) = gram_eigvec(:, nt - n_keep + k)/SQRT(eig(nt - n_keep + k))
            END DO
            IF (nrow_local > 0) THEN
               ALLOCATE (T_reduced(nrow_local, n_keep))
               CALL DGEMM('N', 'N', nrow_local, n_keep, nt, 1.0_dp, &
                          fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1), &
                          gram_eigvec(:, nt - n_keep + 1:nt), nt, 0.0_dp, T_reduced, nrow_local)
               fm_T%local_data(1:nrow_local, first_col:first_col + n_keep - 1) = T_reduced(:, :)
               fm_T%local_data(1:nrow_local, first_col + n_keep:first_col + nt - 1) = 0.0_dp
               DEALLOCATE (T_reduced)
            END IF
            IF (PRESENT(n_dropped)) n_dropped = n_dropped + nt - n_keep
            nt = n_keep
            DEALLOCATE (eig, gram_eigvec, gram, chol_factor)
            ALLOCATE (gram(nt, nt), chol_factor(nt, nt))
            ipass = 0
            CYCLE
         END IF
         ! T_ia,k <- sum_l T_ia,l (R^-1)_lk on the local rows
         CALL para_env%bcast(chol_factor)
         IF (nrow_local > 0) THEN
            CALL DTRSM('R', 'U', 'N', 'N', nrow_local, nt, 1.0_dp, chol_factor, nt, &
                       fm_T%local_data(:, first_col:first_col + nt - 1), SIZE(fm_T%local_data, 1))
         END IF
      END DO

      DEALLOCATE (gram, chol_factor)
      CALL timestop(handle)

   END SUBROUTINE cholqr2

! **************************************************************************************************
!> \brief Largest deviation of the nt new columns from orthonormality against the whole basis,
!>        sum_ia V_ia,m D_ia,n - δ_mn for m = 1..nv+nt, n = nv+1..nv+nt; D = V in the Euclidean case
!> \param fm_V the basis with its nt new columns
!> \param fm_dual D, fm_V itself in the Euclidean case
!> \param nv columns of the basis before the extension
!> \param nt new columns
!> \param para_env communicator of the row-distributed matrices
!> \return the largest deviation
! **************************************************************************************************
   FUNCTION extension_deviation(fm_V, fm_dual, nv, nt, para_env) RESULT(dev)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V, fm_dual
      INTEGER, INTENT(IN)                                :: nv, nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp)                                      :: dev

      INTEGER                                            :: k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gram

      dev = 0.0_dp
      IF (nt == 0) RETURN
      ALLOCATE (gram(nv + nt, nt))
      CALL subspace_gram(fm_V, 1, nv + nt, fm_dual, nv + 1, nt, para_env, gram)
      DO k = 1, nt
         gram(nv + k, k) = gram(nv + k, k) - 1.0_dp
      END DO
      dev = MAXVAL(ABS(gram))
      DEALLOCATE (gram)

   END FUNCTION extension_deviation

! **************************************************************************************************
!> \brief Debug check of a Davidson result against the full diagonalization of the explicit
!>        matrices: energies, and per multiplet the smallest singular value of X_ref^T X for the
!>        TDA, or of X_ref^T X - Y_ref^T Y for ABBA. ABBA diagonalizes the Hermitian reduction
!>        C = (A-B)^1/2 (A+B) (A-B)^1/2 with C T^n = (Ω^n)^2 T^n,
!>        (X+Y)^n = (Ω^n)^-1/2 (A-B)^1/2 T^n, (X-Y)^n = (Ω^n)^1/2 (A-B)^-1/2 T^n.
!>        A failed reference solve aborts: this is debug output with nothing to fall back on
!> \param fm_A_explicit A from create_A_and_B, N_ov x N_ov
!> \param exc_ens Davidson energies in Hartree
!> \param fm_X Davidson X_ia^n, column n
!> \param mp2_env handed to create_hermitian_form_of_ABBA
!> \param unit_nr output unit, positive on the writing rank only
!> \param fm_B_explicit present for an ABBA result, together with fm_Y
!> \param fm_Y Davidson Y_ia^n, with fm_B_explicit
! **************************************************************************************************
   SUBROUTINE bse_davidson_refcheck(fm_A_explicit, exc_ens, fm_X, mp2_env, unit_nr, fm_B_explicit, fm_Y)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A_explicit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_ens
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_B_explicit, fm_Y

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_davidson_refcheck'

      INTEGER                                            :: diag_info, handle, k, n_ov, n_ref, n_want
      REAL(KIND=dp)                                      :: min_overlap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ref_ens
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: diff_ref, sum_ref, X_dav, X_ref, Y_dav, &
                                                            Y_ref
      TYPE(cp_fm_type)                                   :: fm_C, fm_eigvec, fm_inv_sqrt_A_minus_B, &
                                                            fm_sqrt_A_minus_B, fm_work

      CALL timeset(routineN, handle)

      ! B and Y arrive together or not at all
      CPASSERT(PRESENT(fm_B_explicit) .EQV. PRESENT(fm_Y))
      CALL cp_fm_get_info(fm_A_explicit, nrow_global=n_ov)
      n_want = SIZE(exc_ens)
      ALLOCATE (ref_ens(n_ov))

      IF (PRESENT(fm_B_explicit)) THEN
         CALL create_hermitian_form_of_ABBA(fm_A_explicit, fm_B_explicit, fm_C, fm_sqrt_A_minus_B, &
                                            fm_inv_sqrt_A_minus_B, unit_nr, mp2_env, 0.0_dp)
         CALL cp_fm_create(fm_eigvec, fm_C%matrix_struct)
         CALL choose_eigv_solver(fm_C, fm_eigvec, ref_ens, diag_info)
         IF (diag_info /= 0) CPABORT("Reference diagonalization of C failed in the BSE Davidson check")
         IF (ref_ens(1) <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "Reference matrix C has a non-positive eigenvalue in the BSE Davidson check")
         END IF
         ref_ens(:) = SQRT(ref_ens(:))

         n_ref = last_required_root(ref_ens, n_want, n_ov)
         ALLOCATE (sum_ref(n_ov, n_ref), diff_ref(n_ov, n_ref))
         CALL cp_fm_create(fm_work, fm_C%matrix_struct)
         CALL parallel_gemm("N", "N", n_ov, n_ref, n_ov, 1.0_dp, fm_sqrt_A_minus_B, fm_eigvec, 0.0_dp, &
                            fm_work)
         CALL cp_fm_get_submatrix(fm_work, sum_ref, 1, 1, n_ov, n_ref)
         CALL parallel_gemm("N", "N", n_ov, n_ref, n_ov, 1.0_dp, fm_inv_sqrt_A_minus_B, fm_eigvec, &
                            0.0_dp, fm_work)
         CALL cp_fm_get_submatrix(fm_work, diff_ref, 1, 1, n_ov, n_ref)
         CALL cp_fm_release(fm_work)
         CALL cp_fm_release(fm_C)
         CALL cp_fm_release(fm_sqrt_A_minus_B)
         CALL cp_fm_release(fm_inv_sqrt_A_minus_B)

         ALLOCATE (X_ref(n_ov, n_ref), Y_ref(n_ov, n_ref), Y_dav(n_ov, n_want))
         DO k = 1, n_ref
            sum_ref(:, k) = sum_ref(:, k)/SQRT(ref_ens(k))
            diff_ref(:, k) = diff_ref(:, k)*SQRT(ref_ens(k))
         END DO
         X_ref(:, :) = 0.5_dp*(sum_ref(:, :) + diff_ref(:, :))
         Y_ref(:, :) = 0.5_dp*(sum_ref(:, :) - diff_ref(:, :))
         DEALLOCATE (sum_ref, diff_ref)
         CALL cp_fm_get_submatrix(fm_Y, Y_dav)
      ELSE
         ! the eigensolver destroys its input
         CALL cp_fm_create(fm_work, fm_A_explicit%matrix_struct)
         CALL cp_fm_to_fm(fm_A_explicit, fm_work)
         CALL cp_fm_create(fm_eigvec, fm_A_explicit%matrix_struct)
         CALL choose_eigv_solver(fm_work, fm_eigvec, ref_ens, diag_info)
         IF (diag_info /= 0) CPABORT("Reference diagonalization of A failed in the BSE Davidson check")
         CALL cp_fm_release(fm_work)

         n_ref = last_required_root(ref_ens, n_want, n_ov)
         ALLOCATE (X_ref(n_ov, n_ref))
         CALL cp_fm_get_submatrix(fm_eigvec, X_ref, 1, 1, n_ov, n_ref)
      END IF
      CALL cp_fm_release(fm_eigvec)

      ALLOCATE (X_dav(n_ov, n_want))
      CALL cp_fm_get_submatrix(fm_X, X_dav)

      ! Y_ref and Y_dav stay unallocated for the TDA and then count as absent (Fortran 2008)
      CALL min_multiplet_overlap(ref_ens, n_want, n_ref, X_ref, X_dav, min_overlap, Y_ref, Y_dav)
      CALL print_refcheck(MAXVAL(ABS(exc_ens(:) - ref_ens(1:n_want))), min_overlap, unit_nr)

      DEALLOCATE (ref_ens, X_ref, X_dav)
      IF (ALLOCATED(Y_ref)) DEALLOCATE (Y_ref, Y_dav)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_refcheck

! **************************************************************************************************
!> \brief Smallest singular value over the multiplets of O_kl = sum_ia (X_ref,ia^k X_ia^l - Y_ref,ia^k Y_ia^l),
!>        k over a multiplet of the reference and l over the Davidson states inside it. Degenerate
!>        eigenvectors are arbitrary within their multiplet, so the vectors are not compared one by one:
!>        every singular value of O is 1 when the Davidson states lie in the span of the multiplet
!> \param ref_ens reference energies in Hartree, ascending
!> \param n_want Davidson states compared
!> \param n_ref reference states up to the end of the multiplet that contains state n_want
!> \param X_ref reference X, N_ov x n_ref
!> \param X_dav Davidson X, N_ov x n_want
!> \param min_overlap the smallest singular value found, 1 if every multiplet is reproduced
!> \param Y_ref reference Y for ABBA, with Y_dav
!> \param Y_dav Davidson Y for ABBA
! **************************************************************************************************
   SUBROUTINE min_multiplet_overlap(ref_ens, n_want, n_ref, X_ref, X_dav, min_overlap, Y_ref, Y_dav)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ref_ens
      INTEGER, INTENT(IN)                                :: n_want, n_ref
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: X_ref, X_dav
      REAL(KIND=dp), INTENT(OUT)                         :: min_overlap
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: Y_ref, Y_dav

      INTEGER                                            :: info, lwork, mult_first, mult_last, &
                                                            mult_last_dav
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_vals, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap
      REAL(KIND=dp), DIMENSION(1, 1)                     :: dummy

      min_overlap = 1.0_dp
      mult_first = 1
      DO WHILE (mult_first <= n_want)
         mult_last = mult_first
         DO WHILE (mult_last < n_ref)
            IF (ref_ens(mult_last + 1) - ref_ens(mult_last) >= deg_thresh) EXIT
            mult_last = mult_last + 1
         END DO
         mult_last_dav = MIN(mult_last, n_want)
         ALLOCATE (overlap(mult_last - mult_first + 1, mult_last_dav - mult_first + 1), sing_vals(mult_last_dav - mult_first + 1))
         overlap(:, :) = MATMUL(TRANSPOSE(X_ref(:, mult_first:mult_last)), X_dav(:, mult_first:mult_last_dav))
         IF (PRESENT(Y_ref)) THEN
            overlap(:, :) = overlap(:, :) - MATMUL(TRANSPOSE(Y_ref(:, mult_first:mult_last)), Y_dav(:, mult_first:mult_last_dav))
         END IF
         ! DGESVD wants lwork >= MAX(3 n + m, 5 n) for the m x n overlap with n <= m; 5 m + 10 covers both
         lwork = 5*(mult_last - mult_first + 1) + 10
         ALLOCATE (work(lwork))
 CALL DGESVD('N', 'N', mult_last - mult_first + 1, mult_last_dav - mult_first + 1, overlap, mult_last - mult_first + 1, sing_vals, &
                     dummy, 1, dummy, 1, work, lwork, info)
         IF (info /= 0) CPABORT("SVD failed in the BSE Davidson check")
         min_overlap = MIN(min_overlap, MINVAL(sing_vals))
         DEALLOCATE (overlap, sing_vals, work)
         mult_first = mult_last + 1
      END DO

   END SUBROUTINE min_multiplet_overlap

! **************************************************************************************************
!> \brief The two lines of a reference check
!> \param dev_E largest energy deviation in Hartree
!> \param min_overlap smallest singular value from min_multiplet_overlap
!> \param unit_nr output unit, positive on the writing rank only
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

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

   USE bse_matvec,                      ONLY: bse_matvec_apply,&
                                              bse_matvec_diagonal,&
                                              bse_matvec_env_type,&
                                              bse_matvec_vector_struct
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_iter_and_cond,&
                                              bse_iter_en_cond,&
                                              bse_iter_or_cond,&
                                              bse_iter_res_cond
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE physcon,                         ONLY: evolt
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_davidson'

   ! energies closer than this belong to one multiplet (Hartree)
   REAL(KIND=dp), PARAMETER, PRIVATE :: deg_thresh = 1.0E-6_dp

   PUBLIC :: bse_davidson_tda, bse_davidson_refcheck_tda

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
                                                            n_act, n_conv, n_kernel, n_ov, n_req, &
                                                            n_restart, n_want, nt
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected
      LOGICAL                                            :: all_conv, has_prev, ok, stalled
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: conv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, diag, res, theta, theta_prev
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C, C_prev, C_theta, G, G_new, Q
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_AZ, fm_work, fm_Z
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      para_env => mv_env%para_env
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

      ALLOCATE (G(m_max, m_max), C(m_max, m_max), C_prev(m_max, n_act), C_theta(m_max, n_act))
      ALLOCATE (Q(m_max, 2*n_act))
      ALLOCATE (theta(m_max), theta_prev(n_act), res(n_act), dE(n_act), conv(n_act), selected(n_act))
      G(:, :) = 0.0_dp
      C_prev(:, :) = 0.0_dp
      theta_prev(:) = 0.0_dp
      has_prev = .FALSE.
      stalled = .FALSE.
      n_restart = 0

      CALL initial_guess(diag, n_act, m_max, deg_thresh, fm_Z, m)
      CALL bse_matvec_apply(mv_env, fm_Z, 1, m, fm_AZ)
      n_kernel = m
      CALL subspace_gram(fm_Z, 1, m, fm_AZ, 1, m, para_env, G(1:m, 1:m))

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'Block Davidson iterations within the TDA:'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A5,T14,A9,T25,A9,T36,A14,T52,A13,T67,A14)') 'BSE|', &
            'Iter.', 'Subspace', 'Converged', 'Max. residual', 'Max. dE (eV)', 'Lowest E (eV)'
      END IF

      all_conv = .FALSE.
      DO iter = 1, mp2_env%bse%max_iter

         CALL solve_reduced(G, m, para_env, theta, C)

         ! thick restart onto the current and the previous Ritz vectors, no kernel application
         IF (m + block_size > m_max .AND. m > n_act .AND. m < n_ov) THEN
            CALL restart_basis(C, C_prev, m, n_act, has_prev, Q, k)
            CALL rotate_in_place(fm_Z, m, Q, k, fm_work)
            CALL rotate_in_place(fm_AZ, m, Q, k, fm_work)
            CALL rotate_reduced(G, m, Q, k)
            m = k
            n_restart = n_restart + 1
            CALL solve_reduced(G, m, para_env, theta, C)
         END IF

         ! r_ia,k = sum_m (A Z)_ia,m c_mk - θ_k sum_m Z_ia,m c_mk for all tracked roots k
         DO k = 1, n_act
            C_theta(1:m, k) = C(1:m, k)*theta(k)
         END DO
         CALL subspace_rotate(fm_AZ, m, C(1:m, 1:n_act), n_act, fm_work, 1.0_dp, 0.0_dp)
         CALL subspace_rotate(fm_Z, m, C_theta(1:m, 1:n_act), n_act, fm_work, -1.0_dp, 1.0_dp)
         CALL column_norms(fm_work, 1, n_act, para_env, res)

         dE(:) = ABS(theta(1:n_act) - theta_prev(:))
         DO k = 1, n_act
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

         ! a multiplet that is cut by NUM_EXC_EN has to converge as a whole
         n_req = n_want
         DO WHILE (n_req < n_act)
            IF (theta(n_req + 1) - theta(n_req) >= deg_thresh) EXIT
            n_req = n_req + 1
         END DO
         ! a subspace that spans all transitions is exact
         IF (m == n_ov) conv(:) = .TRUE.
         n_conv = COUNT(conv(1:n_want))
         all_conv = ALL(conv(1:n_req))

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,I5,T14,I9,T25,I9,T36,ES14.4,T52,ES13.4,T67,F14.6)') 'BSE|', &
               iter, m, n_conv, MAXVAL(res(1:n_req)), MAXVAL(dE(1:n_req))*evolt, theta(1)*evolt
         END IF
         IF (all_conv) EXIT

         IF (iter == mp2_env%bse%max_iter) EXIT

         ! corrections for the lowest unconverged roots, t_ia,k = r_ia,k/(d_ia - θ_k)
         nt = 0
         DO k = 1, n_req
            IF (conv(k)) CYCLE
            ! a root below EPS_RES only waits for its energy criterion
            IF (mp2_env%bse%convergence_criterion == bse_iter_and_cond .AND. &
                res(k) < mp2_env%bse%eps_res) CYCLE
            IF (nt == MIN(block_size, m_max - m)) EXIT
            nt = nt + 1
            selected(nt) = k
         END DO
         CALL davidson_corrections(fm_work, selected(1:nt), theta, diag, fm_Z, m + 1)
         CALL project_out(fm_Z, m + 1, nt, fm_Z, m, para_env)
         CALL drop_small_columns(fm_Z, m + 1, nt, para_env)
         ! no correction left: the open roots only wait for their energy criterion, allowed once
         IF (nt == 0 .AND. stalled) THEN
            CPABORT("BSE Davidson: correction vectors vanished before convergence")
         END IF
         stalled = nt == 0
         CALL cholqr2(fm_Z, m + 1, nt, para_env, ok)
         IF (.NOT. ok) CPABORT("BSE Davidson: orthonormalisation broke down")

         IF (nt > 0) THEN
            CALL bse_matvec_apply(mv_env, fm_Z, m + 1, nt, fm_AZ)
            n_kernel = n_kernel + nt

            ALLOCATE (G_new(m + nt, nt))
            CALL subspace_gram(fm_Z, 1, m + nt, fm_AZ, m + 1, nt, para_env, G_new)
            G(1:m + nt, m + 1:m + nt) = G_new(:, :)
            G(m + 1:m + nt, 1:m) = TRANSPOSE(G_new(1:m, :))
            DEALLOCATE (G_new)
         END IF

         C_prev(:, :) = 0.0_dp
         C_prev(1:m, 1:n_act) = C(1:m, 1:n_act)
         has_prev = .TRUE.
         theta_prev(:) = theta(1:n_act)
         m = m + nt
      END DO

      IF (.NOT. all_conv) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A12,T30,A11,T52,A8,T67,A14)') 'BSE|', &
               'Excitation n', 'Energy (eV)', 'Residual', 'Converged'
            DO k = 1, n_req
               WRITE (unit_nr, '(T2,A4,T7,I12,T27,F14.6,T46,ES14.4,T67,L14)') 'BSE|', &
                  k, theta(k)*evolt, res(k), conv(k)
            END DO
         END IF
         CPABORT("BSE Davidson: MAX_ITER reached before convergence")
      END IF

      IF (mp2_env%bse%bse_debug_print) THEN
         ALLOCATE (G_new(m, m))
         CALL subspace_gram(fm_Z, 1, m, fm_Z, 1, m, para_env, G_new)
         DO k = 1, m
            G_new(k, k) = G_new(k, k) - 1.0_dp
         END DO
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A10,T13,A,T63,ES18.6)') 'BSE|DEBUG|', &
               'Max deviation of the basis from orthonormality', MAXVAL(ABS(G_new))
         END IF
         DEALLOCATE (G_new)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Davidson iterations', iter
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Kernel applications', n_kernel
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Thick restarts', n_restart
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

      ALLOCATE (exc_ens(n_want))
      exc_ens(:) = theta(1:n_want)
      CALL bse_matvec_vector_struct(mv_env, n_want, fm_struct)
      CALL cp_fm_create(fm_X, fm_struct, name="fm_X_davidson", set_zero=.TRUE.)
      CALL cp_fm_struct_release(fm_struct)
      CALL subspace_rotate(fm_Z, m, C(1:m, 1:n_want), n_want, fm_X, 1.0_dp, 0.0_dp)

      CALL cp_fm_release(fm_Z)
      CALL cp_fm_release(fm_AZ)
      CALL cp_fm_release(fm_work)
      DEALLOCATE (G, C, C_prev, C_theta, Q, theta, theta_prev, res, dE, conv, selected, diag)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_tda

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
      IF (mp2_env%bse%block_size < 1) CPABORT("BSE_ITERAT%BLOCK_SIZE must be at least 1")
      IF (mp2_env%bse%max_iter < 1) CPABORT("BSE_ITERAT%MAX_ITER must be at least 1")

      n_want = mp2_env%bse%num_exc_en
      IF (n_want > n_ov) THEN
         CPWARN("BSE_ITERAT%NUM_EXC_EN exceeds the number of transitions and is reduced to it.")
         n_want = n_ov
      END IF
      n_act = MIN(n_want + mp2_env%bse%num_buffer_states, n_ov)
      block_size = mp2_env%bse%block_size
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
! **************************************************************************************************
   SUBROUTINE subspace_ceiling(mp2_env, n_ov, n_act, n_arrays, para_env, unit_nr, m_max)

      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: n_ov, n_act, n_arrays
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, INTENT(OUT)                               :: m_max

      REAL(KIND=dp), PARAMETER                           :: mem_fraction = 0.5_dp

      CHARACTER(LEN=16)                                  :: mem_str
      LOGICAL                                            :: from_memory
      REAL(KIND=dp)                                      :: GB_per_column, mem_avail_GB, n_columns

      ! one column of length N_ov, spread over all ranks
      GB_per_column = 8.0E-9_dp*REAL(n_ov, dp)/REAL(para_env%num_pe, dp)
      from_memory = .FALSE.

      IF (mp2_env%bse%max_subspace_factor == -1) THEN
         from_memory = .TRUE.
         CALL mp_mem_avail_per_rank_GB(para_env, mem_avail_GB)
         IF (mem_avail_GB > 0.0_dp) THEN
            ! work array of 2*n_act columns besides the n_arrays subspace arrays
            n_columns = mem_fraction*mem_avail_GB/GB_per_column - 2.0_dp*REAL(n_act, dp)
            n_columns = MIN(n_columns/REAL(n_arrays, dp), 20.0_dp*REAL(n_act, dp))
            IF (n_columns < 3.0_dp*REAL(n_act, dp) .AND. 3*n_act <= n_ov) THEN
               WRITE (mem_str, '(F16.3)') GB_per_column*REAL((3*n_arrays + 2)*n_act, dp)/mem_fraction
               CALL cp_abort(__LOCATION__, &
                             "BSE Davidson: not enough memory for the smallest subspace of three "// &
                             "vectors per state. Free memory needed per MPI rank (GB): "// &
                             TRIM(ADJUSTL(mem_str)))
            END IF
            m_max = INT(n_columns)
         ELSE
            ! free memory not detectable on this system
            m_max = 3*n_act
         END IF
      ELSE IF (mp2_env%bse%max_subspace_factor >= 3) THEN
         m_max = mp2_env%bse%max_subspace_factor*n_act
      ELSE
         CPABORT("BSE_ITERAT%MAX_SUBSPACE_FACTOR must be at least 3 (or -1)")
      END IF
      m_max = MIN(MAX(m_max, 3*n_act), n_ov)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (from_memory) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
               'Maximum subspace dimension (from free memory)', m_max
            WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', &
               'Fraction of the free memory made available', mem_fraction
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
               'Maximum subspace dimension (MAX_SUBSPACE_FACTOR)', m_max
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', &
            'Memory of the subspace arrays per MPI rank (GB)', &
            GB_per_column*REAL(n_arrays*m_max + MIN(2*n_act, m_max), dp)
      END IF

   END SUBROUTINE subspace_ceiling

! **************************************************************************************************
!> \brief Unit vectors on the lowest entries of the diagonal, never cutting a degenerate group
!> \param diag ...
!> \param n_guess ...
!> \param n_max largest number of guess vectors
!> \param deg_thresh_diag ...
!> \param fm_Z ...
!> \param n_used ...
! **************************************************************************************************
   SUBROUTINE initial_guess(diag, n_guess, n_max, deg_thresh_diag, fm_Z, n_used)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: diag
      INTEGER, INTENT(IN)                                :: n_guess, n_max
      REAL(KIND=dp), INTENT(IN)                          :: deg_thresh_diag
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(OUT)                               :: n_used

      INTEGER                                            :: iloc, k, n_ov, nrow_local
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

      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      DO k = 1, n_used
         DO iloc = 1, nrow_local
            IF (row_indices(iloc) == idx(k)) fm_Z%local_data(iloc, k) = 1.0_dp
         END DO
      END DO

      DEALLOCATE (diag_sorted, idx)

   END SUBROUTINE initial_guess

! **************************************************************************************************
!> \brief Eigenpairs of the symmetrised reduced matrix, solved on one rank and broadcast
!> \param G ...
!> \param m ...
!> \param para_env ...
!> \param theta ...
!> \param C ...
! **************************************************************************************************
   SUBROUTINE solve_reduced(G, m, para_env, theta, C)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: G
      INTEGER, INTENT(IN)                                :: m
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: theta
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: C

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_reduced'

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

   END SUBROUTINE solve_reduced

! **************************************************************************************************
!> \brief Orthonormal coefficient vectors of the restart basis: current Ritz vectors, then the
!>        previous ones (modified Gram-Schmidt in the replicated m-space)
!> \param C current Ritz coefficients
!> \param C_prev Ritz coefficients of the previous iteration, zero-padded to m rows
!> \param m ...
!> \param n_act ...
!> \param has_prev ...
!> \param Q m x n_new restart basis
!> \param n_new ...
! **************************************************************************************************
   SUBROUTINE restart_basis(C, C_prev, m, n_act, has_prev, Q, n_new)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: C, C_prev
      INTEGER, INTENT(IN)                                :: m, n_act
      LOGICAL, INTENT(IN)                                :: has_prev
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: Q
      INTEGER, INTENT(OUT)                               :: n_new

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'restart_basis'
      REAL(KIND=dp), PARAMETER                           :: norm_drop = 1.0E-8_dp

      INTEGER                                            :: handle, ipass, j, k
      REAL(KIND=dp)                                      :: norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: basis

      CALL timeset(routineN, handle)

      ALLOCATE (basis(m, 2*n_act))
      basis(:, 1:n_act) = C(1:m, 1:n_act)
      n_new = n_act
      IF (has_prev) THEN
         DO k = 1, n_act
            basis(:, n_new + 1) = C_prev(1:m, k)
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
      END IF
      Q(:, :) = 0.0_dp
      Q(1:m, 1:n_new) = basis(:, 1:n_new)
      DEALLOCATE (basis)

      CALL timestop(handle)

   END SUBROUTINE restart_basis

! **************************************************************************************************
!> \brief out_ia,k = alpha sum_m V_ia,m C_mk + beta out_ia,k for m <= nv, k <= nc on the local rows ia
!> \param fm_V ...
!> \param nv ...
!> \param C ...
!> \param nc ...
!> \param fm_out ...
!> \param alpha ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE subspace_rotate(fm_V, nv, C, nc, fm_out, alpha, beta)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: nv
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(IN)                                      :: C
      INTEGER, INTENT(IN)                                :: nc
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_out
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'subspace_rotate'

      INTEGER                                            :: handle, nrow_local

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local)
      IF (nrow_local > 0) THEN
         CALL DGEMM('N', 'N', nrow_local, nc, nv, alpha, fm_V%local_data, SIZE(fm_V%local_data, 1), &
                    C, SIZE(C, 1), beta, fm_out%local_data, SIZE(fm_out%local_data, 1))
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
      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local)
      fm_V%local_data(1:nrow_local, 1:nq) = fm_work%local_data(1:nrow_local, 1:nq)
      fm_V%local_data(1:nrow_local, nq + 1:m) = 0.0_dp

   END SUBROUTINE rotate_in_place

! **************************************************************************************************
!> \brief G_kl <- sum_pq Q_pk G_pq Q_ql for p, q <= m and k, l <= nq
!> \param G ...
!> \param m ...
!> \param Q ...
!> \param nq ...
! **************************************************************************************************
   SUBROUTINE rotate_reduced(G, m, Q, nq)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: G
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Q
      INTEGER, INTENT(IN)                                :: nq

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G_rot

      ALLOCATE (G_rot(nq, nq))
      G_rot(:, :) = MATMUL(TRANSPOSE(Q(1:m, 1:nq)), MATMUL(G(1:m, 1:m), Q(1:m, 1:nq)))
      G(:, :) = 0.0_dp
      G(1:nq, 1:nq) = G_rot(:, :)
      DEALLOCATE (G_rot)

   END SUBROUTINE rotate_reduced

! **************************************************************************************************
!> \brief G_kl = sum_ia U_ia,u0+k-1 V_ia,v0+l-1 for k <= nu, l <= nv, local product and sum over the ranks
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
!> \brief Euclidean norms of nt columns of a row-distributed matrix
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
!> \param fm_R residuals, column k belongs to root k
!> \param selected ...
!> \param theta ...
!> \param diag ...
!> \param fm_Z ...
!> \param first_col ...
! **************************************************************************************************
   SUBROUTINE davidson_corrections(fm_R, selected, theta, diag, fm_Z, first_col)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_R
      INTEGER, DIMENSION(:), INTENT(IN)                  :: selected
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: theta, diag
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(IN)                                :: first_col

      CHARACTER(LEN=*), PARAMETER :: routineN = 'davidson_corrections'
      REAL(KIND=dp), PARAMETER :: eref_scale = 0.99_dp, threshold = 16.0_dp*EPSILON(1.0_dp)

      INTEGER                                            :: handle, iloc, it, k, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp)                                      :: denom

! near-degeneracy guard as in tddfpt_compute_residual_vects

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      DO it = 1, SIZE(selected)
         k = selected(it)
         DO iloc = 1, nrow_local
            denom = diag(row_indices(iloc)) - theta(k)
            IF (ABS(denom) < threshold) denom = denom + (1.0_dp - eref_scale)*theta(k)
            fm_Z%local_data(iloc, first_col + it - 1) = fm_R%local_data(iloc, k)/denom
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE davidson_corrections

! **************************************************************************************************
!> \brief T_ia,k <- T_ia,k - sum_m V_ia,m sum_jb V_jb,m T_jb,k for m <= nv and nt columns k, twice
!> \param fm_T ...
!> \param first_col ...
!> \param nt ...
!> \param fm_V orthonormal columns
!> \param nv ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE project_out(fm_T, first_col, nt, fm_V, nv, para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col, nt
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      INTEGER, INTENT(IN)                                :: nv
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'project_out'

      INTEGER                                            :: handle, ipass, nrow_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coeff

      IF (nt == 0 .OR. nv == 0) RETURN
      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local)
      ALLOCATE (coeff(nv, nt))
      DO ipass = 1, 2
         CALL subspace_gram(fm_V, 1, nv, fm_T, first_col, nt, para_env, coeff)
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
!> \brief Orthonormalises nt columns by two Cholesky QR passes on the replicated Gram matrix; a
!>        failed Cholesky factorization is repeated once with a diagonal shift and a third pass
!> \param fm_T ...
!> \param first_col ...
!> \param nt ...
!> \param para_env ...
!> \param ok ...
! **************************************************************************************************
   SUBROUTINE cholqr2(fm_T, first_col, nt, para_env, ok)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      INTEGER, INTENT(IN)                                :: first_col, nt
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(OUT)                               :: ok

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cholqr2'

      INTEGER                                            :: handle, info, ipass, k, n_ov, n_pass, &
                                                            nrow_local
      REAL(KIND=dp)                                      :: shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G, R

      ok = .TRUE.
      IF (nt == 0) RETURN
      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm_T, nrow_local=nrow_local, nrow_global=n_ov)
      ALLOCATE (G(nt, nt), R(nt, nt))

      n_pass = 2
      ipass = 0
      DO WHILE (ipass < n_pass)
         ipass = ipass + 1
         CALL subspace_gram(fm_T, first_col, nt, fm_T, first_col, nt, para_env, G)
         R(:, :) = G(:, :)
         info = 0
         IF (para_env%is_source()) CALL DPOTRF('U', nt, R, nt, info)
         CALL para_env%bcast(info)
         IF (info /= 0) THEN
            IF (n_pass == 3) THEN
               ok = .FALSE.
               EXIT
            END IF
            ! shifted Cholesky QR, Fukaya et al., SIAM J. Sci. Comput. (2020), doi:10.1137/18M1218212
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
               ok = .FALSE.
               EXIT
            END IF
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
!> \brief Debug check of the Davidson result against the full diagonalization of the explicit A:
!>        energies, and per multiplet the smallest singular value of X_ref^T X
!> \param fm_A_explicit destroyed on exit
!> \param exc_ens ...
!> \param fm_X ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE bse_davidson_refcheck_tda(fm_A_explicit, exc_ens, fm_X, unit_nr)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_A_explicit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_ens
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_davidson_refcheck_tda'

      INTEGER                                            :: diag_info, handle, info, k1, k2, k2_dav, &
                                                            lwork, n_ov, n_ref, n_want
      REAL(KIND=dp)                                      :: dev_E, min_overlap
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ref_ens, sing_vals, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap, X_dav, X_ref
      REAL(KIND=dp), DIMENSION(1, 1)                     :: dummy
      TYPE(cp_fm_type)                                   :: fm_eigvec

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_A_explicit, nrow_global=n_ov)
      n_want = SIZE(exc_ens)

      CALL cp_fm_create(fm_eigvec, fm_A_explicit%matrix_struct)
      ALLOCATE (ref_ens(n_ov))
      CALL choose_eigv_solver(fm_A_explicit, fm_eigvec, ref_ens, diag_info)
      IF (diag_info /= 0) CPABORT("Reference diagonalization of A failed in the BSE Davidson check")

      dev_E = MAXVAL(ABS(exc_ens(:) - ref_ens(1:n_want)))*evolt

      ! reference columns up to the end of the multiplet that contains state n_want
      n_ref = n_want
      DO WHILE (n_ref < n_ov)
         IF (ref_ens(n_ref + 1) - ref_ens(n_ref) >= deg_thresh) EXIT
         n_ref = n_ref + 1
      END DO
      ALLOCATE (X_ref(n_ov, n_ref), X_dav(n_ov, n_want))
      CALL cp_fm_get_submatrix(fm_eigvec, X_ref, 1, 1, n_ov, n_ref)
      CALL cp_fm_get_submatrix(fm_X, X_dav)
      CALL cp_fm_release(fm_eigvec)

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
         lwork = 5*(k2 - k1 + 1) + 10
         ALLOCATE (work(lwork))
         CALL DGESVD('N', 'N', k2 - k1 + 1, k2_dav - k1 + 1, overlap, k2 - k1 + 1, sing_vals, &
                     dummy, 1, dummy, 1, work, lwork, info)
         IF (info /= 0) CPABORT("SVD failed in the BSE Davidson check")
         min_overlap = MIN(min_overlap, MINVAL(sing_vals))
         DEALLOCATE (overlap, sing_vals, work)
         k1 = k2 + 1
      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A10,T13,A,T63,ES18.6)') 'BSE|DEBUG|', &
            'Max energy deviation vs full diagonalization (eV)', dev_E
         WRITE (unit_nr, '(T2,A10,T13,A,T63,F18.12)') 'BSE|DEBUG|', &
            'Min subspace overlap vs full diagonalization', min_overlap
      END IF

      DEALLOCATE (ref_ens, X_ref, X_dav)

      CALL timestop(handle)

   END SUBROUTINE bse_davidson_refcheck_tda

END MODULE bse_davidson

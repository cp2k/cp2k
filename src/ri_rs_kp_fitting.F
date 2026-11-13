!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief AO sampling, D_llp(R), and the periodic RI-RS fit Z_lP(q).
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE ri_rs_kp_fitting
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE input_constants,                 ONLY: do_potential_short,&
                                              rskp_backend_gemm
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_transform,                ONLY: cell_to_k_transform
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_print_mem_per_rank,&
                                              mp_request_type,&
                                              mp_waitall
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE ri_rs_kp_integrals,              ONLY: compute_d_lP_R,&
                                              compute_d_lP_cutoffs,&
                                              get_local_l_range,&
                                              write_ri_rs_progress
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_kp_fitting'

   ! Complex numbers in the 64 MiB batches of q points.
   INTEGER(KIND=int_8), PARAMETER :: q_batch_elements = 64_int_8*1024_int_8**2/16_int_8

   PUBLIC :: compute_Z_lP_q, &
             evaluate_mu_rl_R_values

CONTAINS

! **************************************************************************************************
!> \brief Build explicit-cell AO values and periodic RI-RS coefficients Z_lP(q).
!> \param bs_env Band-structure environment.
!> \param grid_points Reference-cell RI-RS points r_l.
!> \param Z_lP_q Gaussian RI-RS coefficients on the shifted SCF mesh.
!> \param Z_lP_zero ...
! **************************************************************************************************
   SUBROUTINE compute_Z_lP_q(bs_env, grid_points, &
                             Z_lP_q, Z_lP_zero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: Z_lP_q(:, :, :), Z_lP_zero(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Z_lP_q'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: charge_lq(:, :)
      INTEGER                                            :: first_l, handle, last_l, n_D_cells, &
                                                            n_d_lP_cells, n_grid
      INTEGER, ALLOCATABLE                               :: active_cells(:), cells_R_D_llp(:, :), &
                                                            cells_R_d_lP(:, :)
      REAL(KIND=dp)                                      :: elapsed, start_time
      REAL(KIND=dp), ALLOCATABLE                         :: D_llp_R(:, :, :), d_lP_R(:, :, :), &
                                                            qpoints_head(:, :), qpoints_scf(:, :)

      CALL timeset(routineN, handle)

      start_time = m_walltime()
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A)') "RI-RS: Compute Z_lP(k); required: D_ll'^R and d_lP^R"

      IF (bs_env%ri_rs%tikhonov < 0.0_dp) THEN
         CPABORT("RI_RS%TIKHONOV must not be negative")
      END IF

      n_grid = bs_env%ri_rs%n_grid_points
      CPASSERT(n_grid > 0)
      CPASSERT(SIZE(grid_points, 2) == n_grid)
      CPASSERT(bs_env%nkp_scf_desymm > 0)

      CALL build_cells_R_for_d_lP_R(bs_env, grid_points, cells_R_d_lP)
      CALL build_cells_R_for_D_llp(bs_env, grid_points, cells_R_D_llp, active_cells)
      n_d_lP_cells = SIZE(cells_R_d_lP, 2)
      n_D_cells = SIZE(cells_R_D_llp, 2)
      CPASSERT(n_d_lP_cells > 0)
      CPASSERT(n_D_cells > 0)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') 'RI-RS d_lP^R cells R:', n_d_lP_cells
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') "RI-RS D_ll'^R cells R:", n_D_cells
         FLUSH (bs_env%unit_nr)
      END IF

      ! Each rank stores a rectangular (l,l') block for every translation A.
      CALL get_local_l_range(n_grid, bs_env%para_env%num_pe, bs_env%para_env%mepos, first_l, last_l)
      ALLOCATE (D_llp_R(last_l - first_l + 1, n_grid, n_D_cells))
      CALL compute_D_llp(bs_env, grid_points, cells_R_D_llp, active_cells, first_l, last_l, D_llp_R)
      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, "After distribution of D_ll'^R:")
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS three-center projection d_lP^R: starting'
         FLUSH (bs_env%unit_nr)
      END IF
      ALLOCATE (qpoints_scf(3, bs_env%nkp_scf_desymm))
      qpoints_scf(:, :) = bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm)
      ALLOCATE (qpoints_head(3, 1), SOURCE=0.0_dp)
      CALL compute_d_lP_R(bs_env, grid_points, cells_R_d_lP, d_lP_R, qpoints_head, charge_lq)
      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, 'After RI-RS three-center projection d_lP^R')

      CALL solve_for_Z_lP_q(bs_env, cells_R_D_llp, cells_R_d_lP, D_llp_R, d_lP_R, first_l, last_l, qpoints_scf, &
                            Z_lP_q, distributed=.TRUE.)
      CALL solve_for_Z_lP_q(bs_env, cells_R_D_llp, cells_R_d_lP, D_llp_R, d_lP_R, first_l, last_l, &
                            qpoints_head, Z_lP_zero, distributed=.FALSE., charge_lq=charge_lq)

      DEALLOCATE (cells_R_d_lP, cells_R_D_llp, D_llp_R, qpoints_scf, d_lP_R)
      DEALLOCATE (qpoints_head)

      elapsed = m_walltime() - start_time
      CALL bs_env%para_env%max(elapsed)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') &
         'Completed computing Z_lP(q), Exec. time:', elapsed, ' s'
      IF (bs_env%unit_nr > 0) FLUSH (bs_env%unit_nr)

      CALL timestop(handle)

   END SUBROUTINE compute_Z_lP_q

! **************************************************************************************************
!> \brief Evaluate ϕ_μ(r_l-R) locally for one cell R and the supplied points r_l.
!> \param bs_env ...
!> \param grid_points ...
!> \param i_cell Index of R in index_to_cell_3c.
!> \param ao_values AO values indexed (l,μ), without MPI communication.
! **************************************************************************************************
   SUBROUTINE evaluate_mu_rl_R_values(bs_env, grid_points, i_cell, ao_values)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: i_cell
      REAL(dp), INTENT(OUT)                              :: ao_values(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_mu_rl_R_values'

      INTEGER                                            :: first_ao, handle, i_atom, i_kind, last_ao
      REAL(dp)                                           :: cutoff_squared, source_position(3)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particle_set(:)

      CALL timeset(routineN, handle)

      cell => bs_env%ri_rs%cell
      particle_set => bs_env%ri_rs%particle_set
      CPASSERT(SIZE(ao_values, 1) == SIZE(grid_points, 2))
      CPASSERT(SIZE(ao_values, 2) == bs_env%n_ao)
      ao_values = 0.0_dp
      DO i_atom = 1, bs_env%n_atom
         i_kind = particle_set(i_atom)%atomic_kind%kind_number
         basis => bs_env%basis_set_AO(i_kind)%gto_basis_set
         first_ao = bs_env%i_ao_start_from_atom(i_atom)
         last_ao = bs_env%i_ao_end_from_atom(i_atom)
         source_position = pbc(particle_set(i_atom)%r, cell) + &
                           MATMUL(bs_env%hmat, REAL(bs_env%index_to_cell_3c(1:3, i_cell), dp))
         cutoff_squared = bs_env%ri_rs%radius_ao_per_atom(i_atom)**2
         IF (bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) THEN
            cutoff_squared = bs_env%ri_rs%cutoff_radius_ri_ao**2
         END IF
         CALL evaluate_ao_basis_on_points(ao_values(:, first_ao:last_ao), &
                                          grid_points, basis, source_position, cell, &
                                          cutoff_squared=cutoff_squared, minimum_image=.FALSE.)
      END DO

      CALL timestop(handle)

   END SUBROUTINE evaluate_mu_rl_R_values

! **************************************************************************************************
!> \brief Solve the periodic RI-RS equations at the supplied momentum-transfer points.
!> At zero Tikhonov parameter the Hermitian D matrix is solved in the
!> Moore-Penrose sense, yielding the minimum-norm unregularized fit.
!> \param bs_env ...
!> \param cells_R_D_llp ...
!> \param cells_R_d_lP ...
!> \param D_llp_R ...
!> \param d_lP_R ...
!> \param first_l ...
!> \param last_l ...
!> \param qpoints ...
!> \param Z_lP_q ...
!> \param distributed ...
!> \param charge_lq ...
! **************************************************************************************************
   SUBROUTINE solve_for_Z_lP_q(bs_env, cells_R_D_llp, cells_R_d_lP, D_llp_R, d_lP_R, first_l, last_l, qpoints, &
                               Z_lP_q, distributed, charge_lq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: cells_R_D_llp, cells_R_d_lP
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: D_llp_R(:, :, :), d_lP_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_l, last_l
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: qpoints
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: Z_lP_q(:, :, :)
      LOGICAL, INTENT(IN)                                :: distributed
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: charge_lq(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_for_Z_lP_q'

      COMPLEX(KIND=dp)                                   :: empty_matrix(0, 0), phase
      COMPLEX(KIND=dp), ALLOCATABLE                      :: D_llp_q(:, :), D_llp_q_batch(:, :, :), &
                                                            D_llp_q_local(:, :), d_lP_q(:, :)
      INTEGER :: backend, first_l_rank, first_q, handle, i_batch, i_rank, info, iq, iq_batch, &
         iq_local, iq_store, l, last_l_rank, last_q, lp, n_grid, n_q, n_q_local, n_RI, n_RI_solve, &
         n_sub, next_progress
      INTEGER(KIND=int_8)                                :: n_done, n_local
      INTEGER, ALLOCATABLE                               :: displs(:), recvcounts(:)
      LOGICAL, ALLOCATABLE                               :: active_grid_point(:)
      REAL(KIND=dp)                                      :: max_scale, metric_head, progress_start, &
                                                            volume
      REAL(KIND=dp), ALLOCATABLE                         :: scale(:)

      CALL timeset(routineN, handle)

      n_grid = bs_env%ri_rs%n_grid_points
      n_RI = bs_env%n_RI
      n_RI_solve = n_RI
      IF (PRESENT(charge_lq)) THEN
         IF (bs_env%ri_metric%potential_type /= do_potential_short) THEN
            CPABORT('Unsupported RI metric')
         END IF
         n_RI_solve = n_RI + 1
      END IF
      n_q = SIZE(qpoints, 2)
      backend = MERGE(bs_env%kpoints_scf_desymm%rskp_transform, rskp_backend_gemm, distributed)
      ! q points are cyclic: rank MOD(q,P) for the SCF mesh, MOD(q-1,P) otherwise.
      iq = MODULO(bs_env%para_env%mepos - MERGE(1, 0, distributed), bs_env%para_env%num_pe) + 1
      n_q_local = n_q
      IF (distributed) n_q_local = (n_q + bs_env%para_env%num_pe - iq)/bs_env%para_env%num_pe
      ALLOCATE (Z_lP_q(n_grid, n_RI_solve, n_q_local))
      ! Only ranks solving a q point need a complete D(q).
      IF (iq <= n_q) THEN
         ALLOCATE (D_llp_q(n_grid, n_grid))
      ELSE
         ALLOCATE (D_llp_q(0, 0))
      END IF
      ! D(q) of up to n_sub q points at a time, bounded to 64 MiB unless one q is larger.
      n_sub = INT(MAX(1_int_8, MIN(INT(bs_env%para_env%num_pe, int_8), q_batch_elements/ &
                                   MAX(1_int_8, SIZE(D_llp_R, 1, KIND=int_8)*n_grid))))
      ALLOCATE (D_llp_q_local(n_grid, last_l - first_l + 1), &
                D_llp_q_batch(last_l - first_l + 1, n_grid, n_sub), &
                recvcounts(bs_env%para_env%num_pe), displs(bs_env%para_env%num_pe), &
                d_lP_q(n_grid, n_RI_solve))
      DO i_rank = 0, bs_env%para_env%num_pe - 1
         CALL get_local_l_range(n_grid, bs_env%para_env%num_pe, i_rank, first_l_rank, last_l_rank)
         recvcounts(i_rank + 1) = (last_l_rank - first_l_rank + 1)*n_grid
         displs(i_rank + 1) = SUM(recvcounts(1:i_rank))
      END DO
      ALLOCATE (SCALE(n_grid), active_grid_point(n_grid))
      Z_lP_q = z_zero

      progress_start = m_walltime()
      next_progress = 1
      n_done = 0_int_8
      n_local = 0_int_8
      IF (bs_env%para_env%mepos == 0) THEN
         DO iq = MERGE(bs_env%para_env%num_pe, 1, distributed), n_q, bs_env%para_env%num_pe
            n_local = n_local + 1_int_8
         END DO
      END IF
      iq_local = 0
      DO i_batch = 0, (n_q - 1)/bs_env%para_env%num_pe
         ! D_ll'(q) = sum_A exp(i q.A) D_ll'^A; gather only on the q rank.
         DO first_q = i_batch*bs_env%para_env%num_pe + 1, &
            MIN(n_q, (i_batch + 1)*bs_env%para_env%num_pe), n_sub
            last_q = MIN(first_q + n_sub - 1, n_q, (i_batch + 1)*bs_env%para_env%num_pe)
            CALL cell_to_k_transform(D_llp_R, cells_R_D_llp, -qpoints(:, first_q:last_q), &
                                     bs_env%kpoints_scf_desymm%nkp_grid, &
                                     D_llp_q_batch(:, :, 1:last_q - first_q + 1), backend, &
                                     profile_fft=.TRUE., profile_label='RI-RS D_R -> D_q')
            DO iq_batch = first_q, last_q
               D_llp_q_local(:, :) = TRANSPOSE(D_llp_q_batch(:, :, iq_batch - first_q + 1))
               i_rank = MODULO(iq_batch - 1, bs_env%para_env%num_pe)
               IF (distributed) i_rank = MODULO(iq_batch, bs_env%para_env%num_pe)
               ! Gather conjugate-transposed row blocks directly into contiguous columns.
               IF (bs_env%para_env%mepos == i_rank) THEN
                  CALL bs_env%para_env%gatherv(D_llp_q_local, D_llp_q, recvcounts, displs, i_rank)
               ELSE
                  CALL bs_env%para_env%gatherv(D_llp_q_local, empty_matrix, recvcounts, displs, i_rank)
               END IF
            END DO
         END DO
         IF (PRESENT(charge_lq)) THEN
            CALL fourier_transform_d_lP_R_to_q(bs_env, cells_R_d_lP, d_lP_R, qpoints, i_batch, &
                                               distributed, d_lP_q(:, 2:), iq)
         ELSE
            CALL fourier_transform_d_lP_R_to_q(bs_env, cells_R_d_lP, d_lP_R, qpoints, i_batch, &
                                               distributed, d_lP_q, iq)
         END IF
         IF (iq > n_q) CYCLE
         iq_local = iq_local + 1
         IF (distributed) THEN
            iq_store = iq_local
         ELSE
            iq_store = iq
         END IF
         IF (PRESENT(charge_lq)) THEN
            ! The Gamma-correction scalar head is separate. Keep the full Gaussian RI.
            volume = bs_env%ri_rs%cell%deth
            IF (ALL(bs_env%periodic == [1, 1, 0])) THEN
               ! The separated 2D charge channel is normalized by layer area,
               ! independent of the vacuum-cell height.
               volume = NORM2([bs_env%hmat(2, 1)*bs_env%hmat(3, 2) - &
                               bs_env%hmat(3, 1)*bs_env%hmat(2, 2), &
                               bs_env%hmat(3, 1)*bs_env%hmat(1, 2) - &
                               bs_env%hmat(1, 1)*bs_env%hmat(3, 2), &
                               bs_env%hmat(1, 1)*bs_env%hmat(2, 2) - &
                               bs_env%hmat(2, 1)*bs_env%hmat(1, 2)])
            END IF
            CPASSERT(volume > 0.0_dp)
            metric_head = pi/bs_env%ri_metric%omega**2
            d_lP_q(:, 1) = metric_head/SQRT(volume)*charge_lq(:, iq)
         END IF

         DO lp = 1, n_grid
            DO l = 1, lp
               phase = 0.5_dp*(D_llp_q(l, lp) + CONJG(D_llp_q(lp, l)))
               D_llp_q(l, lp) = phase
               D_llp_q(lp, l) = CONJG(phase)
            END DO
         END DO

         DO l = 1, n_grid
            SCALE(l) = SQRT(MAX(REAL(D_llp_q(l, l), KIND=dp), TINY(1.0_dp)))
         END DO
         max_scale = MAXVAL(scale)
         active_grid_point(:) = scale > &
                                SQRT(MAX(bs_env%eps_filter, EPSILON(1.0_dp)))*max_scale
         IF (.NOT. ANY(active_grid_point)) THEN
            CPABORT("Periodic RI-RS D(q) has no active grid points")
         END IF

         ! Form D'_ll'(q) = D_ll'(q)/(s_l s_l') + eta delta_ll'.
         DO lp = 1, n_grid
            DO l = 1, n_grid
               IF (active_grid_point(l) .AND. active_grid_point(lp)) THEN
                  D_llp_q(l, lp) = D_llp_q(l, lp)/(SCALE(l)*SCALE(lp))
               ELSE
                  D_llp_q(l, lp) = z_zero
               END IF
            END DO
         END DO
         DO l = 1, n_grid
            D_llp_q(l, l) = D_llp_q(l, l) + bs_env%ri_rs%tikhonov*z_one
            IF (active_grid_point(l)) THEN
               d_lP_q(l, :) = d_lP_q(l, :)/SCALE(l)
            ELSE
               d_lP_q(l, :) = z_zero
               SCALE(l) = 1.0_dp
            END IF
         END DO

         IF (bs_env%ri_rs%tikhonov <= 0.0_dp) THEN
            CALL solve_unregularized_D_llp(D_llp_q, d_lP_q)
         ELSE
            CALL zpotrf('L', n_grid, D_llp_q, n_grid, info)
            IF (info == 0) CALL zpotrs('L', n_grid, n_RI_solve, D_llp_q, n_grid, d_lP_q, n_grid, info)
            IF (info /= 0) CPABORT("Periodic RI-RS regularized solve failed")
         END IF

         DO l = 1, n_grid
            Z_lP_q(l, :, iq_store) = d_lP_q(l, :)/SCALE(l)
         END DO
         IF (bs_env%para_env%mepos == 0 .AND. n_local > 1) THEN
            n_done = n_done + 1_int_8
            CALL write_ri_rs_progress(bs_env%unit_nr, 'Z_lP(q) q-point', &
                                      n_done, n_local, next_progress, progress_start)
         END IF
      END DO
      IF (distributed) THEN
         CPASSERT(iq_local == n_q_local)
      END IF
      IF (.NOT. distributed) CALL bs_env%para_env%sum(Z_lP_q)
      DEALLOCATE (D_llp_q, D_llp_q_batch, D_llp_q_local, recvcounts, displs, d_lP_q, scale, active_grid_point)

      CALL timestop(handle)

   END SUBROUTINE solve_for_Z_lP_q

! **************************************************************************************************
!> \brief Solve a possibly rank-deficient Hermitian D matrix without Tikhonov regularization.
!> Only eigenvalues indistinguishable from zero at matrix-roundoff precision
!> are excluded from the Moore-Penrose inverse.
!> \param D_llp_q ...
!> \param rhs ...
! **************************************************************************************************
   SUBROUTINE solve_unregularized_D_llp(D_llp_q, rhs)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: D_llp_q, rhs

      COMPLEX(KIND=dp), ALLOCATABLE                      :: projected_rhs(:, :), work(:)
      INTEGER                                            :: ieigenvalue, info, ngrid
      REAL(KIND=dp)                                      :: rank_tolerance
      REAL(KIND=dp), ALLOCATABLE                         :: eigenvalues(:), rwork(:)

      ngrid = SIZE(D_llp_q, 1)
      CPASSERT(ngrid > 0)
      CPASSERT(SIZE(D_llp_q, 2) == ngrid)
      CPASSERT(SIZE(rhs, 1) == ngrid)
      ALLOCATE (eigenvalues(ngrid), rwork(MAX(1, 3*ngrid - 2)), &
                work(MAX(1, 2*ngrid - 1)), projected_rhs(ngrid, SIZE(rhs, 2)))
      CALL zheev('V', 'L', ngrid, D_llp_q, ngrid, eigenvalues, work, SIZE(work), rwork, info)
      IF (info /= 0) CPABORT('Unregularized periodic RI-RS eigensolver failed')
      rank_tolerance = EPSILON(1.0_dp)*REAL(ngrid, dp)*MAXVAL(ABS(eigenvalues))
      IF (.NOT. ANY(eigenvalues > rank_tolerance)) THEN
         CPABORT('Unregularized RI-RS D matrix has zero numerical rank')
      END IF

      CALL zgemm('C', 'N', ngrid, SIZE(rhs, 2), ngrid, z_one, D_llp_q, ngrid, &
                 rhs, ngrid, z_zero, projected_rhs, ngrid)
      DO ieigenvalue = 1, ngrid
         IF (eigenvalues(ieigenvalue) > rank_tolerance) THEN
            projected_rhs(ieigenvalue, :) = &
               projected_rhs(ieigenvalue, :)/eigenvalues(ieigenvalue)
         ELSE
            projected_rhs(ieigenvalue, :) = z_zero
         END IF
      END DO
      CALL zgemm('N', 'N', ngrid, SIZE(rhs, 2), ngrid, z_one, D_llp_q, ngrid, &
                 projected_rhs, ngrid, z_zero, rhs, ngrid)

   END SUBROUTINE solve_unregularized_D_llp

! **************************************************************************************************
!> \brief Enumerate translations T that survive the R_cut^(lP) criterion for at least one (l,P).
!> \param bs_env ...
!> \param grid_points ...
!> \param cells_R_d_lP ...
! **************************************************************************************************
   SUBROUTINE build_cells_R_for_d_lP_R(bs_env, grid_points, cells_R_d_lP)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells_R_d_lP(:, :)

      INTEGER                                            :: cell_index(3), i_atom, i_l, ix, iy, iz, &
                                                            n_max(3), n_RI, n_work
      INTEGER, ALLOCATABLE                               :: cells_work(:, :)
      LOGICAL                                            :: keep
      REAL(KIND=dp)                                      :: center(3), cutoff, max_cutoff
      REAL(KIND=dp), ALLOCATABLE                         :: automatic_cutoffs(:)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      particle_set => bs_env%ri_rs%particle_set
      CALL compute_d_lP_cutoffs(bs_env, automatic_cutoffs)
      max_cutoff = MAXVAL(automatic_cutoffs)
      DO ix = 1, 3
         IF (bs_env%periodic(ix) == 1) THEN
            ! If |h*T + delta| <= R, then every fractional component obeys
            ! |T_i| <= R*|(h^-1)_i| + |(h^-1 delta)_i|.  Both the atom and the
            ! grid point are mapped into the reference cell, so the last term
            ! is smaller than one.  This reciprocal-vector bound remains safe
            ! for skewed cells, unlike a bound based on direct-vector lengths.
            n_max(ix) = CEILING(max_cutoff*NORM2(bs_env%ri_rs%cell%h_inv(ix, :))) + 1
         ELSE
            n_max(ix) = 0
         END IF
      END DO

      n_work = PRODUCT(2*n_max + 1)
      ALLOCATE (cells_work(3, n_work))
      n_RI = 0
      DO ix = -n_max(1), n_max(1)
         DO iy = -n_max(2), n_max(2)
            DO iz = -n_max(3), n_max(3)
               cell_index = [ix, iy, iz]
               keep = .FALSE.
               DO i_atom = 1, bs_env%n_atom
                  cutoff = automatic_cutoffs(i_atom)
                  center = pbc(particle_set(i_atom)%r, bs_env%ri_rs%cell) + &
                           MATMUL(bs_env%hmat, REAL(cell_index, dp))
                  DO i_l = 1, bs_env%ri_rs%n_grid_points
                     IF (NORM2(center - grid_points(:, i_l)) <= cutoff) THEN
                        keep = .TRUE.
                        EXIT
                     END IF
                  END DO
                  IF (keep) EXIT
               END DO
               IF (.NOT. keep) CYCLE
               n_RI = n_RI + 1
               cells_work(:, n_RI) = cell_index
            END DO
         END DO
      END DO

      ALLOCATE (cells_R_d_lP(3, n_RI))
      cells_R_d_lP(:, :) = cells_work(:, 1:n_RI)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') &
            'Periodic RI-RS automatic R_cut^(lP):', &
            max_cutoff*angstrom, ' Å'
      END IF
      DEALLOCATE (automatic_cutoffs, cells_work)

   END SUBROUTINE build_cells_R_for_d_lP_R

! **************************************************************************************************
!> \brief Enumerate the translations A required by the AO-product D matrix D_ll'^A.
!>
!> Unlike the right-hand side d_lP^T, D_ll'^A is independent of R_cut^(lP). A translation can
!> contribute only if it is the difference of two image cells with nonzero AO values on the
!> reference-cell grid. Keeping a separate cell set prevents Eq. (10) screening from truncating
!> the least-squares metric and avoids enumerating empty three-center image cells.
!> \param bs_env ...
!> \param grid_points Reference-cell RI-RS points.
!> \param cells_R_D_llp ...
!> \param active_cells Indices of cells with nonzero AO values.
! **************************************************************************************************
   SUBROUTINE build_cells_R_for_D_llp(bs_env, grid_points, cells_R_D_llp, active_cells)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells_R_D_llp(:, :), active_cells(:)

      INTEGER                                            :: cell_index(3), i, i_R, i_R_prime, &
                                                            n_active, n_D_cells
      INTEGER, ALLOCATABLE                               :: cells_work(:, :), nonzero_cells(:)
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE                              :: phi(:, :)

      ALLOCATE (nonzero_cells(SIZE(bs_env%index_to_cell_3c, 2)))
      ALLOCATE (phi(bs_env%ri_rs%n_grid_points, bs_env%n_ao))
      n_active = 0
      DO i_R = 1, SIZE(bs_env%index_to_cell_3c, 2)
         CALL evaluate_mu_rl_R_values(bs_env, grid_points, i_R, phi)
         IF (.NOT. ANY(ABS(phi) > 0.0_dp)) CYCLE
         n_active = n_active + 1
         nonzero_cells(n_active) = i_R
      END DO
      CPASSERT(n_active > 0)
      ALLOCATE (active_cells(n_active))
      active_cells(:) = nonzero_cells(:n_active)

      ALLOCATE (cells_work(3, n_active*n_active))
      n_D_cells = 0
      DO i_R = 1, n_active
         DO i_R_prime = 1, n_active
            cell_index = bs_env%index_to_cell_3c(:, active_cells(i_R)) - &
                         bs_env%index_to_cell_3c(:, active_cells(i_R_prime))
            found = .FALSE.
            DO i = 1, n_D_cells
               IF (ALL(cell_index == cells_work(:, i))) THEN
                  found = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. found) THEN
               n_D_cells = n_D_cells + 1
               cells_work(:, n_D_cells) = cell_index
            END IF
         END DO
      END DO

      ALLOCATE (cells_R_D_llp(3, n_D_cells))
      cells_R_D_llp(:, :) = cells_work(:, 1:n_D_cells)
      DEALLOCATE (cells_work, phi)

   END SUBROUTINE build_cells_R_for_D_llp

! **************************************************************************************************
!> \brief Construct D_ll'^R = [Σ_{S-S'=R,μ} ϕ_μ^S(r_L) ϕ_μ^S'(r_L')]².
!> The progress count enumerates pairs of active cell images (S,S'), not Gaussian pairs.
!> \param bs_env ...
!> \param grid_points Reference-cell RI-RS points.
!> \param cells_R_D_llp ...
!> \param active_cells Indices of cells with nonzero AO values.
!> \param first_l ...
!> \param last_l ...
!> \param D_llp_R ...
! **************************************************************************************************
   SUBROUTINE compute_D_llp(bs_env, grid_points, cells_R_D_llp, active_cells, first_l, last_l, D_llp_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: cells_R_D_llp
      INTEGER, INTENT(IN)                                :: active_cells(:), first_l, last_l
      REAL(KIND=dp), CONTIGUOUS, INTENT(OUT)             :: D_llp_R(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_D_llp'

      INTEGER                                            :: cell_A(3), handle, i_A, i_active, &
                                                            i_active_prime, i_R, i_R_prime, &
                                                            next_progress, nrows
      INTEGER(KIND=int_8)                                :: n_done, n_total
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: progress_start
      REAL(KIND=dp), ALLOCATABLE                         :: phi_full(:, :), phi_local(:, :)

      CALL timeset(routineN, handle)

      progress_start = m_walltime()
      next_progress = 1
      n_done = 0_int_8
      n_total = INT(SIZE(active_cells), int_8)**2
      ! Sum the AO images for each owned (l,l') rectangle before squaring.
      D_llp_R = 0.0_dp
      nrows = MAX(0, last_l - first_l + 1)
      IF (first_l <= last_l) THEN
         ALLOCATE (phi_local(nrows, bs_env%n_ao), phi_full(bs_env%ri_rs%n_grid_points, bs_env%n_ao))
         DO i_active_prime = 1, SIZE(active_cells)
            i_R_prime = active_cells(i_active_prime)
            CALL evaluate_mu_rl_R_values(bs_env, grid_points, i_R_prime, phi_full)
            DO i_active = 1, SIZE(active_cells)
               i_R = active_cells(i_active)
               CALL evaluate_mu_rl_R_values(bs_env, grid_points(:, first_l:last_l), i_R, phi_local)
               IF (n_done > 0) CALL write_ri_rs_progress(bs_env%unit_nr, "D_ll'^R Cell-image pairs", &
                                                         n_done, n_total, next_progress, progress_start)
               n_done = n_done + 1
               cell_A = bs_env%index_to_cell_3c(:, i_R) - &
                        bs_env%index_to_cell_3c(:, i_R_prime)
               CALL find_cell(cell_A, cells_R_D_llp, found, i_A)
               IF (.NOT. found) CYCLE
               CALL dgemm('N', 'T', nrows, bs_env%ri_rs%n_grid_points, bs_env%n_ao, &
                          1.0_dp, phi_local, MAX(1, nrows), phi_full, &
                          MAX(1, bs_env%ri_rs%n_grid_points), 1.0_dp, D_llp_R(:, :, i_A), &
                          MAX(1, nrows))
            END DO
         END DO
      END IF
      IF (first_l <= last_l) CALL write_ri_rs_progress(bs_env%unit_nr, "D_ll'^R Cell-image pairs", &
                                                       n_done, n_total, next_progress, progress_start)
      DO i_A = 1, SIZE(D_llp_R, 3)
         D_llp_R(:, :, i_A) = D_llp_R(:, :, i_A)*D_llp_R(:, :, i_A)
      END DO
      IF (ALLOCATED(phi_local)) DEALLOCATE (phi_local, phi_full)

      CALL timestop(handle)

   END SUBROUTINE compute_D_llp

! **************************************************************************************************
!> \brief Find a lattice triplet in an explicit image-cell list.
!> \param cell ...
!> \param index_to_cell ...
!> \param found ...
!> \param index ...
! **************************************************************************************************
   SUBROUTINE find_cell(cell, index_to_cell, found, index)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      LOGICAL, INTENT(OUT)                               :: found
      INTEGER, INTENT(OUT)                               :: index

      INTEGER                                            :: i

      found = .FALSE.
      index = 0
      DO i = 1, SIZE(index_to_cell, 2)
         IF (ALL(cell == index_to_cell(:, i))) THEN
            found = .TRUE.
            index = i
            RETURN
         END IF
      END DO

   END SUBROUTINE find_cell

! **************************************************************************************************
!> \brief Fourier-transform d_lP(R) to d_lP(q) and gather l rows on each q owner.
!> d_lP(q)=Σ_T exp(iqT)d_lP(T); two send buffers overlap transfer with the next q-point transform.
!> Each rank holds all P and T for its l interval; each destination receives one complete q.
!> \param bs_env ...
!> \param cells_R_d_lP ...
!> \param d_lP_R ...
!> \param qpoints ...
!> \param i_batch ...
!> \param distributed ...
!> \param d_lP_q ...
!> \param iq ...
! **************************************************************************************************
   SUBROUTINE fourier_transform_d_lP_R_to_q(bs_env, cells_R_d_lP, d_lP_R, qpoints, i_batch, distributed, &
                                            d_lP_q, iq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: d_lP_R(:, :, :)
      REAL(dp), INTENT(IN)                               :: qpoints(:, :)
      INTEGER, INTENT(IN)                                :: i_batch
      LOGICAL, INTENT(IN)                                :: distributed
      COMPLEX(dp), ASYNCHRONOUS, CONTIGUOUS, INTENT(OUT) :: d_lP_q(:, :)
      INTEGER, INTENT(OUT)                               :: iq

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fourier_transform_d_lP_R_to_q'

      COMPLEX(dp), ALLOCATABLE                           :: result_q(:, :, :)
      COMPLEX(dp), ALLOCATABLE, ASYNCHRONOUS             :: receive_q(:), send_q(:, :, :)
      INTEGER :: backend, first_l, first_q, handle, i_dest, i_l, i_P, i_q, i_rank, i_slot, last_l, &
         last_q, local_first_l, local_last_l, n_active_q, n_l, n_P, n_ranks, n_recv, n_sub, &
         position, q_column, q_offset
      INTEGER, ALLOCATABLE                               :: q_dest(:), q_indices(:)
      LOGICAL                                            :: receiving, sending(2)
      TYPE(mp_request_type)                              :: send_requests(2)
      TYPE(mp_request_type), ALLOCATABLE                 :: recv_requests(:)

      CALL timeset(routineN, handle)

      n_ranks = bs_env%para_env%num_pe
      backend = MERGE(bs_env%kpoints_scf_desymm%rskp_transform, rskp_backend_gemm, distributed)
      i_rank = bs_env%para_env%mepos
      q_offset = 0
      IF (distributed) q_offset = 1
      iq = i_batch*n_ranks + MODULO(i_rank - q_offset, n_ranks) + 1
      CALL get_local_l_range(bs_env%ri_rs%n_grid_points, n_ranks, i_rank, local_first_l, local_last_l)
      CPASSERT(SIZE(d_lP_R, 1) == local_last_l - local_first_l + 1)
      CPASSERT(SIZE(d_lP_R, 2) == SIZE(d_lP_q, 2))
      CPASSERT(SIZE(d_lP_q, KIND=int_8) <= INT(HUGE(0), int_8))
      n_P = bs_env%n_RI
      ALLOCATE (recv_requests(n_ranks), send_q(SIZE(d_lP_R, 1), n_P, 2))
      d_lP_q = z_zero
      n_recv = 0
      receiving = iq <= SIZE(qpoints, 2)
      IF (receiving) THEN
         ! Packed rank blocks: offset=(first_l-1)*n_P, followed by local (l,P) order.
         ! Receive into contiguous storage before assigning the noncontiguous output rows.
         ALLOCATE (receive_q(SIZE(d_lP_q)))
         DO i_dest = 0, n_ranks - 1
            IF (i_dest == i_rank) CYCLE
            CALL get_local_l_range(bs_env%ri_rs%n_grid_points, n_ranks, i_dest, first_l, last_l)
            IF (last_l < first_l .OR. n_P == 0) CYCLE
            n_recv = n_recv + 1
            CALL bs_env%para_env%irecv(receive_q((first_l - 1)*n_P + 1:last_l*n_P), &
                                       i_dest, recv_requests(n_recv), tag=29571)
         END DO
      END IF
      n_active_q = 0
      ALLOCATE (q_indices(n_ranks), q_dest(n_ranks))
      DO i_dest = 0, n_ranks - 1
         i_q = i_batch*n_ranks + MODULO(i_dest - q_offset, n_ranks) + 1
         IF (i_q > SIZE(qpoints, 2)) CYCLE
         n_active_q = n_active_q + 1
         q_indices(n_active_q) = i_q
         q_dest(n_active_q) = i_dest
      END DO
      sending = .FALSE.
      IF (SIZE(d_lP_R, 1) > 0 .AND. n_P > 0) THEN
         ! Sums of up to n_sub q points at a time, bounded to 64 MiB unless one q is larger.
         n_sub = INT(MAX(1_int_8, MIN(INT(n_active_q, int_8), q_batch_elements/ &
                                      (SIZE(d_lP_R, 1, KIND=int_8)*n_P))))
         ALLOCATE (result_q(SIZE(d_lP_R, 1), n_P, n_sub))
         DO first_q = 1, n_active_q, n_sub
            last_q = MIN(first_q + n_sub - 1, n_active_q)
            CALL cell_to_k_transform(d_lP_R, cells_R_d_lP, qpoints(:, q_indices(first_q:last_q)), &
                                     bs_env%kpoints_scf_desymm%nkp_grid, &
                                     result_q(:, :, 1:last_q - first_q + 1), backend, &
                                     profile_fft=.TRUE., profile_label='RI-RS d_lP_R -> d_lP_q')
            DO q_column = first_q, last_q
               i_dest = q_dest(q_column)
               i_slot = MODULO(i_dest, 2) + 1
               IF (sending(i_slot)) CALL mp_waitall(send_requests(i_slot:i_slot))
               sending(i_slot) = .FALSE.
               send_q(:, :, i_slot) = result_q(:, :, q_column - first_q + 1)
               IF (i_dest == i_rank) THEN
                  d_lP_q(local_first_l:local_last_l, :) = send_q(:, :, i_slot)
               ELSE
                  CALL bs_env%para_env%isend(send_q(:, :, i_slot), i_dest, send_requests(i_slot), tag=29571)
                  sending(i_slot) = .TRUE.
               END IF
            END DO
         END DO
         DEALLOCATE (result_q)
      END IF
      DO i_slot = 1, 2
         IF (sending(i_slot)) CALL mp_waitall(send_requests(i_slot:i_slot))
      END DO
      IF (n_recv > 0) CALL mp_waitall(recv_requests(:n_recv))
      IF (receiving) THEN
         DO i_dest = 0, n_ranks - 1
            IF (i_dest == i_rank) CYCLE
            CALL get_local_l_range(bs_env%ri_rs%n_grid_points, n_ranks, i_dest, first_l, last_l)
            n_l = last_l - first_l + 1
            DO i_P = 1, n_P
               position = (first_l - 1)*n_P + (i_P - 1)*n_l
               DO i_l = 1, n_l
                  d_lP_q(first_l + i_l - 1, i_P) = receive_q(position + i_l)
               END DO
            END DO
         END DO
         DEALLOCATE (receive_q)
      END IF
      DEALLOCATE (send_q, recv_requests)

      CALL timestop(handle)

   END SUBROUTINE fourier_transform_d_lP_R_to_q

END MODULE ri_rs_kp_fitting

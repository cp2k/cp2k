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
MODULE ri_rs_kp_compute_Z_lP
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE input_constants,                 ONLY: do_potential_short,&
                                              rskp_backend_gemm
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_lattice_fft,              ONLY: fft_profile_enabled
   USE kpoint_transform,                ONLY: cell_to_k_transform
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_print_mem_per_rank,&
                                              mp_request_type,&
                                              mp_waitall
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE ri_rs_kp_compute_d_lP,           ONLY: compute_d_lP_R
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points,&
                                              get_local_l_range,&
                                              write_ri_rs_progress
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_kp_compute_Z_lP'

   ! Complex numbers in the 64 MiB batches of q points.
   INTEGER(KIND=int_8), PARAMETER :: q_batch_elements = 64_int_8*1024_int_8**2/16_int_8

   TYPE :: ao_cell_type
      INTEGER :: index = 0
      REAL(dp), ALLOCATABLE :: phi_local(:, :)
   END TYPE ao_cell_type

   PUBLIC :: compute_Z_lP_q, &
             evaluate_mu_rl_R

CONTAINS

! **************************************************************************************************
!> \brief Build explicit-cell AO values and periodic RI-RS coefficients Z_lP(q).
!> \param bs_env Band-structure environment.
!> \param grid_points Reference-cell RI-RS points r_l.
!> \param q_slot Compact owner-local index of each retained q, zero for omitted points.
!> \param Z_lP_q Gaussian RI-RS coefficients at the retained q points.
!> \param Z_lP_zero ...
! **************************************************************************************************
   SUBROUTINE compute_Z_lP_q(bs_env, grid_points, q_slot, &
                             Z_lP_q, Z_lP_zero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, INTENT(IN)                                :: q_slot(:)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: Z_lP_q(:, :, :), Z_lP_zero(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Z_lP_q'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: charge_lq(:, :)
      INTEGER                                            :: first_l, handle, last_l, n_grid
      INTEGER, ALLOCATABLE                               :: cells_R_D_llp(:, :), cells_R_d_lP(:, :)
      REAL(KIND=dp)                                      :: elapsed, start_time
      REAL(KIND=dp), ALLOCATABLE                         :: D_llp_R(:, :, :), d_lP_R(:, :, :), &
                                                            qpoints_head(:, :), qpoints_scf(:, :)

      CALL timeset(routineN, handle)

      start_time = m_walltime()
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A)') "RI-RS: Compute Z_lP(k); required: D_ll'^R and d_lP^R"

      n_grid = bs_env%ri_rs%n_grid_points
      CPASSERT(n_grid > 0)
      CPASSERT(SIZE(grid_points, 2) == n_grid)
      CPASSERT(bs_env%nkp_scf_desymm > 0)

      ! Each rank stores a rectangular (l,l') block for every translation A.
      CALL get_local_l_range(n_grid, bs_env%para_env%num_pe, bs_env%para_env%mepos, first_l, last_l)
      CALL build_D_llp_R(bs_env, grid_points, first_l, last_l, cells_R_D_llp, D_llp_R)
      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, "After distribution of D_ll'^R:")
      ALLOCATE (qpoints_scf(3, bs_env%nkp_scf_desymm))
      qpoints_scf(:, :) = bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm)
      ALLOCATE (qpoints_head(3, 1), SOURCE=0.0_dp)
      CALL compute_d_lP_R(bs_env, grid_points, cells_R_d_lP, d_lP_R, qpoints_head, charge_lq)

      CALL solve_for_Z_lP_q(bs_env, cells_R_D_llp, cells_R_d_lP, D_llp_R, d_lP_R, first_l, last_l, qpoints_scf, &
                            q_slot, Z_lP_q, distributed=.TRUE.)
      CALL solve_for_Z_lP_q(bs_env, cells_R_D_llp, cells_R_d_lP, D_llp_R, d_lP_R, first_l, last_l, &
                            qpoints_head, [1], Z_lP_zero, distributed=.FALSE., charge_lq=charge_lq)

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
   SUBROUTINE evaluate_mu_rl_R(bs_env, grid_points, i_cell, ao_values)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: i_cell
      REAL(dp), INTENT(OUT)                              :: ao_values(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'evaluate_mu_rl_R'

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

   END SUBROUTINE evaluate_mu_rl_R

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
!> \param q_slot compact local indices; zero omits a q point
!> \param Z_lP_q ...
!> \param distributed ...
!> \param charge_lq ...
! **************************************************************************************************
   SUBROUTINE solve_for_Z_lP_q(bs_env, cells_R_D_llp, cells_R_d_lP, D_llp_R, d_lP_R, first_l, last_l, qpoints, &
                               q_slot, Z_lP_q, distributed, charge_lq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: cells_R_D_llp, cells_R_d_lP
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: D_llp_R(:, :, :), d_lP_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_l, last_l
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: qpoints
      INTEGER, INTENT(IN)                                :: q_slot(:)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: Z_lP_q(:, :, :)
      LOGICAL, INTENT(IN)                                :: distributed
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: charge_lq(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_for_Z_lP_q'

      COMPLEX(KIND=dp)                                   :: phase
      COMPLEX(KIND=dp), ALLOCATABLE                      :: D_llp_q_batch(:, :, :), d_lP_q(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: D_llp_q(:, :), D_lpl_q_batch(:, :, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: receive_D(:), send_D(:)
      INTEGER :: backend, first_l_rank, first_q, handle, i_batch, i_rank, info, iq, iq_batch, &
         iq_store, l, last_l_rank, last_q, lp, n_active, n_grid, n_q, n_q_local, n_rank, n_RI, &
         n_RI_solve, n_sub, next_progress
      INTEGER(KIND=int_8)                                :: n_done, n_local
      INTEGER, ALLOCATABLE                               :: displs(:), q_batch(:), &
                                                            receive_counts(:), recvcounts(:), &
                                                            send_counts(:), send_offsets(:)
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
      n_rank = bs_env%para_env%num_pe
      CPASSERT(SIZE(q_slot) == n_q)
      backend = MERGE(bs_env%kpoints_scf_desymm%rskp_transform, rskp_backend_gemm, distributed)
      ! q points are cyclic: rank MOD(q,P) for the SCF mesh, MOD(q-1,P) otherwise.
      iq = MODULO(bs_env%para_env%mepos - MERGE(1, 0, distributed), n_rank) + 1
      n_q_local = COUNT(q_slot(iq:n_q:n_rank) > 0)
      ALLOCATE (Z_lP_q(n_grid, n_RI_solve, MERGE(n_q_local, n_q, distributed)))
      ! Only ranks solving a q point need a complete D(q).
      IF (n_q_local > 0) THEN
         ALLOCATE (D_llp_q(n_grid, n_grid))
      ELSE
         ALLOCATE (D_llp_q(0, 0))
      END IF
      ! D(q) of up to n_sub q points at a time, bounded to 64 MiB unless one q is larger. n_sub is
      ! the same on all ranks (largest row block), since each batch is one collective all-to-all.
      n_sub = INT(MAX(1_int_8, MIN(INT(bs_env%para_env%num_pe, int_8), q_batch_elements/ &
                                   MAX(1_int_8, INT((n_grid + bs_env%para_env%num_pe - 1)/ &
                                                    bs_env%para_env%num_pe, int_8)*n_grid))))
      ALLOCATE (D_lpl_q_batch(n_grid, last_l - first_l + 1, n_sub), &
                D_llp_q_batch(last_l - first_l + 1, n_grid, n_sub), &
                recvcounts(bs_env%para_env%num_pe), displs(bs_env%para_env%num_pe), &
                receive_counts(bs_env%para_env%num_pe), send_counts(bs_env%para_env%num_pe), &
                send_offsets(bs_env%para_env%num_pe), d_lP_q(n_grid, n_RI_solve), q_batch(n_rank))
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
      n_local = INT(n_q_local, int_8)
      DO i_batch = 0, (n_q - 1)/n_rank
         n_active = 0
         DO iq_batch = i_batch*n_rank + 1, MIN(n_q, (i_batch + 1)*n_rank)
            IF (q_slot(iq_batch) == 0) CYCLE
            n_active = n_active + 1
            q_batch(n_active) = iq_batch
         END DO
         IF (n_active == 0) CYCLE
         ! D_ll'(q) = sum_A exp(i q.A) D_ll'^A; gather only retained q points.
         DO first_q = 1, n_active, n_sub
            last_q = MIN(first_q + n_sub - 1, n_active)
            CALL cell_to_k_transform(D_llp_R, cells_R_D_llp, -qpoints(:, q_batch(first_q:last_q)), &
                                     bs_env%kpoints_scf_desymm%nkp_grid, &
                                     D_llp_q_batch(:, :, 1:last_q - first_q + 1), backend, &
                                     profile_fft=fft_profile_enabled, profile_label='RI-RS D_R -> D_q')
            ! One all-to-all sends the transposed row blocks of every q to its owner, where they
            ! are contiguous columns of D(q); all owners of the batch receive at the same time.
            send_counts(:) = 0
            send_offsets(:) = 0
            receive_counts(:) = 0
            DO iq_batch = first_q, last_q
               D_lpl_q_batch(:, :, iq_batch - first_q + 1) = TRANSPOSE(D_llp_q_batch(:, :, iq_batch - first_q + 1))
               i_rank = MODULO(q_batch(iq_batch) - 1 + MERGE(1, 0, distributed), n_rank)
               send_counts(i_rank + 1) = SIZE(D_lpl_q_batch(:, :, 1))
               send_offsets(i_rank + 1) = (iq_batch - first_q)*SIZE(D_lpl_q_batch(:, :, 1))
               IF (bs_env%para_env%mepos == i_rank) receive_counts(:) = recvcounts
            END DO
            send_D(1:SIZE(D_lpl_q_batch)) => D_lpl_q_batch
            receive_D(1:SIZE(D_llp_q)) => D_llp_q
            CALL bs_env%para_env%alltoall(send_D, send_counts, send_offsets, &
                                          receive_D, receive_counts, displs)
         END DO
         IF (PRESENT(charge_lq)) THEN
            CALL fourier_transform_d_lP_R_to_q(bs_env, cells_R_d_lP, d_lP_R, qpoints, i_batch, &
                                               distributed, q_slot, d_lP_q(:, 2:), iq)
         ELSE
            CALL fourier_transform_d_lP_R_to_q(bs_env, cells_R_d_lP, d_lP_R, qpoints, i_batch, &
                                               distributed, q_slot, d_lP_q, iq)
         END IF
         IF (iq > n_q) CYCLE
         iq_store = q_slot(iq)
         IF (iq_store == 0) CYCLE
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
      IF (.NOT. distributed) CALL bs_env%para_env%sum(Z_lP_q)
      DEALLOCATE (D_llp_q, D_llp_q_batch, D_lpl_q_batch, recvcounts, displs, receive_counts, send_counts, &
                  send_offsets, d_lP_q, q_batch, scale, active_grid_point)

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
!> \brief Build the AO-product D matrix D_ll'^A on its active translations A.
!>
!> Unlike the right-hand side d_lP^T, D_ll'^A is independent of R_cut^(lP). A translation can
!> contribute only if it is the difference of two image cells with nonzero AO values on the
!> reference-cell grid. Keeping a separate cell set prevents Eq. (10) screening from truncating
!> the least-squares metric and avoids enumerating empty three-center image cells.
!> \param bs_env ...
!> \param grid_points Reference-cell RI-RS points.
!> \param first_l First grid point owned by this rank.
!> \param last_l Last grid point owned by this rank.
!> \param cells_R_D_llp ...
!> \param D_llp_R Local grid rows of D for each translation.
! **************************************************************************************************
   SUBROUTINE build_D_llp_R(bs_env, grid_points, first_l, last_l, cells_R_D_llp, D_llp_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: first_l, last_l
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells_R_D_llp(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: D_llp_R(:, :, :)

      INTEGER                                            :: cell_index(3), extent(3), i_R, &
                                                            i_R_prime, n_active, n_D_cells
      INTEGER, ALLOCATABLE                               :: cell_to_index(:, :, :), cells_work(:, :)
      REAL(dp), ALLOCATABLE                              :: phi(:, :)
      TYPE(ao_cell_type), ALLOCATABLE                    :: active_cells(:)

      ALLOCATE (active_cells(SIZE(bs_env%index_to_cell_3c, 2)))
      ALLOCATE (phi(bs_env%ri_rs%n_grid_points, bs_env%n_ao))
      n_active = 0
      ! Keep the owned rows from the full-grid activity test for the contractions below.
      DO i_R = 1, SIZE(bs_env%index_to_cell_3c, 2)
         CALL evaluate_mu_rl_R(bs_env, grid_points, i_R, phi)
         IF (.NOT. ANY(ABS(phi) > 0.0_dp)) CYCLE
         n_active = n_active + 1
         active_cells(n_active)%index = i_R
         ALLOCATE (active_cells(n_active)%phi_local, SOURCE=phi(first_l:last_l, :))
      END DO
      DEALLOCATE (phi)
      CPASSERT(n_active > 0)

      extent = MAXVAL(bs_env%index_to_cell_3c(:, active_cells(:n_active)%index), DIM=2) - &
               MINVAL(bs_env%index_to_cell_3c(:, active_cells(:n_active)%index), DIM=2)
      ALLOCATE (cell_to_index(-extent(1):extent(1), -extent(2):extent(2), -extent(3):extent(3)), SOURCE=0)
      ALLOCATE (cells_work(3, MIN(SIZE(cell_to_index, KIND=int_8), INT(n_active, int_8)**2)))
      n_D_cells = 0
      DO i_R = 1, n_active
         DO i_R_prime = 1, n_active
            cell_index = bs_env%index_to_cell_3c(:, active_cells(i_R)%index) - &
                         bs_env%index_to_cell_3c(:, active_cells(i_R_prime)%index)
            IF (cell_to_index(cell_index(1), cell_index(2), cell_index(3)) /= 0) CYCLE
            n_D_cells = n_D_cells + 1
            cell_to_index(cell_index(1), cell_index(2), cell_index(3)) = n_D_cells
            cells_work(:, n_D_cells) = cell_index
         END DO
      END DO

      ALLOCATE (cells_R_D_llp(3, n_D_cells))
      cells_R_D_llp(:, :) = cells_work(:, 1:n_D_cells)
      DEALLOCATE (cells_work)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') "RI-RS D_ll'^R cells R:", n_D_cells
         FLUSH (bs_env%unit_nr)
      END IF
      ALLOCATE (D_llp_R(last_l - first_l + 1, SIZE(grid_points, 2), n_D_cells))
      CALL compute_D_llp(bs_env, grid_points, cell_to_index, active_cells(:n_active), D_llp_R)

   END SUBROUTINE build_D_llp_R

! **************************************************************************************************
!> \brief Construct D_ll'^R = [Σ_{S-S'=R,μ} ϕ_μ^S(r_L) ϕ_μ^S'(r_L')]².
!> The progress count enumerates pairs of active cell images (S,S'), not Gaussian pairs.
!> \param bs_env ...
!> \param grid_points Reference-cell RI-RS points.
!> \param cell_to_index Index of each active translation in D_llp_R.
!> \param active_cells Active image cells and their local AO rows.
!> \param D_llp_R ...
! **************************************************************************************************
   SUBROUTINE compute_D_llp(bs_env, grid_points, cell_to_index, active_cells, D_llp_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, ALLOCATABLE, INTENT(IN)                   :: cell_to_index(:, :, :)
      TYPE(ao_cell_type), INTENT(IN)                     :: active_cells(:)
      REAL(KIND=dp), CONTIGUOUS, INTENT(OUT)             :: D_llp_R(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_D_llp'

      INTEGER                                            :: cell_A(3), handle, i_A, i_active, &
                                                            i_active_prime, i_R, i_R_prime, &
                                                            next_progress, nrows
      INTEGER(KIND=int_8)                                :: n_done, n_total
      REAL(KIND=dp)                                      :: progress_start
      REAL(KIND=dp), ALLOCATABLE                         :: phi_full(:, :)

      CALL timeset(routineN, handle)

      progress_start = m_walltime()
      next_progress = 1
      n_done = 0_int_8
      n_total = INT(SIZE(active_cells), int_8)**2
      ! Sum the AO images for each owned (l,l') rectangle before squaring.
      D_llp_R = 0.0_dp
      nrows = SIZE(D_llp_R, 1)
      IF (nrows > 0) THEN
         ALLOCATE (phi_full(bs_env%ri_rs%n_grid_points, bs_env%n_ao))
         DO i_active_prime = 1, SIZE(active_cells)
            i_R_prime = active_cells(i_active_prime)%index
            IF (nrows == SIZE(grid_points, 2)) THEN
               phi_full(:, :) = active_cells(i_active_prime)%phi_local
            ELSE
               CALL evaluate_mu_rl_R(bs_env, grid_points, i_R_prime, phi_full)
            END IF
            DO i_active = 1, SIZE(active_cells)
               i_R = active_cells(i_active)%index
               IF (n_done > 0) CALL write_ri_rs_progress(bs_env%unit_nr, "D_ll'^R Cell-image pairs", &
                                                         n_done, n_total, next_progress, progress_start)
               n_done = n_done + 1
               cell_A = bs_env%index_to_cell_3c(:, i_R) - &
                        bs_env%index_to_cell_3c(:, i_R_prime)
               i_A = cell_to_index(cell_A(1), cell_A(2), cell_A(3))
               CALL dgemm('N', 'T', nrows, bs_env%ri_rs%n_grid_points, bs_env%n_ao, &
                          1.0_dp, active_cells(i_active)%phi_local, MAX(1, nrows), phi_full, &
                          MAX(1, bs_env%ri_rs%n_grid_points), 1.0_dp, D_llp_R(:, :, i_A), &
                          MAX(1, nrows))
            END DO
         END DO
         CALL write_ri_rs_progress(bs_env%unit_nr, "D_ll'^R Cell-image pairs", &
                                   n_done, n_total, next_progress, progress_start)
         DEALLOCATE (phi_full)
      END IF
      DO i_A = 1, SIZE(D_llp_R, 3)
         D_llp_R(:, :, i_A) = D_llp_R(:, :, i_A)*D_llp_R(:, :, i_A)
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_D_llp

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
!> \param q_slot compact local indices; zero omits a q point
!> \param d_lP_q ...
!> \param iq ...
! **************************************************************************************************
   SUBROUTINE fourier_transform_d_lP_R_to_q(bs_env, cells_R_d_lP, d_lP_R, qpoints, i_batch, distributed, &
                                            q_slot, d_lP_q, iq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: d_lP_R(:, :, :)
      REAL(dp), INTENT(IN)                               :: qpoints(:, :)
      INTEGER, INTENT(IN)                                :: i_batch
      LOGICAL, INTENT(IN)                                :: distributed
      INTEGER, INTENT(IN)                                :: q_slot(:)
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
      IF (receiving) receiving = q_slot(iq) > 0
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
         IF (q_slot(i_q) == 0) CYCLE
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
                                     profile_fft=fft_profile_enabled, profile_label='RI-RS d_lP_R -> d_lP_q')
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

END MODULE ri_rs_kp_compute_Z_lP

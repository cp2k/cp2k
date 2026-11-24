!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Analytic d_lP computation for periodic RI-RS with k-points.
!>
!> d_lP^T = Σ_(μR,νS) φ_μR(r_l) φ_νS(r_l) (μR νS|P T) with the erfc RI metric.
!> Translation invariance (μR νS|PT) = (μ0 ν(S-R)|P(T-R)) lets every unique integral block
!> (μ0 νS'|P T') be computed once and contracted for all image cells R and all points r_l.
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE ri_rs_kp_compute_d_lP_analytic
   USE ai_overlap,                      ONLY: overlap_ab
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              real_to_scaled
   USE cryssym,                         ONLY: crys_sym_gen,&
                                              csym_type,&
                                              kpoint_gen,&
                                              release_csym_type
   USE gw_utils_compute_integrals,      ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              rootpi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_mem_avail_per_rank_GB
   USE omp_lib,                         ONLY: omp_destroy_lock,&
                                              omp_get_max_threads,&
                                              omp_get_thread_num,&
                                              omp_init_lock,&
                                              omp_lock_kind,&
                                              omp_set_lock,&
                                              omp_unset_lock
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset,&
                                              nso
   USE orbital_transformation_matrices, ONLY: calculate_rotmat,&
                                              orbrotmat_type,&
                                              release_rotmat
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE ri_rs_kp_compute_d_lP_numeric,   ONLY: compute_d_lP_cutoffs,&
                                              d_lP_integration_accuracy,&
                                              d_lP_quadrature_resolution,&
                                              prepare_ao_images
   USE ri_rs_kp_types,                  ONLY: analytic_pair_type,&
                                              ao_image_type,&
                                              space_group_type
   USE ri_rs_utils,                     ONLY: get_local_l_range,&
                                              write_ri_rs_progress
   USE spglib_f08,                      ONLY: spg_get_major_version
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_kp_compute_d_lP_analytic'

   ! A product φ_μ φ_ν (μ0 νS'|P T') is skipped below this fraction of the d_lP accuracy.
   REAL(dp), PARAMETER :: screening_fraction = 1.0E-3_dp
   ! Shares of the available memory per rank for the accumulator columns (the final group sum
   ! needs a second copy) and for the column slices exchanged within a group per chunk.
   REAL(dp), PARAMETER :: accumulator_share = 0.25_dp, round_share = 0.1_dp
   ! Primitive bounds |x^n exp(-a x^2)| <= C exp(-envelope_fraction a x^2).
   REAL(dp), PARAMETER :: envelope_fraction = 0.75_dp
   ! Tolerance of the space-group search in fractional coordinates.
   REAL(dp), PARAMETER :: eps_symmetry = 1.0E-6_dp
   ! Below this distance*rho, envelope_bound uses the distance -> 0 limit of the erf difference.
   REAL(dp), PARAMETER :: eps_rho_distance = 1.0E-6_dp

   PUBLIC :: analytic_d_lP_is_faster, &
             compute_d_lP_R_analytic

CONTAINS

! **************************************************************************************************
!> \brief Analytic d_lP(T) and the charges Q_l = INT f_l(r)^2 dr, same layout as compute_d_lP_R_numeric.
!> Translation, exchange of μ and ν and the space group of the crystal map every integral block
!> (μ0 νS'|c T') onto a symmetry-unique block. Every
!> symmetry-unique block is assigned to one rank, which computes its integrals once, rotates them
!> into every block of its orbit and contracts each block at once with φ_μ(r_l) φ_ν(r_l);
!> integrals are neither kept nor sent between groups. The accumulator d_lP(P, T, l) is held in
!> full by a group of ranks, the smallest group whose share fits into memory (one rank if the
!> full accumulator fits); within a group of several ranks it is split over the RI columns P and
!> the column slices of a block go to the group rank holding them. Pair data φ_μ(r_l)
!> φ_ν(r_l) of all pairs are on every rank. Finally the groups are summed and the columns
!> gathered onto the l rows.
!> \param bs_env ...
!> \param grid_points ...
!> \param cells_R_d_lP Integer cell vectors T of d_lP^T, indexed (3,T).
!> \param d_lP_R d_lP^T indexed by local point l, RI function P and cell T.
!> \param qpoints Must hold the single point q = 0.
!> \param charge_lq Q_l at q = 0.
! **************************************************************************************************
   SUBROUTINE compute_d_lP_R_analytic(bs_env, grid_points, cells_R_d_lP, d_lP_R, qpoints, charge_lq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: d_lP_R(:, :, :)
      REAL(dp), INTENT(IN)                               :: qpoints(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: charge_lq(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_d_lP_R_analytic'

      INTEGER :: block(9), cell_shift(3), first_l, first_rows, handle, high, i, ichunk, il, irank, &
         it, k, last_l, low, mepos, n_chunks, n_chunks_team, n_cols, n_done, n_keys, n_l, n_owned, &
         n_pairs_kept, n_recv, n_RI, n_send, n_T, n_team, n_teams, next_progress, num_pe, team, &
         team_rank
      INTEGER(int_8)                                     :: n_blocks, n_products, n_unique, n_used, &
                                                            n_values
      INTEGER(int_8), ALLOCATABLE                        :: owned_code(:), pair_key(:)
      INTEGER(KIND=omp_lock_kind)                        :: send_lock
      INTEGER(KIND=omp_lock_kind), ALLOCATABLE           :: locks(:)
      INTEGER, ALLOCATABLE :: ao_first(:), cell_of_T(:, :, :), chunk_last(:), col_first(:), &
         image_of(:, :, :, :), key_list(:, :), order(:), pair_index(:), recv_counts(:), &
         recv_displ(:), recv_head(:, :), ri_first(:), row_of(:, :), send_counts(:), send_displ(:), &
         send_head(:, :)
      REAL(dp)                                           :: accuracy, available_GB, chunk_volume, &
                                                            max_chunk, pair_GB, round_size, &
                                                            start_time, t_loop, threshold
      REAL(dp), ALLOCATABLE :: accumulator(:, :), ao_envelope(:, :), ao_max(:, :), &
         atom_position(:, :), buffer_out(:), buffer_P(:), buffer_ref(:), buffer_sel(:), &
         buffer_unique(:), buffer_w(:), charge(:), charge_thread(:), cutoffs(:), gathered(:), &
         received(:), recv_data(:), ri_bound(:), ri_centers(:, :, :), ri_envelope(:, :), &
         send_data(:), U_atom(:, :, :)
      TYPE(analytic_pair_type), ALLOCATABLE              :: pairs(:)
      TYPE(ao_image_type), ALLOCATABLE                   :: images(:)
      TYPE(gw_3c_ctx_type)                               :: ctx
      TYPE(gw_3c_ws_type), ALLOCATABLE                   :: ws(:)
      TYPE(mp_comm_type)                                 :: cross_comm, team_comm
      TYPE(space_group_type)                             :: group

      CALL timeset(routineN, handle)

      start_time = m_walltime()
      CPASSERT(SIZE(qpoints, 2) == 1)
      CPASSERT(ALL(qpoints == 0.0_dp))
      n_l = bs_env%ri_rs%n_grid_points
      n_RI = bs_env%n_RI
      n_T = SIZE(cells_R_d_lP, 2)
      num_pe = bs_env%para_env%num_pe
      mepos = bs_env%para_env%mepos
      CALL get_local_l_range(n_l, num_pe, mepos, first_l, last_l)
      ALLOCATE (d_lP_R(last_l - first_l + 1, n_RI, n_T), charge_lq(n_l, 1))

      accuracy = d_lP_integration_accuracy(bs_env)
      threshold = screening_fraction*accuracy
      CALL prepare_ao_images(bs_env, grid_points, images)
      CALL compute_d_lP_cutoffs(bs_env, cutoffs)
      CALL index_images_and_cells(bs_env, images, cells_R_d_lP, image_of, row_of, cell_of_T)
      cell_shift = MINVAL(cells_R_d_lP, DIM=2) - 1
      CALL contracted_ri_coefficients(bs_env, U_atom)
      CALL integral_bounds(bs_env, U_atom, ao_max, ri_bound)
      CALL primitive_envelopes(bs_env, U_atom, ao_envelope, ao_first, ri_envelope, ri_first)
      CALL crystal_space_group(bs_env, cells_R_d_lP, group)
      ALLOCATE (atom_position(3, bs_env%n_atom), ri_centers(3, bs_env%n_atom, n_T))
      DO k = 1, bs_env%n_atom
         atom_position(:, k) = pbc(bs_env%ri_rs%particle_set(k)%r, bs_env%ri_rs%cell)
      END DO
      DO it = 1, n_T
         DO k = 1, bs_env%n_atom
            ri_centers(:, k, it) = atom_position(:, k) + MATMUL(bs_env%hmat, REAL(cells_R_d_lP(:, it), dp))
         END DO
      END DO

      ! Smallest group of consecutive ranks whose column share of the accumulator fits into memory;
      ! group rank t holds the RI columns col_first(t) + 1, ..., col_first(t + 1).
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, available_GB)
      DO n_team = 1, num_pe
         IF (MODULO(num_pe, n_team) /= 0) CYCLE
         IF (8.0_dp*REAL(n_l, dp)*REAL(n_T, dp)*REAL((n_RI + n_team - 1)/n_team, dp) <= &
             accumulator_share*available_GB*1.0E9_dp) EXIT
      END DO
      IF (n_team > num_pe) THEN
         n_team = num_pe
         CPWARN("d_lP of the analytic evaluation may exceed the available memory")
      END IF
      n_teams = num_pe/n_team
      team = mepos/n_team
      team_rank = MODULO(mepos, n_team)
      CALL team_comm%from_split(bs_env%para_env, team, key=team_rank)
      CALL cross_comm%from_split(bs_env%para_env, team_rank, key=team)
      ALLOCATE (col_first(0:n_team))
      DO irank = 0, n_team
         col_first(irank) = (n_RI*irank)/n_team
      END DO
      n_cols = col_first(team_rank + 1) - col_first(team_rank)

      ! Pair data of all pairs (μ0 νS') on every rank; this rank keeps its share of the candidate
      ! blocks (c, T') to find the symmetry-unique blocks, and adds its share of the charges.
      CALL list_unique_ao_pairs(bs_env, images, image_of, key_list)
      n_keys = SIZE(key_list, 2)
      ALLOCATE (pairs(n_keys), charge(n_l))
      charge = 0.0_dp
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, grid_points, images, image_of, row_of, cells_R_d_lP, ri_centers, cutoffs, ao_max, &
      !$OMP ri_bound, ao_envelope, ao_first, atom_position, key_list, n_keys, mepos, num_pe, threshold, pairs, &
      !$OMP charge, n_l) &
      !$OMP PRIVATE(k, charge_thread)
      ALLOCATE (charge_thread(n_l))
      charge_thread(:) = 0.0_dp
      !$OMP DO SCHEDULE(DYNAMIC)
      DO k = 1, n_keys
         CALL prepare_ao_pair(bs_env, grid_points, images, image_of, row_of, cells_R_d_lP, ri_centers, &
                              cutoffs, ao_max, ri_bound, ao_envelope, ao_first, atom_position, key_list(:, k), &
                              mepos, num_pe, MODULO(k - 1, num_pe) == mepos, threshold, pairs(k), charge_thread)
      END DO
      !$OMP END DO
      !$OMP CRITICAL(analytic_d_lP_charge_sum)
      charge(:) = charge(:) + charge_thread(:)
      !$OMP END CRITICAL(analytic_d_lP_charge_sum)
      DEALLOCATE (charge_thread)
      !$OMP END PARALLEL
      ! Pairs sorted by the code of (a, b, S'), for the lookup of the pair of a block.
      ALLOCATE (pair_key(n_keys), pair_index(n_keys))
      n_pairs_kept = 0
      pair_GB = 0.0_dp
      DO k = 1, n_keys
         IF (pairs(k)%n_seg == 0) CYCLE
         n_pairs_kept = n_pairs_kept + 1
         pair_key(n_pairs_kept) = block_code(group, [key_list(:, k), 1, 0, 0, 0])
         pair_index(n_pairs_kept) = k
         pair_GB = pair_GB + 4.0_dp*REAL(5*pairs(k)%n_seg, dp)/1.0E9_dp + &
                   8.0_dp*REAL(SIZE(pairs(k)%amax_a) + SIZE(pairs(k)%amax_b) + &
                               SIZE(pairs(k)%pair_envelope) + 7*pairs(k)%n_seg, dp)/1.0E9_dp
      END DO
      ALLOCATE (order(n_pairs_kept))
      CALL sort(pair_key, n_pairs_kept, order)
      pair_index(1:n_pairs_kept) = pair_index(order)
      DEALLOCATE (order)

      CALL owned_unique_blocks(bs_env, group, pairs, owned_code)
      n_owned = SIZE(owned_code)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') 'd_lP analytic space-group operations:', group%n_op
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') "d_lP analytic pairs (mu0 nuS'):", n_pairs_kept
         WRITE (bs_env%unit_nr, '(T2,A,T72,F9.2)') 'd_lP analytic pair data per rank [GB]:', pair_GB
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') 'd_lP analytic ranks per group holding d_lP:', n_team
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') 'd_lP analytic symmetry-unique blocks of rank 0:', n_owned
         FLUSH (bs_env%unit_nr)
      END IF

      ! In a group of several ranks, the owned blocks are processed in chunks whose column slices
      ! for the other group ranks (at most 2 n_op blocks per orbit) fit into round_size.
      round_size = round_share*available_GB*1.0E9_dp/8.0_dp
      ALLOCATE (chunk_last(0:n_owned + 1))
      chunk_last(0) = 0
      n_chunks = 0
      chunk_volume = 0.0_dp
      max_chunk = 0.0_dp
      DO i = 1, n_owned
         block = decoded_block(group, owned_code(i))
         IF (n_team > 1 .AND. chunk_volume > round_size) THEN
            n_chunks = n_chunks + 1
            chunk_last(n_chunks) = i - 1
            chunk_volume = 0.0_dp
         END IF
         chunk_volume = chunk_volume + 2.0_dp*REAL(group%n_op, dp)*REAL(nsgf_atom(bs_env, block(1)), dp)* &
                        REAL(nsgf_atom(bs_env, block(2)), dp)* &
                        REAL(bs_env%i_RI_end_from_atom(block(6)) - bs_env%i_RI_start_from_atom(block(6)) + 1, dp)
         max_chunk = MAX(max_chunk, chunk_volume)
      END DO
      n_chunks = n_chunks + 1
      chunk_last(n_chunks) = n_owned
      n_chunks_team = n_chunks
      IF (n_team > 1) CALL team_comm%max(n_chunks_team)

      ALLOCATE (accumulator(n_cols*n_T, n_l), locks(n_T), send_head(10, 1024), &
                send_data(MERGE(INT(MIN(max_chunk, REAL(HUGE(0), dp))), 1, n_team > 1)))
      accumulator = 0.0_dp
      DO it = 1, n_T
         CALL omp_init_lock(locks(it))
      END DO
      CALL omp_init_lock(send_lock)
      CALL gw_3c_ctx_create(ctx, bs_env, bs_env%ri_metric, basis_j=bs_env%basis_set_AO, &
                            basis_k=bs_env%basis_set_AO, basis_i=bs_env%basis_set_RI)
      n_blocks = 0
      n_unique = 0
      n_used = 0
      n_products = 0
      ALLOCATE (ws(0:omp_get_max_threads() - 1))
      !$OMP PARALLEL DEFAULT(NONE) SHARED(ctx, ws)
      CALL gw_3c_ws_create(ws(omp_get_thread_num()), ctx)
      !$OMP END PARALLEL
      next_progress = 1
      n_done = 0
      t_loop = m_walltime()
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') 'd_lP analytic setup before the first block, time:', &
            t_loop - start_time, ' s'
         FLUSH (bs_env%unit_nr)
      END IF
      DO ichunk = 1, n_chunks_team
         n_send = 0
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
         !$OMP SHARED(bs_env, ctx, ws, group, owned_code, chunk_last, ichunk, n_chunks, pairs, pair_key, &
         !$OMP pair_index, n_pairs_kept, images, row_of, grid_points, cell_of_T, cell_shift, ri_centers, cutoffs, U_atom, &
         !$OMP ri_bound, ri_envelope, ri_first, atom_position, threshold, team_rank, n_team, col_first, locks, &
         !$OMP accumulator, send_lock, n_send, send_head, send_data, n_done, n_owned, next_progress, &
         !$OMP t_loop) &
         !$OMP PRIVATE(i, buffer_unique, buffer_ref, buffer_P, buffer_sel, buffer_w, buffer_out) &
         !$OMP REDUCTION(+:n_blocks, n_unique, n_used, n_products)
         DO i = chunk_last(MIN(ichunk - 1, n_chunks)) + 1, chunk_last(MIN(ichunk, n_chunks))
            CALL contract_orbit_of_unique_block(bs_env, ctx, ws(omp_get_thread_num()), group, owned_code(i), pairs, &
                                                pair_key(1:n_pairs_kept), pair_index(1:n_pairs_kept), images, row_of, &
                                                grid_points, cell_of_T, cell_shift, ri_centers, cutoffs, U_atom, ri_bound, &
                                                ri_envelope, ri_first, atom_position, threshold, team_rank, n_team, col_first, &
                                                locks, accumulator, send_lock, n_send, send_head, send_data, buffer_unique, &
                                                buffer_ref, buffer_P, buffer_sel, buffer_w, buffer_out, n_blocks, n_unique, &
                                                n_used, n_products)
            !$OMP ATOMIC UPDATE
            n_done = n_done + 1
            ! Progress of rank 0, printed by its master thread.
            IF (omp_get_thread_num() == 0 .AND. bs_env%unit_nr > 0) THEN
               CALL write_ri_rs_progress(bs_env%unit_nr, 'd_lP unique blocks done (rank 0)', &
                                         INT(n_done, int_8), INT(n_owned, int_8), next_progress, t_loop)
            END IF
         END DO
         !$OMP END PARALLEL DO

         IF (n_team > 1) THEN
            CALL exchange_slices(team_comm, n_send, send_head, send_data, recv_head, recv_data)
            n_recv = SIZE(recv_head, 2)
            ! Received slices: (pair, atom c, cell T', first and last column, data offset).
            !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
            !$OMP SHARED(bs_env, images, row_of, n_recv, recv_head, recv_data, pairs, grid_points, cell_of_T, &
            !$OMP cell_shift, ri_centers, cutoffs, ri_bound, threshold, col_first, team_rank, &
            !$OMP locks, accumulator) &
            !$OMP PRIVATE(i, buffer_sel, buffer_w, buffer_out) REDUCTION(+:n_used, n_products)
            DO i = 1, n_recv
               CALL contract_block(pairs(recv_head(1, i)), images, row_of, recv_head(2, i), recv_head(3:5, i), &
                                   recv_head(6, i), recv_head(7, i), recv_data(recv_head(8, i) + 1:), &
                                   grid_points, cell_of_T, cell_shift, ri_centers, cutoffs, ri_bound, &
                                   threshold, col_first(team_rank), locks, accumulator, buffer_sel, &
                                   buffer_w, buffer_out, n_used, n_products)
            END DO
            !$OMP END PARALLEL DO
            DEALLOCATE (recv_head, recv_data)
         END IF
      END DO
      CALL write_ri_rs_progress(bs_env%unit_nr, 'd_lP unique blocks done (rank 0)', INT(n_owned, int_8), &
                                INT(n_owned, int_8), next_progress, t_loop)
      !$OMP PARALLEL DEFAULT(NONE) SHARED(ws)
      CALL gw_3c_ws_release(ws(omp_get_thread_num()))
      !$OMP END PARALLEL
      CALL gw_3c_ctx_release(ctx)
      DO it = 1, n_T
         CALL omp_destroy_lock(locks(it))
      END DO
      CALL omp_destroy_lock(send_lock)
      DEALLOCATE (send_head, send_data, pairs)

      ! Sum the groups: member j of cross_comm receives the l rows of the ranks j*n_team, ...,
      ! (j+1)*n_team - 1 in the columns of this group rank.
      ALLOCATE (recv_counts(n_teams))
      DO irank = 0, n_teams - 1
         CALL get_local_l_range(n_l, num_pe, irank*n_team, low, k)
         CALL get_local_l_range(n_l, num_pe, (irank + 1)*n_team - 1, k, high)
         n_values = INT(n_cols*n_T, int_8)*INT(high - low + 1, int_8)
         CPASSERT(n_values <= INT(HUGE(0), int_8))
         recv_counts(irank + 1) = INT(n_values)
      END DO
      ALLOCATE (received(recv_counts(team + 1)))
      CALL cross_comm%sum_scatter(accumulator, received, recv_counts)
      CALL cross_comm%free()
      DEALLOCATE (accumulator)
      ! Within the group, every rank receives the columns of all group ranks for its own l rows.
      CALL get_local_l_range(n_l, num_pe, team*n_team, first_rows, k)
      ALLOCATE (send_counts(n_team), send_displ(n_team), recv_displ(n_team))
      DEALLOCATE (recv_counts)
      ALLOCATE (recv_counts(n_team))
      DO irank = 0, n_team - 1
         CALL get_local_l_range(n_l, num_pe, team*n_team + irank, low, high)
         send_counts(irank + 1) = n_cols*n_T*(high - low + 1)
         send_displ(irank + 1) = n_cols*n_T*(low - first_rows)
         recv_counts(irank + 1) = (col_first(irank + 1) - col_first(irank))*n_T*(last_l - first_l + 1)
      END DO
      recv_displ(1) = 0
      DO irank = 2, n_team
         recv_displ(irank) = recv_displ(irank - 1) + recv_counts(irank - 1)
      END DO
      ALLOCATE (gathered(SUM(recv_counts)))
      CALL team_comm%alltoall(received, send_counts, send_displ, gathered, recv_counts, recv_displ)
      CALL team_comm%free()
      DO irank = 0, n_team - 1
         k = recv_displ(irank + 1)
         DO il = 1, last_l - first_l + 1
            DO it = 1, n_T
               d_lP_R(il, col_first(irank) + 1:col_first(irank + 1), it) = &
                  gathered(k + 1:k + col_first(irank + 1) - col_first(irank))
               k = k + col_first(irank + 1) - col_first(irank)
            END DO
         END DO
      END DO
      CALL bs_env%para_env%sum(charge)
      charge_lq(:, 1) = CMPLX(charge, 0.0_dp, dp)

      CALL bs_env%para_env%sum(n_blocks)
      CALL bs_env%para_env%sum(n_unique)
      CALL bs_env%para_env%sum(n_used)
      CALL bs_env%para_env%sum(n_products)
      start_time = m_walltime() - start_time
      CALL bs_env%para_env%max(start_time)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') 'd_lP analytic symmetry-unique blocks computed:', &
            REAL(n_unique, dp)
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') "d_lP analytic integral blocks (a b S'|c T'):", &
            REAL(n_blocks, dp)
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') 'd_lP analytic column slices used in contraction:', &
            REAL(n_used, dp)
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') 'd_lP analytic contraction flops:', &
            2.0_dp*REAL(n_products, dp)
         WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') &
            'd_lP completed analytic integration, Exec. time:', start_time, ' s'
         FLUSH (bs_env%unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_d_lP_R_analytic

! **************************************************************************************************
!> \brief Pair data of (a 0, b S'): segments, i.e. the image cells R with AO rows r_l shared by the
!> images (a, R) and (b, R+S'), their AO values and screening data, the candidate integral
!> blocks (c, T') of this rank and the contribution to Q_l.
!> \param bs_env ...
!> \param grid_points ...
!> \param images ...
!> \param image_of ...
!> \param row_of ...
!> \param cells_R_d_lP ...
!> \param ri_centers Centers of RI atom c in cell T, ri_centers(:, c, T).
!> \param cutoffs ...
!> \param ao_max ...
!> \param ri_bound ...
!> \param ao_envelope ...
!> \param ao_first ...
!> \param atom_position ...
!> \param key (a, b, S')
!> \param part Only the candidate blocks with MOD(block - 1, n_parts) = part are kept.
!> \param n_parts ...
!> \param with_charge Add the contribution to Q_l.
!> \param threshold ...
!> \param pair ...
!> \param charge ...
! **************************************************************************************************
   SUBROUTINE prepare_ao_pair(bs_env, grid_points, images, image_of, row_of, cells_R_d_lP, ri_centers, &
                              cutoffs, ao_max, ri_bound, ao_envelope, ao_first, atom_position, key, part, &
                              n_parts, with_charge, threshold, pair, charge)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: image_of(:, :, :, :), row_of(:, :), &
                                                            cells_R_d_lP(:, :)
      REAL(dp), INTENT(IN)                               :: ri_centers(:, :, :), cutoffs(:), &
                                                            ao_max(:, :), ri_bound(:), &
                                                            ao_envelope(:, :)
      INTEGER, INTENT(IN)                                :: ao_first(:)
      REAL(dp), INTENT(IN)                               :: atom_position(:, :)
      INTEGER, INTENT(IN)                                :: key(5), part, n_parts
      LOGICAL, INTENT(IN)                                :: with_charge
      REAL(dp), INTENT(IN)                               :: threshold
      TYPE(analytic_pair_type), INTENT(OUT)              :: pair
      REAL(dp), INTENT(INOUT)                            :: charge(:)

      INTEGER :: block_cell(3), c, i, index_3d(3), iseg, it, j, k, lower_T(3), mu, n_a, n_b, &
         n_candidate, n_pair, n_rows, n_seg, nu, shift(3), upper(3), upper_T(3)
      INTEGER, ALLOCATABLE                               :: block_of(:, :, :, :), grown(:, :), &
                                                            kept_cell(:, :), kept_image(:, :), &
                                                            row_l(:), rows_a(:), rows_b(:)
      REAL(dp)                                           :: distance_AB(3)
      REAL(dp), ALLOCATABLE                              :: overlap(:, :), work(:, :)

      pair%atom_a = key(1)
      pair%atom_b = key(2)
      pair%cell = key(3:5)
      n_a = nsgf_atom(bs_env, pair%atom_a)
      n_b = nsgf_atom(bs_env, pair%atom_b)
      pair%n_a = n_a
      pair%n_b = n_b
      pair%fac = 2.0_dp
      IF (pair%atom_a == pair%atom_b .AND. ALL(pair%cell == 0)) pair%fac = 1.0_dp
      upper = [SIZE(image_of, 2), SIZE(image_of, 3), SIZE(image_of, 4)]
      shift = MINVAL(bs_env%index_to_cell_3c, DIM=2) - 1

      ! Segments: image cells R whose images (a, R) and (b, R + S') share grid rows.
      n_seg = COUNT(images(:)%atom == pair%atom_a)
      ALLOCATE (pair%seg_cell(3, n_seg), pair%seg_image(2, n_seg))
      n_seg = 0
      DO i = 1, SIZE(images)
         IF (images(i)%atom /= pair%atom_a) CYCLE
         index_3d = images(i)%cell + pair%cell - shift
         IF (ANY(index_3d < 1) .OR. ANY(index_3d > upper)) CYCLE
         j = image_of(pair%atom_b, index_3d(1), index_3d(2), index_3d(3))
         IF (j == 0) CYCLE
         IF (.NOT. ANY(row_of(:, i) /= 0 .AND. row_of(:, j) /= 0)) CYCLE
         n_seg = n_seg + 1
         pair%seg_image(:, n_seg) = [i, j]
         pair%seg_cell(:, n_seg) = images(i)%cell
      END DO
      pair%n_seg = n_seg
      ALLOCATE (kept_cell(3, n_seg), kept_image(2, n_seg))
      kept_cell(:, :) = pair%seg_cell(:, 1:n_seg)
      kept_image(:, :) = pair%seg_image(:, 1:n_seg)
      CALL MOVE_ALLOC(kept_cell, pair%seg_cell)
      CALL MOVE_ALLOC(kept_image, pair%seg_image)
      IF (n_seg == 0) RETURN

      ! Bounding spheres of the segment rows, the largest AO coefficients per function and
      ! Q_l = Σ fac φ_μ(r_l) S_(μ0,νS') φ_ν(r_l). The rows are not kept.
      ALLOCATE (pair%seg_center(3, n_seg), pair%seg_radius(n_seg), pair%seg_weight(n_seg), &
                pair%seg_amax(n_seg), pair%amax_a(n_a, n_seg), pair%amax_b(n_b, n_seg))
      IF (with_charge) CALL contracted_overlap(bs_env, pair%atom_a, pair%atom_b, pair%cell, overlap)
      DO iseg = 1, n_seg
         CALL shared_grid_rows(images, row_of, pair%seg_image(1, iseg), pair%seg_image(2, iseg), n_rows, &
                               row_l, rows_a, rows_b)
         pair%seg_center(:, iseg) = SUM(grid_points(:, row_l(1:n_rows)), DIM=2)/REAL(n_rows, dp)
         pair%seg_radius(iseg) = 0.0_dp
         DO k = 1, n_rows
            pair%seg_radius(iseg) = MAX(pair%seg_radius(iseg), &
                                        NORM2(grid_points(:, row_l(k)) - pair%seg_center(:, iseg)))
         END DO
         ASSOCIATE (primary_a => images(pair%seg_image(1, iseg))%primary(rows_a(1:n_rows), 1:n_a), &
                    primary_b => images(pair%seg_image(2, iseg))%primary(rows_b(1:n_rows), 1:n_b))
            pair%amax_a(:, iseg) = MAXVAL(ABS(primary_a), DIM=1)
            pair%amax_b(:, iseg) = MAXVAL(ABS(primary_b), DIM=1)
            IF (with_charge) THEN
               ALLOCATE (work(n_rows, n_b))
               CALL dgemm('N', 'N', n_rows, n_b, n_a, 1.0_dp, primary_a, n_rows, overlap, n_a, 0.0_dp, &
                          work, n_rows)
               work(:, :) = work(:, :)*primary_b
               DO k = 1, n_rows
                  charge(row_l(k)) = charge(row_l(k)) + pair%fac*SUM(work(k, :))
               END DO
               DEALLOCATE (work)
            END IF
         END ASSOCIATE
         pair%seg_amax(iseg) = pair%fac*MAXVAL(pair%amax_a(:, iseg))*MAXVAL(pair%amax_b(:, iseg))
         ! |φ_μ φ_ν (μ ν|P)| <= amax_a amax_b min(|μ|_inf, |ν|_inf) ri_bound
         pair%seg_weight(iseg) = 0.0_dp
         DO nu = 1, n_b
            DO mu = 1, n_a
               pair%seg_weight(iseg) = MAX(pair%seg_weight(iseg), pair%amax_a(mu, iseg)*pair%amax_b(nu, iseg)* &
                                           MIN(ao_max(mu, pair%atom_a), ao_max(nu, pair%atom_b)))
            END DO
         END DO
      END DO

      ! Blocks (c, T') with T' = T - R whose RI center T lies within the d_lP cutoff of a segment.
      lower_T = MINVAL(cells_R_d_lP, DIM=2) - MAXVAL(pair%seg_cell, DIM=2)
      upper_T = MAXVAL(cells_R_d_lP, DIM=2) - MINVAL(pair%seg_cell, DIM=2)
      ALLOCATE (block_of(bs_env%n_atom, lower_T(1):upper_T(1), lower_T(2):upper_T(2), lower_T(3):upper_T(3)))
      block_of = 0
      n_candidate = 0
      ALLOCATE (pair%block_list(4, 64))
      DO iseg = 1, n_seg
         DO it = 1, SIZE(cells_R_d_lP, 2)
            block_cell = cells_R_d_lP(:, it) - pair%seg_cell(:, iseg)
            DO c = 1, bs_env%n_atom
               IF (pair%fac*ri_bound(c)*pair%seg_weight(iseg) <= threshold) CYCLE
               IF (block_of(c, block_cell(1), block_cell(2), block_cell(3)) /= 0) CYCLE
               IF (NORM2(ri_centers(:, c, it) - pair%seg_center(:, iseg)) - pair%seg_radius(iseg) > cutoffs(c)) CYCLE
               block_of(c, block_cell(1), block_cell(2), block_cell(3)) = 1
               n_candidate = n_candidate + 1
               IF (MODULO(n_candidate - 1, n_parts) /= part) CYCLE
               pair%n_block = pair%n_block + 1
               IF (pair%n_block > SIZE(pair%block_list, 2)) THEN
                  ALLOCATE (grown(4, 2*SIZE(pair%block_list, 2)))
                  grown(:, 1:pair%n_block - 1) = pair%block_list(:, 1:pair%n_block - 1)
                  CALL MOVE_ALLOC(grown, pair%block_list)
               END IF
               pair%block_list(:, pair%n_block) = [c, block_cell]
            END DO
         END DO
      END DO

      ! Gaussian product envelopes of the primitive pairs: exponent, center, prefactor.
      distance_AB = atom_position(:, pair%atom_b) + MATMUL(bs_env%hmat, REAL(pair%cell, dp)) - &
                    atom_position(:, pair%atom_a)
      ALLOCATE (pair%pair_envelope(5, (ao_first(pair%atom_a + 1) - ao_first(pair%atom_a))* &
                                   (ao_first(pair%atom_b + 1) - ao_first(pair%atom_b))))
      n_pair = 0
      DO i = ao_first(pair%atom_a), ao_first(pair%atom_a + 1) - 1
         DO j = ao_first(pair%atom_b), ao_first(pair%atom_b + 1) - 1
            n_pair = n_pair + 1
            pair%pair_envelope(1, n_pair) = ao_envelope(1, i) + ao_envelope(1, j)
            pair%pair_envelope(2:4, n_pair) = atom_position(:, pair%atom_a) + &
                                              ao_envelope(1, j)*distance_AB/pair%pair_envelope(1, n_pair)
            pair%pair_envelope(5, n_pair) = ao_envelope(2, i)*ao_envelope(2, j)* &
                                            EXP(-ao_envelope(1, i)*ao_envelope(1, j)/ &
                                                pair%pair_envelope(1, n_pair)*SUM(distance_AB**2))
         END DO
      END DO

   END SUBROUTINE prepare_ao_pair

! **************************************************************************************************
!> \brief Integral block (μ0 νS'|P_c T') of a pair in the contracted RI basis, rotated from its
!> symmetry-unique block, split into the column slices of the group ranks: the own slice is
!> contracted, the others are buffered.
!> \param bs_env ...
!> \param group ...
!> \param code Code of the symmetry-unique block.
!> \param op Operation from unique_block for this block.
!> \param unique_integrals Integrals of the symmetry-unique block in the reference RI basis.
!> \param pair ...
!> \param images ...
!> \param row_of ...
!> \param k_pair Index of the pair.
!> \param c RI atom of the block.
!> \param block_cell T'
!> \param amax Largest AO product of the segments of the pair that need the block.
!> \param grid_points ...
!> \param cell_of_T ...
!> \param cell_shift cell_of_T(T - cell_shift) is the index of cell T.
!> \param ri_centers ...
!> \param cutoffs ...
!> \param U_atom ...
!> \param ri_bound ...
!> \param threshold ...
!> \param team_rank ...
!> \param n_team ...
!> \param col_first Team rank t holds the RI columns col_first(t) + 1, ..., col_first(t + 1).
!> \param locks ...
!> \param accumulator ...
!> \param send_lock ...
!> \param n_send ...
!> \param send_head Per slice: team rank, pair, c, T', first and last column, data offset, length.
!> \param send_data ...
!> \param buffer_ref Work buffer of the calling thread.
!> \param buffer_P Work buffer of the calling thread.
!> \param buffer_sel Work buffer of the calling thread.
!> \param buffer_w Work buffer of the calling thread.
!> \param buffer_out Work buffer of the calling thread.
!> \param n_blocks ...
!> \param n_used ...
!> \param n_products ...
! **************************************************************************************************
   SUBROUTINE rotate_and_contract_block(bs_env, group, code, op, unique_integrals, pair, images, row_of, k_pair, c, &
                                        block_cell, &
                                        amax, grid_points, cell_of_T, &
                                        cell_shift, ri_centers, cutoffs, U_atom, ri_bound, threshold, team_rank, &
                                        n_team, col_first, locks, accumulator, &
                                        send_lock, n_send, send_head, send_data, buffer_ref, buffer_P, buffer_sel, &
                                        buffer_w, buffer_out, n_blocks, n_used, n_products)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER(int_8), INTENT(IN)                         :: code
      INTEGER, INTENT(IN)                                :: op
      REAL(dp), INTENT(IN)                               :: unique_integrals(:)
      TYPE(analytic_pair_type), INTENT(IN)               :: pair
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: row_of(:, :), k_pair, c, block_cell(3)
      REAL(dp), INTENT(IN)                               :: amax, grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cell_of_T(:, :, :), cell_shift(3)
      REAL(dp), INTENT(IN)                               :: ri_centers(:, :, :), cutoffs(:), &
                                                            U_atom(:, :, :), ri_bound(:), threshold
      INTEGER, INTENT(IN)                                :: team_rank, n_team, col_first(0:)
      INTEGER(KIND=omp_lock_kind), INTENT(INOUT)         :: locks(:)
      REAL(dp), INTENT(INOUT)                            :: accumulator(:, :)
      INTEGER(KIND=omp_lock_kind), INTENT(INOUT)         :: send_lock
      INTEGER, INTENT(INOUT)                             :: n_send
      INTEGER, ALLOCATABLE, INTENT(INOUT)                :: send_head(:, :)
      REAL(dp), INTENT(INOUT)                            :: send_data(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: buffer_ref(:), buffer_P(:), &
                                                            buffer_sel(:), buffer_w(:), &
                                                            buffer_out(:)
      INTEGER(int_8), INTENT(INOUT)                      :: n_blocks, n_used, n_products

      INTEGER                                            :: first, first_P, irank, last, length, &
                                                            n_ab, n_data, n_P, n_ref
      INTEGER, ALLOCATABLE                               :: grown(:, :)

      n_ab = pair%n_a*pair%n_b
      n_ref = nsgf_ri(bs_env, c)
      first_P = bs_env%i_RI_start_from_atom(c)
      n_P = bs_env%i_RI_end_from_atom(c) - first_P + 1
      CALL ensure_buffer_size(buffer_ref, n_ab*n_ref)
      CALL ensure_buffer_size(buffer_P, n_ab*n_P)
      CALL rotate_block(bs_env, group, code, op, unique_integrals, buffer_ref(1:n_ab*n_ref))
      CALL dgemm('N', 'N', n_ab, n_P, n_ref, 1.0_dp, buffer_ref, n_ab, U_atom(:, :, c), &
                 SIZE(U_atom, 1), 0.0_dp, buffer_P, n_ab)
      n_blocks = n_blocks + 1

      DO irank = 0, n_team - 1
         first = MAX(first_P, col_first(irank) + 1)
         last = MIN(first_P + n_P - 1, col_first(irank + 1))
         IF (first > last) CYCLE
         ASSOCIATE (slice => buffer_P((first - first_P)*n_ab + 1:(last - first_P + 1)*n_ab))
            IF (amax*MAXVAL(ABS(slice)) <= threshold) CYCLE
            IF (irank == team_rank) THEN
               CALL contract_block(pair, images, row_of, c, block_cell, first, last, slice, grid_points, &
                                   cell_of_T, cell_shift, ri_centers, cutoffs, ri_bound, threshold, &
                                   col_first(team_rank), locks, accumulator, buffer_sel, buffer_w, &
                                   buffer_out, n_used, n_products)
            ELSE
               length = SIZE(slice)
               CALL omp_set_lock(send_lock)
               n_data = 0
               IF (n_send > 0) n_data = send_head(9, n_send) + send_head(10, n_send)
               IF (n_send == SIZE(send_head, 2)) THEN
                  ALLOCATE (grown(10, 2*SIZE(send_head, 2)))
                  grown(:, 1:n_send) = send_head(:, 1:n_send)
                  CALL MOVE_ALLOC(grown, send_head)
               END IF
               n_send = n_send + 1
               send_head(:, n_send) = [irank, k_pair, c, block_cell, first, last, n_data, length]
               send_data(n_data + 1:n_data + length) = slice
               CALL omp_unset_lock(send_lock)
            END IF
         END ASSOCIATE
      END DO

   END SUBROUTINE rotate_and_contract_block

! **************************************************************************************************
!> \brief Symmetry-unique blocks owned by this rank: the candidate blocks (c, T') this rank keeps in
!> its pairs are mapped onto their symmetry-unique blocks, whose codes are sent to their owners.
!> \param bs_env ...
!> \param group ...
!> \param pairs Candidate blocks are released on output.
!> \param owned_code Sorted codes of the symmetry-unique blocks of this rank.
! **************************************************************************************************
   SUBROUTINE owned_unique_blocks(bs_env, group, pairs, owned_code)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(space_group_type), INTENT(IN)                 :: group
      TYPE(analytic_pair_type), INTENT(INOUT)            :: pairs(:)
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: owned_code(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'owned_unique_blocks'

      INTEGER                                            :: handle, i, irank, k, n, num_pe, op
      INTEGER(int_8), ALLOCATABLE                        :: code(:), received(:), sent(:)
      INTEGER, ALLOCATABLE                               :: first_code(:), next(:), recv_count(:), &
                                                            recv_displ(:), send_count(:), &
                                                            send_displ(:)

      CALL timeset(routineN, handle)

      num_pe = bs_env%para_env%num_pe
      ALLOCATE (code(MAX(1, SUM(pairs(:)%n_block))), first_code(SIZE(pairs)))
      n = 0
      DO k = 1, SIZE(pairs)
         first_code(k) = n
         n = n + pairs(k)%n_block
      END DO
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) SHARED(group, pairs, code, first_code) PRIVATE(k, i, op)
      DO k = 1, SIZE(pairs)
         DO i = 1, pairs(k)%n_block
            CALL unique_block(group, [pairs(k)%atom_a, pairs(k)%atom_b, pairs(k)%cell, pairs(k)%block_list(:, i)], &
                              code(first_code(k) + i), op)
         END DO
         IF (ALLOCATED(pairs(k)%block_list)) DEALLOCATE (pairs(k)%block_list)
         pairs(k)%n_block = 0
      END DO
      !$OMP END PARALLEL DO
      CALL unique_sorted(code, n)

      ALLOCATE (send_count(num_pe), recv_count(num_pe), send_displ(num_pe), recv_displ(num_pe), next(num_pe))
      send_count = 0
      DO i = 1, n
         irank = owner_of_block(code(i), num_pe) + 1
         send_count(irank) = send_count(irank) + 1
      END DO
      CALL bs_env%para_env%alltoall(send_count, recv_count, 1)
      send_displ(1) = 0
      recv_displ(1) = 0
      DO irank = 2, num_pe
         send_displ(irank) = send_displ(irank - 1) + send_count(irank - 1)
         recv_displ(irank) = recv_displ(irank - 1) + recv_count(irank - 1)
      END DO
      ALLOCATE (sent(MAX(1, n)), received(MAX(1, SUM(recv_count))))
      next(:) = send_displ(:)
      DO i = 1, n
         irank = owner_of_block(code(i), num_pe) + 1
         next(irank) = next(irank) + 1
         sent(next(irank)) = code(i)
      END DO
      CALL bs_env%para_env%alltoall(sent, send_count, send_displ, received, recv_count, recv_displ)
      n = SUM(recv_count)
      CALL unique_sorted(received, n)
      ALLOCATE (owned_code(n))
      owned_code(:) = received(1:n)

      CALL timestop(handle)

   END SUBROUTINE owned_unique_blocks

! **************************************************************************************************
!> \brief Sorts the first n entries of a list and removes repeated entries; n is the new length.
!> \param list ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE unique_sorted(list, n)
      INTEGER(int_8), INTENT(INOUT)                      :: list(:)
      INTEGER, INTENT(INOUT)                             :: n

      INTEGER                                            :: i, m
      INTEGER, ALLOCATABLE                               :: order(:)

      IF (n == 0) RETURN
      ALLOCATE (order(n))
      CALL sort(list, n, order)
      m = 1
      DO i = 2, n
         IF (list(i) == list(m)) CYCLE
         m = m + 1
         list(m) = list(i)
      END DO
      n = m

   END SUBROUTINE unique_sorted

! **************************************************************************************************
!> \brief All blocks (a 0, b S'|c T') of the orbit of one symmetry-unique block: the blocks of pairs
!> that pass the screening are collected, the integrals of the unique block are computed once,
!> rotated into each of these blocks and contracted (rotate_and_contract_block).
!> \param bs_env ...
!> \param ctx ...
!> \param ws Integral workspace of the calling thread.
!> \param group ...
!> \param code Code of the symmetry-unique block.
!> \param pairs ...
!> \param pair_key Sorted codes of the pairs (a, b, S') with segments.
!> \param pair_index Index in pairs of each entry of pair_key.
!> \param images ...
!> \param row_of ...
!> \param grid_points ...
!> \param cell_of_T ...
!> \param cell_shift cell_of_T(T - cell_shift) is the index of cell T.
!> \param ri_centers ...
!> \param cutoffs ...
!> \param U_atom ...
!> \param ri_bound ...
!> \param ri_envelope ...
!> \param ri_first ...
!> \param atom_position ...
!> \param threshold ...
!> \param team_rank ...
!> \param n_team ...
!> \param col_first ...
!> \param locks ...
!> \param accumulator ...
!> \param send_lock ...
!> \param n_send ...
!> \param send_head ...
!> \param send_data ...
!> \param buffer_unique Work buffer of the calling thread.
!> \param buffer_ref Work buffer of the calling thread.
!> \param buffer_P Work buffer of the calling thread.
!> \param buffer_sel Work buffer of the calling thread.
!> \param buffer_w Work buffer of the calling thread.
!> \param buffer_out Work buffer of the calling thread.
!> \param n_blocks ...
!> \param n_unique ...
!> \param n_used ...
!> \param n_products ...
! **************************************************************************************************
   SUBROUTINE contract_orbit_of_unique_block(bs_env, ctx, ws, group, code, pairs, pair_key, pair_index, images, row_of, &
                                             grid_points, cell_of_T, cell_shift, ri_centers, cutoffs, U_atom, ri_bound, &
                                             ri_envelope, ri_first, atom_position, threshold, team_rank, n_team, col_first, &
                                             locks, accumulator, send_lock, n_send, send_head, send_data, buffer_unique, &
                                             buffer_ref, buffer_P, buffer_sel, buffer_w, buffer_out, n_blocks, n_unique, &
                                             n_used, n_products)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(gw_3c_ws_type), INTENT(INOUT)                 :: ws
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER(int_8), INTENT(IN)                         :: code
      TYPE(analytic_pair_type), INTENT(IN)               :: pairs(:)
      INTEGER(int_8), INTENT(IN)                         :: pair_key(:)
      INTEGER, INTENT(IN)                                :: pair_index(:)
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: row_of(:, :)
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cell_of_T(:, :, :), cell_shift(3)
      REAL(dp), INTENT(IN)                               :: ri_centers(:, :, :), cutoffs(:), &
                                                            U_atom(:, :, :), ri_bound(:), &
                                                            ri_envelope(:, :)
      INTEGER, INTENT(IN)                                :: ri_first(:)
      REAL(dp), INTENT(IN)                               :: atom_position(:, :), threshold
      INTEGER, INTENT(IN)                                :: team_rank, n_team, col_first(0:)
      INTEGER(KIND=omp_lock_kind), INTENT(INOUT)         :: locks(:)
      REAL(dp), INTENT(INOUT)                            :: accumulator(:, :)
      INTEGER(KIND=omp_lock_kind), INTENT(INOUT)         :: send_lock
      INTEGER, INTENT(INOUT)                             :: n_send
      INTEGER, ALLOCATABLE, INTENT(INOUT)                :: send_head(:, :)
      REAL(dp), INTENT(INOUT)                            :: send_data(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: buffer_unique(:), buffer_ref(:), &
                                                            buffer_P(:), buffer_sel(:), &
                                                            buffer_w(:), buffer_out(:)
      INTEGER(int_8), INTENT(INOUT)                      :: n_blocks, n_unique, n_used, n_products

      INTEGER                                            :: g, i, k, n_image, n_kept, n_valid, op, &
                                                            position, sgn, unique(9)
      INTEGER(int_8)                                     :: image_unique
      INTEGER, ALLOCATABLE                               :: order(:), valid_cell(:), valid_seg(:)
      LOGICAL                                            :: screened
      REAL(dp)                                           :: bound
      REAL(dp), ALLOCATABLE                              :: work(:, :, :)
      REAL(dp)                                           :: amax(2*group%n_op)
      INTEGER(int_8)                                     :: image_code(2*group%n_op)
      INTEGER                                            :: kept(2*group%n_op), &
                                                            image(9, 2*group%n_op)

      unique = decoded_block(group, code)
      ! The orbit: g(unique block) and its exchange for all operations g.
      n_image = 0
      DO sgn = 1, -1, -2
         DO g = 1, group%n_op
            n_image = n_image + 1
            image(:, n_image) = transformed_block(group, g, unique)
            IF (sgn < 0) image(:, n_image) = [image(2, n_image), image(1, n_image), -image(3:5, n_image), &
                                              image(6, n_image), image(7:9, n_image) - image(3:5, n_image)]
            image_code(n_image) = block_code(group, image(:, n_image))
         END DO
      END DO
      ALLOCATE (order(n_image))
      CALL sort(image_code, n_image, order)
      ! Blocks of pairs with segments within the cutoff whose envelope bound is not negligible; the
      ! bound depends only on the geometry and is the same for all blocks of the orbit.
      n_kept = 0
      bound = -1.0_dp
      DO i = 1, n_image
         IF (i > 1) THEN
            IF (image_code(i) == image_code(i - 1)) CYCLE
         END IF
         ASSOCIATE (b => image(:, order(i)))
            position = FINDLOC(pair_key, block_code(group, [b(1:5), 1, 0, 0, 0]), DIM=1)
            IF (position == 0) CYCLE
            k = pair_index(position)
            CALL segments_needing_block(pairs(k), b(6), b(7:9), cell_of_T, cell_shift, ri_centers, cutoffs, ri_bound, &
                                        threshold, n_valid, valid_seg, valid_cell)
            IF (n_valid == 0) CYCLE
            IF (bound < 0.0_dp) bound = envelope_bound(pairs(k)%pair_envelope, &
                                                       ri_envelope(:, ri_first(b(6)):ri_first(b(6) + 1) - 1), &
                                                       atom_position(:, b(6)) + &
                                                       MATMUL(bs_env%hmat, REAL(b(7:9), dp)), &
                                                       bs_env%ri_metric%omega)
            amax(i) = MAXVAL(pairs(k)%seg_amax(valid_seg(1:n_valid)))
            IF (bound*amax(i) <= threshold) CYCLE
            n_kept = n_kept + 1
            kept(n_kept) = i
         END ASSOCIATE
      END DO
      IF (n_kept == 0) RETURN

      ALLOCATE (work(nsgf_atom(bs_env, unique(1)), nsgf_atom(bs_env, unique(2)), nsgf_ri(bs_env, unique(6))))
      work = 0.0_dp
      CALL build_3c_integral_block_ctx(work, ctx, ws, unique(1), unique(2), unique(6), cell_j=[0, 0, 0], &
                                       cell_k=unique(3:5), cell_i=unique(7:9), screened=screened)
      n_unique = n_unique + 1
      IF (screened) RETURN
      CALL ensure_buffer_size(buffer_unique, SIZE(work))
      buffer_unique(1:SIZE(work)) = RESHAPE(work, [SIZE(work)])

      DO i = 1, n_kept
         ASSOCIATE (b => image(:, order(kept(i))))
            k = pair_index(FINDLOC(pair_key, block_code(group, [b(1:5), 1, 0, 0, 0]), DIM=1))
            CALL unique_block(group, b, image_unique, op)
            CPASSERT(image_unique == code)
            CALL rotate_and_contract_block(bs_env, group, code, op, buffer_unique(1:SIZE(work)), pairs(k), images, row_of, k, &
                                           b(6), b(7:9), amax(kept(i)), grid_points, cell_of_T, cell_shift, ri_centers, &
                                           cutoffs, U_atom, ri_bound, threshold, team_rank, n_team, col_first, locks, &
                                           accumulator, send_lock, n_send, send_head, send_data, buffer_ref, buffer_P, &
                                           buffer_sel, buffer_w, buffer_out, n_blocks, n_used, n_products)
         END ASSOCIATE
      END DO

   END SUBROUTINE contract_orbit_of_unique_block

! **************************************************************************************************
!> \brief Segments of a pair that need the block (c, T'), with the index of their cell T = T' + R.
!> \param pair ...
!> \param c ...
!> \param block_cell T'
!> \param cell_of_T ...
!> \param cell_shift cell_of_T(T - cell_shift) is the index of cell T.
!> \param ri_centers ...
!> \param cutoffs ...
!> \param ri_bound ...
!> \param threshold ...
!> \param n_valid ...
!> \param valid_seg ...
!> \param valid_cell ...
! **************************************************************************************************
   SUBROUTINE segments_needing_block(pair, c, block_cell, cell_of_T, cell_shift, ri_centers, cutoffs, ri_bound, &
                                     threshold, n_valid, valid_seg, valid_cell)
      TYPE(analytic_pair_type), INTENT(IN)               :: pair
      INTEGER, INTENT(IN)                                :: c, block_cell(3), cell_of_T(:, :, :), &
                                                            cell_shift(3)
      REAL(dp), INTENT(IN)                               :: ri_centers(:, :, :), cutoffs(:), &
                                                            ri_bound(:), threshold
      INTEGER, INTENT(OUT)                               :: n_valid
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: valid_seg(:), valid_cell(:)

      INTEGER                                            :: index_3d(3), iseg, it

      ALLOCATE (valid_seg(pair%n_seg), valid_cell(pair%n_seg))
      n_valid = 0
      IF (pair%fac*ri_bound(c)*MAXVAL(pair%seg_weight) <= threshold) RETURN
      DO iseg = 1, pair%n_seg
         index_3d = block_cell + pair%seg_cell(:, iseg) - cell_shift
         IF (ANY(index_3d < 1) .OR. ANY(index_3d > SHAPE(cell_of_T))) CYCLE
         it = cell_of_T(index_3d(1), index_3d(2), index_3d(3))
         IF (it == 0) CYCLE
         IF (NORM2(ri_centers(:, c, it) - pair%seg_center(:, iseg)) - pair%seg_radius(iseg) > cutoffs(c)) CYCLE
         IF (pair%fac*ri_bound(c)*pair%seg_weight(iseg) <= threshold) CYCLE
         n_valid = n_valid + 1
         valid_seg(n_valid) = iseg
         valid_cell(n_valid) = it
      END DO

   END SUBROUTINE segments_needing_block

! **************************************************************************************************
!> \brief Grid rows shared by the AO images i and j: point l and the rows of the AO values in
!> their primary arrays.
!> \param images ...
!> \param row_of Row of point l in the primary array of an image, zero if absent.
!> \param i ...
!> \param j ...
!> \param n_rows ...
!> \param row_l ...
!> \param rows_a Rows in the primary array of image i.
!> \param rows_b Rows in the primary array of image j.
! **************************************************************************************************
   SUBROUTINE shared_grid_rows(images, row_of, i, j, n_rows, row_l, rows_a, rows_b)
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: row_of(:, :), i, j
      INTEGER, INTENT(OUT)                               :: n_rows
      INTEGER, ALLOCATABLE, INTENT(INOUT)                :: row_l(:), rows_a(:), rows_b(:)

      INTEGER                                            :: il, k, n

      ! Rows of image i: packed rows of primary_l, or all rows when the image is not packed.
      IF (ALLOCATED(images(i)%primary_l)) THEN
         n = SIZE(images(i)%primary_l)
      ELSE
         n = SIZE(row_of, 1)
      END IF
      IF (ALLOCATED(row_l)) THEN
         IF (SIZE(row_l) < n) DEALLOCATE (row_l, rows_a, rows_b)
      END IF
      IF (.NOT. ALLOCATED(row_l)) ALLOCATE (row_l(n), rows_a(n), rows_b(n))
      n_rows = 0
      DO k = 1, n
         il = k
         IF (ALLOCATED(images(i)%primary_l)) il = images(i)%primary_l(k)
         IF (row_of(il, i) == 0 .OR. row_of(il, j) == 0) CYCLE
         n_rows = n_rows + 1
         row_l(n_rows) = il
         rows_a(n_rows) = row_of(il, i)
         rows_b(n_rows) = row_of(il, j)
      END DO

   END SUBROUTINE shared_grid_rows

! **************************************************************************************************
!> \brief Adds Σ_(μ ν) fac φ_μ(r_l) φ_ν(r_l) (μ0 νS'|P T') to d_lP(P, T' + R, l) for the
!> RI columns P = first_col, ..., last_col of a block, all segments R and all points r_l inside
!> the d_lP cutoff, with one DGEMM per segment and a per-cell lock.
!> \param pair ...
!> \param images ...
!> \param row_of ...
!> \param c ...
!> \param block_cell ...
!> \param first_col ...
!> \param last_col ...
!> \param integrals (μ0 νS'|P T'), indexed (μ ν, P).
!> \param grid_points ...
!> \param cell_of_T ...
!> \param cell_shift cell_of_T(T - cell_shift) is the index of cell T.
!> \param ri_centers ...
!> \param cutoffs ...
!> \param ri_bound ...
!> \param threshold ...
!> \param col_offset Columns of this rank start after col_offset.
!> \param locks ...
!> \param accumulator ...
!> \param buffer_sel Work buffer of the calling thread.
!> \param buffer_w Work buffer of the calling thread.
!> \param buffer_out Work buffer of the calling thread.
!> \param n_used ...
!> \param n_products ...
! **************************************************************************************************
   SUBROUTINE contract_block(pair, images, row_of, c, block_cell, first_col, last_col, integrals, grid_points, &
                             cell_of_T, cell_shift, ri_centers, cutoffs, ri_bound, threshold, &
                             col_offset, locks, accumulator, buffer_sel, buffer_w, buffer_out, n_used, &
                             n_products)
      TYPE(analytic_pair_type), INTENT(IN)               :: pair
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: row_of(:, :), c, block_cell(3), &
                                                            first_col, last_col
      REAL(dp), INTENT(IN)                               :: integrals(:), grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cell_of_T(:, :, :), cell_shift(3)
      REAL(dp), INTENT(IN)                               :: ri_centers(:, :, :), cutoffs(:), &
                                                            ri_bound(:), threshold
      INTEGER, INTENT(IN)                                :: col_offset
      INTEGER(KIND=omp_lock_kind), INTENT(INOUT)         :: locks(:)
      REAL(dp), INTENT(INOUT)                            :: accumulator(:, :)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: buffer_sel(:), buffer_w(:), buffer_out(:)
      INTEGER(int_8), INTENT(INOUT)                      :: n_used, n_products

      INTEGER                                            :: ik, il, ir, iseg, it, iv, k, mu, n_a, &
                                                            n_ab, n_cols, n_keep, n_rows, n_sel, &
                                                            n_valid, nu, row
      INTEGER, ALLOCATABLE                               :: keep(:), row_l(:), rows_a(:), rows_b(:), &
                                                            sel_rows(:), valid_cell(:), &
                                                            valid_seg(:)
      LOGICAL                                            :: used
      REAL(dp)                                           :: block_max
      REAL(dp), ALLOCATABLE                              :: pair_max(:)

      n_a = pair%n_a
      n_ab = n_a*pair%n_b
      n_cols = last_col - first_col + 1
      CALL segments_needing_block(pair, c, block_cell, cell_of_T, cell_shift, ri_centers, cutoffs, ri_bound, &
                                  threshold, n_valid, valid_seg, valid_cell)
      ALLOCATE (pair_max(n_ab), keep(n_ab))
      CALL ensure_buffer_size(buffer_sel, n_ab*n_cols)
      DO k = 1, n_ab
         pair_max(k) = MAXVAL(ABS(integrals(k:(n_cols - 1)*n_ab + k:n_ab)))
      END DO
      block_max = MAXVAL(pair_max)
      used = .FALSE.
      DO iv = 1, n_valid
         iseg = valid_seg(iv)
         it = valid_cell(iv)
         IF (pair%seg_amax(iseg)*block_max <= threshold) CYCLE
         ! Products φ_μ φ_ν whose bound exceeds the threshold.
         n_keep = 0
         DO nu = 1, pair%n_b
            DO mu = 1, n_a
               k = mu + (nu - 1)*n_a
               IF (pair%fac*pair%amax_a(mu, iseg)*pair%amax_b(nu, iseg)*pair_max(k) <= threshold) CYCLE
               n_keep = n_keep + 1
               keep(n_keep) = k
            END DO
         END DO
         IF (n_keep == 0) CYCLE
         ! Rows inside the d_lP cutoff of the RI center, as in the numerical quadrature.
         CALL shared_grid_rows(images, row_of, pair%seg_image(1, iseg), pair%seg_image(2, iseg), n_rows, &
                               row_l, rows_a, rows_b)
         ALLOCATE (sel_rows(n_rows))
         n_sel = 0
         DO k = 1, n_rows
            IF (SUM((grid_points(:, row_l(k)) - ri_centers(:, c, it))**2) > cutoffs(c)**2) CYCLE
            n_sel = n_sel + 1
            sel_rows(n_sel) = k
         END DO
         IF (n_sel > 0) THEN
            ! w(keep, row) = fac φ_μ(r_l) φ_ν(r_l), out(P, row) = Σ_keep integrals w
            CALL ensure_buffer_size(buffer_w, n_keep*n_sel)
            CALL ensure_buffer_size(buffer_out, n_cols*n_sel)
            ASSOCIATE (primary_a => images(pair%seg_image(1, iseg))%primary, &
                       primary_b => images(pair%seg_image(2, iseg))%primary)
               DO ir = 1, n_sel
                  row = sel_rows(ir)
                  DO ik = 1, n_keep
                     nu = (keep(ik) - 1)/n_a + 1
                     mu = keep(ik) - (nu - 1)*n_a
                     buffer_w(ik + (ir - 1)*n_keep) = pair%fac*primary_a(rows_a(row), mu)* &
                                                      primary_b(rows_b(row), nu)
                  END DO
               END DO
            END ASSOCIATE
            IF (n_keep == n_ab) THEN
               CALL dgemm('T', 'N', n_cols, n_sel, n_ab, 1.0_dp, integrals, n_ab, buffer_w, n_ab, 0.0_dp, &
                          buffer_out, n_cols)
            ELSE
               DO ik = 1, n_keep
                  buffer_sel(ik:(n_cols - 1)*n_ab + ik:n_ab) = integrals(keep(ik):(n_cols - 1)*n_ab + keep(ik):n_ab)
               END DO
               CALL dgemm('T', 'N', n_cols, n_sel, n_keep, 1.0_dp, buffer_sel, n_ab, buffer_w, n_keep, 0.0_dp, &
                          buffer_out, n_cols)
            END IF
            n_products = n_products + INT(n_cols, int_8)*INT(n_sel, int_8)*INT(n_keep, int_8)
            used = .TRUE.
            row = (it - 1)*(SIZE(accumulator, 1)/SIZE(locks)) + first_col - col_offset - 1
            CALL omp_set_lock(locks(it))
            DO ir = 1, n_sel
               il = row_l(sel_rows(ir))
               accumulator(row + 1:row + n_cols, il) = accumulator(row + 1:row + n_cols, il) + &
                                                       buffer_out((ir - 1)*n_cols + 1:ir*n_cols)
            END DO
            CALL omp_unset_lock(locks(it))
         END IF
         DEALLOCATE (sel_rows)
      END DO
      IF (used) n_used = n_used + 1

   END SUBROUTINE contract_block

! **************************************************************************************************
!> \brief Sends the buffered integral slices to the team ranks holding their columns.
!> \param team_comm ...
!> \param n_send ...
!> \param send_head Per slice: team rank, pair, c, T', first and last column, data offset, length.
!> \param send_data ...
!> \param recv_head Per received slice: pair, c, T', first and last column, data offset.
!> \param recv_data ...
! **************************************************************************************************
   SUBROUTINE exchange_slices(team_comm, n_send, send_head, send_data, recv_head, recv_data)
      TYPE(mp_comm_type), INTENT(IN)                     :: team_comm
      INTEGER, INTENT(IN)                                :: n_send, send_head(:, :)
      REAL(dp), INTENT(IN)                               :: send_data(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: recv_head(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: recv_data(:)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'exchange_slices'

      INTEGER                                            :: handle, i, irank, k, n_recv
      INTEGER, ALLOCATABLE :: data_count(:), data_displ(:), head_count(:), head_displ(:), &
         next_data(:), next_head(:), received_head(:), recv_data_count(:), recv_data_displ(:), &
         recv_head_count(:), recv_head_displ(:), sorted_head(:)
      REAL(dp), ALLOCATABLE                              :: sorted_data(:)

      CALL timeset(routineN, handle)

      ALLOCATE (head_count(team_comm%num_pe), data_count(team_comm%num_pe), &
                recv_head_count(team_comm%num_pe), recv_data_count(team_comm%num_pe))
      head_count = 0
      data_count = 0
      DO i = 1, n_send
         irank = send_head(1, i) + 1
         head_count(irank) = head_count(irank) + 8
         data_count(irank) = data_count(irank) + send_head(10, i)
      END DO
      CALL team_comm%alltoall(head_count, recv_head_count, 1)
      CALL team_comm%alltoall(data_count, recv_data_count, 1)
      ALLOCATE (head_displ(team_comm%num_pe), data_displ(team_comm%num_pe), &
                recv_head_displ(team_comm%num_pe), recv_data_displ(team_comm%num_pe))
      head_displ(1) = 0
      data_displ(1) = 0
      recv_head_displ(1) = 0
      recv_data_displ(1) = 0
      DO irank = 2, team_comm%num_pe
         head_displ(irank) = head_displ(irank - 1) + head_count(irank - 1)
         data_displ(irank) = data_displ(irank - 1) + data_count(irank - 1)
         recv_head_displ(irank) = recv_head_displ(irank - 1) + recv_head_count(irank - 1)
         recv_data_displ(irank) = recv_data_displ(irank - 1) + recv_data_count(irank - 1)
      END DO

      ! Sort by team rank: pair, c, T', first and last column, length.
      ALLOCATE (sorted_head(SUM(head_count)), sorted_data(SUM(data_count)), &
                next_head(team_comm%num_pe), next_data(team_comm%num_pe))
      next_head(:) = head_displ(:)
      next_data(:) = data_displ(:)
      DO i = 1, n_send
         irank = send_head(1, i) + 1
         sorted_head(next_head(irank) + 1:next_head(irank) + 8) = [send_head(2:8, i), send_head(10, i)]
         sorted_data(next_data(irank) + 1:next_data(irank) + send_head(10, i)) = &
            send_data(send_head(9, i) + 1:send_head(9, i) + send_head(10, i))
         next_head(irank) = next_head(irank) + 8
         next_data(irank) = next_data(irank) + send_head(10, i)
      END DO
      ALLOCATE (received_head(SUM(recv_head_count)), recv_data(SUM(recv_data_count)))
      CALL team_comm%alltoall(sorted_head, head_count, head_displ, received_head, recv_head_count, &
                              recv_head_displ)
      CALL team_comm%alltoall(sorted_data, data_count, data_displ, recv_data, recv_data_count, &
                              recv_data_displ)

      n_recv = SIZE(received_head)/8
      ALLOCATE (recv_head(8, n_recv))
      k = 0
      DO i = 1, n_recv
         recv_head(1:7, i) = received_head(8*i - 7:8*i - 1)
         recv_head(8, i) = k
         k = k + received_head(8*i)
      END DO

      CALL timestop(handle)

   END SUBROUTINE exchange_slices

! **************************************************************************************************
!> \brief Whether the analytic d_lP (Eq. 7) is estimated to be faster than the uniform quadrature
!> (Eq. 10). Per RI-RS point r_l, n_AO(l) AO functions of atom images are nonzero (radius of each
!> basis set from its smallest exponent) and n_PT(l) RI functions P of cells T lie within the
!> d_lP cutoff. Quadrature: n_PT(l) products for each of the (4 pi/3) padding^3/h^3 integration
!> points around r_l (spacing h from the largest, padding from the smallest AO exponent). Analytic:
!> n_AO(l) (n_AO(l) + 1) n_PT(l) AO pairs times RI functions. This is a cost measure, not the
!> number of flops: screening skips most of these products, so the contraction executes 10-80x
!> fewer flops. time_per_product and time_per_ao_pair_ri are measured core times per quadrature
!> product and per AO pair times RI function of the whole analytic d_lP (integrals, rotation,
!> screened contraction and exchange between ranks), fitted to Si, GaAs, AlAs, MgO, WS2 and TiO2
!> with SZ2 and TZ2 basis sets on 1 and 4 nodes.
!> \param bs_env ...
!> \param grid_points ...
!> \param cells_R_d_lP Cells T of d_lP.
!> \return ...
! **************************************************************************************************
   FUNCTION analytic_d_lP_is_faster(bs_env, grid_points, cells_R_d_lP) RESULT(analytic)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      LOGICAL                                            :: analytic

      CHARACTER(LEN=*), PARAMETER :: routineN = 'analytic_d_lP_is_faster'
      REAL(dp), PARAMETER :: time_per_ao_pair_ri = 6.0E-12_dp, time_per_product = 1.24E-10_dp

      INTEGER                                            :: a, c, first_l, handle, il, ir, iset, it, &
                                                            last_l, n_ao, n_PT
      REAL(dp)                                           :: ao_pairs_ri, cores, distance, padding, &
                                                            products, set_radius, spacing, &
                                                            time_analytic, time_quadrature
      REAL(dp), ALLOCATABLE                              :: cutoffs(:), positions(:, :)

      CALL timeset(routineN, handle)

      CALL d_lP_quadrature_resolution(bs_env, spacing, padding)
      CALL compute_d_lP_cutoffs(bs_env, cutoffs)
      ALLOCATE (positions(3, bs_env%n_atom))
      DO a = 1, bs_env%n_atom
         positions(:, a) = pbc(bs_env%ri_rs%particle_set(a)%r, bs_env%ri_rs%cell)
      END DO

      CALL get_local_l_range(SIZE(grid_points, 2), bs_env%para_env%num_pe, bs_env%para_env%mepos, &
                             first_l, last_l)
      products = 0.0_dp
      ao_pairs_ri = 0.0_dp
      DO il = first_l, last_l
         n_ao = 0
         DO ir = 1, SIZE(bs_env%index_to_cell_3c, 2)
            DO a = 1, bs_env%n_atom
               distance = NORM2(grid_points(:, il) - positions(:, a) - &
                                MATMUL(bs_env%hmat, REAL(bs_env%index_to_cell_3c(:, ir), dp)))
               ASSOCIATE (basis => bs_env%basis_set_AO(bs_env%ri_rs%particle_set(a)%atomic_kind%kind_number)%gto_basis_set)
                  DO iset = 1, basis%nset
                     set_radius = SQRT(LOG(1.0_dp/d_lP_integration_accuracy(bs_env))/ &
                                       MINVAL(basis%zet(1:basis%npgf(iset), iset)))
                     IF (distance < set_radius) n_ao = n_ao + basis%nsgf_set(iset)
                  END DO
               END ASSOCIATE
            END DO
         END DO
         n_PT = 0
         DO it = 1, SIZE(cells_R_d_lP, 2)
            DO c = 1, bs_env%n_atom
               IF (NORM2(grid_points(:, il) - positions(:, c) - &
                         MATMUL(bs_env%hmat, REAL(cells_R_d_lP(:, it), dp))) < cutoffs(c)) THEN
                  n_PT = n_PT + bs_env%i_RI_end_from_atom(c) - bs_env%i_RI_start_from_atom(c) + 1
               END IF
            END DO
         END DO
         products = products + REAL(n_PT, dp)*4.0_dp*pi/3.0_dp*(padding/spacing)**3
         ao_pairs_ri = ao_pairs_ri + REAL(n_ao, dp)*REAL(n_ao + 1, dp)*REAL(n_PT, dp)
      END DO
      CALL bs_env%para_env%sum(products)
      CALL bs_env%para_env%sum(ao_pairs_ri)

      cores = REAL(bs_env%para_env%num_pe*omp_get_max_threads(), dp)
      time_quadrature = time_per_product*products/cores
      time_analytic = time_per_ao_pair_ri*ao_pairs_ri/cores
      analytic = time_analytic < time_quadrature

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') 'd_lP estimated quadrature products', products
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') 'd_lP estimated analytic AO pairs times RI functions', &
            ao_pairs_ri
         WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') 'd_lP estimated time of the quadrature', &
            time_quadrature, ' s'
         WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') 'd_lP estimated time of the analytic integrals', &
            time_analytic, ' s'
      END IF

      CALL timestop(handle)

   END FUNCTION analytic_d_lP_is_faster

! **************************************************************************************************
!> \brief Space group of the crystal: point-group rotations with fractional translations from K290
!> (as for the k-point symmetry), combined with the pure lattice translations of a supercell.
!> For every operation the atom permutation, the lattice cells of the mapped atoms, the integer
!> rotation of lattice vectors and the rotation matrices of the real solid harmonics.
!> \param bs_env ...
!> \param cells_R_d_lP Cells T of d_lP, bound the cell components of the block codes.
!> \param group ...
! **************************************************************************************************
   SUBROUTINE crystal_space_group(bs_env, cells_R_d_lP, group)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      TYPE(space_group_type), INTENT(OUT)                :: group

      CHARACTER(LEN=*), PARAMETER :: routineN = 'crystal_space_group'

      INTEGER                                            :: a, b, cell_max, g, handle, ikind, it, k, &
                                                            l, lmax, n_atom, n_rot, n_trans, op
      INTEGER, ALLOCATABLE                               :: kinds(:)
      LOGICAL                                            :: from_spglib
      REAL(dp)                                           :: cartesian(3, 3), frac(3, 3), shift(3)
      REAL(dp), ALLOCATABLE                              :: rotations(:, :, :), scaled(:, :), &
                                                            shifts(:, :), translations(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(csym_type)                                    :: csym
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot

      CALL timeset(routineN, handle)

      n_atom = bs_env%n_atom
      cell => bs_env%ri_rs%cell
      ALLOCATE (scaled(3, n_atom), kinds(n_atom))
      ! Atoms at pbc(r), as in the three-center integrals.
      DO a = 1, n_atom
         CALL real_to_scaled(scaled(:, a), pbc(bs_env%ri_rs%particle_set(a)%r, cell), cell)
         kinds(a) = bs_env%ri_rs%particle_set(a)%atomic_kind%kind_number
      END DO
      ! Space group from SPGLIB if linked (all operations, also the pure lattice translations of a
      ! supercell), otherwise point-group rotations with fractional translations from K290.
      CALL crys_sym_gen(csym, scaled, kinds, cell%hmat, delta=eps_symmetry, &
                        use_spglib=spg_get_major_version() > 0)
      IF (csym%symlib) THEN
         n_rot = csym%n_operations
         ALLOCATE (rotations(3, 3, n_rot), shifts(3, n_rot))
         DO g = 1, n_rot
            rotations(:, :, g) = MATMUL(cell%hmat, MATMUL(TRANSPOSE(REAL(csym%rotations(:, :, g), dp)), &
                                                          cell%h_inv))
         END DO
         shifts(:, :) = csym%translations(:, 1:n_rot)
      ELSE
         CALL kpoint_gen(csym, [1, 1, 1], symm=.TRUE., shift=[0.0_dp, 0.0_dp, 0.0_dp])
         n_rot = MAX(csym%nrtot, 1)
         ALLOCATE (rotations(3, 3, n_rot), shifts(3, n_rot))
         rotations = 0.0_dp
         shifts = 0.0_dp
         IF (csym%nrtot > 0) THEN
            rotations(:, :, :) = csym%rt(:, :, 1:n_rot)
            shifts(:, :) = csym%vt(:, 1:n_rot)
         ELSE
            DO k = 1, 3
               rotations(k, k, 1) = 1.0_dp
            END DO
         END IF
      END IF
      from_spglib = csym%symlib
      CALL release_csym_type(csym)

      ! Pure lattice translations that map the crystal onto itself (more than one in a supercell).
      ALLOCATE (translations(3, n_atom))
      n_trans = 0
      frac = 0.0_dp
      DO k = 1, 3
         frac(k, k) = 1.0_dp
      END DO
      DO b = 1, n_atom
         IF (kinds(b) /= kinds(1) .OR. (from_spglib .AND. b > 1)) CYCLE
         shift = scaled(:, b) - scaled(:, 1)
         shift = shift - ANINT(shift)
         IF (.NOT. maps_crystal(scaled, kinds, frac, shift)) CYCLE
         n_trans = n_trans + 1
         translations(:, n_trans) = shift
      END DO

      lmax = 0
      DO ikind = 1, SIZE(bs_env%basis_set_AO)
         lmax = MAX(lmax, MAXVAL(bs_env%basis_set_AO(ikind)%gto_basis_set%lmax), &
                    MAXVAL(bs_env%basis_set_RI(ikind)%gto_basis_set%lmax))
      END DO
      group%n_atom = n_atom
      group%n_op = n_rot*n_trans
      ALLOCATE (group%atom_map(n_atom, group%n_op), group%atom_cell(3, n_atom, group%n_op), &
                group%frac_rot(3, 3, group%n_op), group%sph_rot(nso(lmax), nso(lmax), 0:lmax, group%n_op))
      group%sph_rot = 0.0_dp
      NULLIFY (orbrot)
      op = 0
      DO g = 1, n_rot
         frac = MATMUL(cell%h_inv, MATMUL(rotations(:, :, g), cell%hmat))
         CPASSERT(MAXVAL(ABS(frac - ANINT(frac))) < SQRT(eps_symmetry))
         frac = ANINT(frac)
         ! Cartesian rotation of the lattice operation, exact for the integer matrix frac.
         cartesian = MATMUL(cell%hmat, MATMUL(frac, cell%h_inv))
         CALL calculate_rotmat(orbrot, cartesian, lmax)
         DO it = 1, n_trans
            op = op + 1
            group%frac_rot(:, :, op) = NINT(frac)
            DO l = 0, lmax
               group%sph_rot(1:nso(l), 1:nso(l), l, op) = orbrot(l)%mat(1:nso(l), 1:nso(l))
            END DO
            DO a = 1, n_atom
               CALL image_of_atom(scaled, kinds, a, frac, shifts(:, g) + translations(:, it), &
                                  group%atom_map(a, op), group%atom_cell(:, a, op))
               CPASSERT(group%atom_map(a, op) > 0)
            END DO
            IF (ALL(group%frac_rot(:, :, op) == RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])) .AND. &
                ALL(group%atom_cell(:, :, op) == 0) .AND. &
                ALL(group%atom_map(:, op) == [(a, a=1, n_atom)])) group%identity = op
         END DO
      END DO
      CALL release_rotmat(orbrot)
      CPASSERT(group%identity > 0)

      ! Cell components of transformed blocks: |S'|, |T'| of a pair and its exchange, rotated.
      cell_max = MAX(MAXVAL(ABS(cells_R_d_lP)), MAXVAL(ABS(bs_env%index_to_cell_3c)))
      group%cell_range = 0
      DO op = 1, group%n_op
         group%cell_range = MAX(group%cell_range, 4*cell_max*MAXVAL(SUM(ABS(group%frac_rot(:, :, op)), DIM=2)) + &
                                2*MAXVAL(ABS(group%atom_cell(:, :, op))))
      END DO
      IF (3.0_dp*LOG(REAL(n_atom, dp)) + 6.0_dp*LOG(REAL(2*group%cell_range + 1, dp)) >= &
          LOG(REAL(HUGE(0_int_8), dp))) THEN
         CPABORT("Block codes of the analytic d_lP exceed 64-bit integers")
      END IF

      CALL timestop(handle)

   END SUBROUTINE crystal_space_group

! **************************************************************************************************
!> \brief Atom b and lattice cell L with frac s_a + shift = s_b + L for an atom b of the kind of a.
!> \param scaled Fractional atom positions.
!> \param kinds ...
!> \param a ...
!> \param frac Rotation in fractional coordinates.
!> \param shift Fractional translation.
!> \param b Zero if no atom matches.
!> \param cell ...
! **************************************************************************************************
   PURE SUBROUTINE image_of_atom(scaled, kinds, a, frac, shift, b, cell)
      REAL(dp), INTENT(IN)                               :: scaled(:, :)
      INTEGER, INTENT(IN)                                :: kinds(:), a
      REAL(dp), INTENT(IN)                               :: frac(3, 3), shift(3)
      INTEGER, INTENT(OUT)                               :: b, cell(3)

      INTEGER                                            :: j
      REAL(dp)                                           :: image(3)

      image = MATMUL(frac, scaled(:, a)) + shift
      b = 0
      cell = 0
      DO j = 1, SIZE(kinds)
         IF (kinds(j) /= kinds(a)) CYCLE
         IF (MAXVAL(ABS(image - scaled(:, j) - ANINT(image - scaled(:, j)))) > SQRT(eps_symmetry)) CYCLE
         b = j
         cell = NINT(image - scaled(:, j))
         RETURN
      END DO

   END SUBROUTINE image_of_atom

! **************************************************************************************************
!> \brief Whether s -> frac s + shift maps every atom onto an atom of the same kind.
!> \param scaled ...
!> \param kinds ...
!> \param frac ...
!> \param shift ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION maps_crystal(scaled, kinds, frac, shift) RESULT(maps)
      REAL(dp), INTENT(IN)                               :: scaled(:, :)
      INTEGER, INTENT(IN)                                :: kinds(:)
      REAL(dp), INTENT(IN)                               :: frac(3, 3), shift(3)
      LOGICAL                                            :: maps

      INTEGER                                            :: a, b, cell(3)

      maps = .TRUE.
      DO a = 1, SIZE(kinds)
         CALL image_of_atom(scaled, kinds, a, frac, shift, b, cell)
         IF (b == 0) THEN
            maps = .FALSE.
            RETURN
         END IF
      END DO

   END FUNCTION maps_crystal

! **************************************************************************************************
!> \brief Code of the symmetry-unique block of (a 0, b S'|c T') and the operation g with
!> unique block = g(block), or g(exchanged block) for op = -g.
!> \param group ...
!> \param block (a, b, S', c, T')
!> \param code ...
!> \param op ...
! **************************************************************************************************
   PURE SUBROUTINE unique_block(group, block, code, op)
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER, INTENT(IN)                                :: block(9)
      INTEGER(int_8), INTENT(OUT)                        :: code
      INTEGER, INTENT(OUT)                               :: op

      INTEGER                                            :: g, sgn, source(9)
      INTEGER(int_8)                                     :: atoms_ab, best_atoms_ab, trial

      code = HUGE(code)
      best_atoms_ab = HUGE(best_atoms_ab)
      op = 0
      DO sgn = 1, -1, -2
         source = block
         IF (sgn < 0) source = [block(2), block(1), -block(3:5), block(6), block(7:9) - block(3:5)]
         DO g = 1, group%n_op
            ! The atoms (a, b) are the leading digits of the block code: an operation that maps
            ! them onto a larger pair than the best one so far gives a larger code.
            atoms_ab = INT(group%atom_map(source(1), g) - 1, int_8)*group%n_atom + group%atom_map(source(2), g) - 1
            IF (atoms_ab > best_atoms_ab) CYCLE
            trial = block_code(group, transformed_block(group, g, source))
            IF (trial < code) THEN
               code = trial
               best_atoms_ab = atoms_ab
               op = sgn*g
            END IF
         END DO
      END DO

   END SUBROUTINE unique_block

! **************************************************************************************************
!> \brief Block g(a 0, b S'|c T') shifted back so that its first atom lies in cell 0.
!> \param group ...
!> \param g ...
!> \param block ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION transformed_block(group, g, block) RESULT(image)
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER, INTENT(IN)                                :: g, block(9)
      INTEGER                                            :: image(9)

      image(1) = group%atom_map(block(1), g)
      image(2) = group%atom_map(block(2), g)
      image(6) = group%atom_map(block(6), g)
      image(3:5) = group%atom_cell(:, block(2), g) + MATMUL(group%frac_rot(:, :, g), block(3:5)) - &
                   group%atom_cell(:, block(1), g)
      image(7:9) = group%atom_cell(:, block(6), g) + MATMUL(group%frac_rot(:, :, g), block(7:9)) - &
                   group%atom_cell(:, block(1), g)

   END FUNCTION transformed_block

! **************************************************************************************************
!> \brief One integer for the block (a, b, S', c, T'); the inverse is decoded_block.
!> \param group ...
!> \param block ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION block_code(group, block) RESULT(code)
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER, INTENT(IN)                                :: block(9)
      INTEGER(int_8)                                     :: code

      INTEGER                                            :: i

      code = INT(block(1) - 1, int_8)*group%n_atom + block(2) - 1
      code = code*group%n_atom + block(6) - 1
      DO i = 3, 9
         IF (i == 6) CYCLE
         code = code*(2*group%cell_range + 1) + block(i) + group%cell_range
      END DO

   END FUNCTION block_code

! **************************************************************************************************
!> \brief Block (a, b, S', c, T') of a code from block_code.
!> \param group ...
!> \param code ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION decoded_block(group, code) RESULT(block)
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER(int_8), INTENT(IN)                         :: code
      INTEGER                                            :: block(9)

      INTEGER                                            :: i
      INTEGER(int_8)                                     :: rest

      rest = code
      DO i = 9, 3, -1
         IF (i == 6) CYCLE
         block(i) = INT(MODULO(rest, INT(2*group%cell_range + 1, int_8))) - group%cell_range
         rest = rest/(2*group%cell_range + 1)
      END DO
      block(6) = INT(MODULO(rest, INT(group%n_atom, int_8))) + 1
      rest = rest/group%n_atom
      block(2) = INT(MODULO(rest, INT(group%n_atom, int_8))) + 1
      block(1) = INT(rest/group%n_atom) + 1

   END FUNCTION decoded_block

! **************************************************************************************************
!> \brief Rank that computes and keeps the integrals of a symmetry-unique block (scrambled code).
!> \param code ...
!> \param n_ranks ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION owner_of_block(code, n_ranks) RESULT(owner)
      INTEGER(int_8), INTENT(IN)                         :: code
      INTEGER, INTENT(IN)                                :: n_ranks
      INTEGER                                            :: owner

      INTEGER(int_8), PARAMETER                          :: modulus = 2147483647_int_8, &
                                                            multiplier = 48271_int_8

      owner = INT(MODULO(MODULO(MODULO(code, modulus)*multiplier, modulus), INT(n_ranks, int_8)))

   END FUNCTION owner_of_block

! **************************************************************************************************
!> \brief Integrals of the block (a 0, b S'|c T') in the reference RI basis from those of its
!> symmetry-unique block: rotation of all three function indices with the inverse operation,
!> then exchange of the two AO indices for op < 0.
!> \param bs_env ...
!> \param group ...
!> \param code Code of the symmetry-unique block.
!> \param op Operation from unique_block.
!> \param source Integrals of the symmetry-unique block, indexed (μ, ν, p).
!> \param rotated Integrals of the block.
! **************************************************************************************************
   SUBROUTINE rotate_block(bs_env, group, code, op, source, rotated)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(space_group_type), INTENT(IN)                 :: group
      INTEGER(int_8), INTENT(IN)                         :: code
      INTEGER, INTENT(IN)                                :: op
      REAL(dp), INTENT(IN)                               :: source(:)
      REAL(dp), INTENT(OUT)                              :: rotated(:)

      INTEGER                                            :: block(9), g, mu, n_a, n_b, n_p, nu, p
      REAL(dp), ALLOCATABLE                              :: first(:), second(:)

      g = ABS(op)
      block = decoded_block(group, code)
      n_a = nsgf_atom(bs_env, block(1))
      n_b = nsgf_atom(bs_env, block(2))
      n_p = nsgf_ri(bs_env, block(6))
      IF (op == group%identity) THEN
         rotated(:) = source(:)
         RETURN
      END IF
      ALLOCATE (first(SIZE(source)), second(SIZE(source)))
      IF (ALL(group%frac_rot(:, :, g) == RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))) THEN
         ! Pure lattice translation: the functions keep their orientation.
         second(:) = source(:)
      ELSE
         ASSOCIATE (kinds => bs_env%ri_rs%particle_set)
            CALL rotate_solid_harmonic_index(bs_env%basis_set_AO(kinds(block(1))%atomic_kind%kind_number)%gto_basis_set, &
                                             group%sph_rot(:, :, :, g), 1, n_a, n_b*n_p, source, first)
            CALL rotate_solid_harmonic_index(bs_env%basis_set_AO(kinds(block(2))%atomic_kind%kind_number)%gto_basis_set, &
                                             group%sph_rot(:, :, :, g), n_a, n_b, n_p, first, second)
            CALL rotate_solid_harmonic_index(bs_env%basis_set_RI(kinds(block(6))%atomic_kind%kind_number)%gto_basis_set, &
                                             group%sph_rot(:, :, :, g), n_a*n_b, n_p, 1, second, first)
         END ASSOCIATE
         second(:) = first(:)
      END IF
      IF (op > 0) THEN
         rotated(:) = second(:)
      ELSE
         DO p = 1, n_p
            DO mu = 1, n_a
               DO nu = 1, n_b
                  rotated(nu + (mu - 1)*n_b + (p - 1)*n_a*n_b) = second(mu + (nu - 1)*n_a + (p - 1)*n_a*n_b)
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE rotate_block

! **************************************************************************************************
!> \brief rotated(i, μ, k) = Σ_μ' D(μ', μ) source(i, μ', k) with the block-diagonal rotation
!> D of the real solid harmonics of the shells of a basis (calculate_rotmat of the operation).
!> \param basis ...
!> \param sph_rot Rotation matrices per angular momentum.
!> \param n_before ...
!> \param n_index ...
!> \param n_after ...
!> \param source ...
!> \param rotated ...
! **************************************************************************************************
   SUBROUTINE rotate_solid_harmonic_index(basis, sph_rot, n_before, n_index, n_after, source, rotated)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: sph_rot(:, :, 0:)
      INTEGER, INTENT(IN)                                :: n_before, n_index, n_after
      REAL(dp), INTENT(IN) :: source(n_before, n_index, n_after)
      REAL(dp), INTENT(OUT) :: rotated(n_before, n_index, n_after)

      REAL(dp), PARAMETER                                :: eps_rotation = 1.0E-14_dp

      INTEGER                                            :: f, iset, ishell, k, l, m, m_source
      REAL(dp)                                           :: d

      ! Rotations of crystal symmetry are mostly signed permutations: zero entries are skipped.
      rotated = 0.0_dp
      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            f = basis%first_sgf(ishell, iset) - 1
            DO m = 1, nso(l)
               DO m_source = 1, nso(l)
                  d = sph_rot(m_source, m, l)
                  IF (ABS(d) < eps_rotation) CYCLE
                  DO k = 1, n_after
                     rotated(:, f + m, k) = rotated(:, f + m, k) + d*source(:, f + m_source, k)
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE rotate_solid_harmonic_index

! **************************************************************************************************
!> \brief Lookup tables: image index of (atom, cell), packed AO row of (l, image), index of cell T.
!> \param bs_env ...
!> \param images ...
!> \param cells_R_d_lP ...
!> \param image_of Image of atom a in cell R, image_of(a, R1, R2, R3), zero if absent.
!> \param row_of Row of r_l in images(i)%primary, row_of(l, i), zero if the AO values vanish.
!> \param cell_of_T Index of T in cells_R_d_lP, zero if absent.
! **************************************************************************************************
   SUBROUTINE index_images_and_cells(bs_env, images, cells_R_d_lP, image_of, row_of, cell_of_T)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: image_of(:, :, :, :), row_of(:, :), &
                                                            cell_of_T(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'index_images_and_cells'

      INTEGER                                            :: handle, i, il, it, low(3), up(3)

      CALL timeset(routineN, handle)

      low = MINVAL(bs_env%index_to_cell_3c, DIM=2)
      up = MAXVAL(bs_env%index_to_cell_3c, DIM=2)
      ALLOCATE (image_of(bs_env%n_atom, low(1):up(1), low(2):up(2), low(3):up(3)))
      ALLOCATE (row_of(bs_env%ri_rs%n_grid_points, SIZE(images)))
      image_of = 0
      row_of = 0
      DO i = 1, SIZE(images)
         ! prepare_ao_images leaves the slots of vanishing images empty (atom 0).
         IF (images(i)%atom == 0) CYCLE
         image_of(images(i)%atom, images(i)%cell(1), images(i)%cell(2), images(i)%cell(3)) = i
         IF (ALLOCATED(images(i)%primary_l)) THEN
            DO il = 1, SIZE(images(i)%primary_l)
               row_of(images(i)%primary_l(il), i) = il
            END DO
         ELSE
            DO il = 1, SIZE(images(i)%primary, 1)
               IF (ANY(images(i)%primary(il, :) /= 0.0_dp)) row_of(il, i) = il
            END DO
         END IF
      END DO
      low = MINVAL(cells_R_d_lP, DIM=2)
      up = MAXVAL(cells_R_d_lP, DIM=2)
      ALLOCATE (cell_of_T(low(1):up(1), low(2):up(2), low(3):up(3)))
      cell_of_T = 0
      DO it = 1, SIZE(cells_R_d_lP, 2)
         cell_of_T(cells_R_d_lP(1, it), cells_R_d_lP(2, it), cells_R_d_lP(3, it)) = it
      END DO

      CALL timestop(handle)

   END SUBROUTINE index_images_and_cells

! **************************************************************************************************
!> \brief AUTO_RI contraction U_Pp of every atom, identity without AUTO_RI.
!> \param bs_env ...
!> \param U_atom U_atom(reference function, contracted function, atom).
! **************************************************************************************************
   SUBROUTINE contracted_ri_coefficients(bs_env, U_atom)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: U_atom(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contracted_ri_coefficients'

      INTEGER                                            :: first_U, handle, iatom, ip, max_ref, np, &
                                                            nref

      CALL timeset(routineN, handle)

      max_ref = 0
      DO iatom = 1, bs_env%n_atom
         max_ref = MAX(max_ref, nsgf_ri(bs_env, iatom))
      END DO
      ALLOCATE (U_atom(max_ref, MAXVAL(bs_env%i_RI_end_from_atom - bs_env%i_RI_start_from_atom + 1), &
                       bs_env%n_atom))
      U_atom = 0.0_dp
      DO iatom = 1, bs_env%n_atom
         nref = nsgf_ri(bs_env, iatom)
         np = bs_env%i_RI_end_from_atom(iatom) - bs_env%i_RI_start_from_atom(iatom) + 1
         IF (bs_env%auto_ri%enabled) THEN
            CPASSERT(bs_env%auto_ri%sizes_ref_RI(iatom) == nref)
            CPASSERT(bs_env%auto_ri%sizes_opt_RI(iatom) == np)
            first_U = bs_env%auto_ri%U_Pp_atom_offset(iatom)
            U_atom(1:nref, 1:np, iatom) = RESHAPE(bs_env%auto_ri%U_Pp_atom(first_U:first_U + nref*np - 1), [nref, np])
         ELSE
            CPASSERT(nref == np)
            DO ip = 1, np
               U_atom(ip, ip, iatom) = 1.0_dp
            END DO
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE contracted_ri_coefficients

! **************************************************************************************************
!> \brief Factors of the bound |(μ ν|P)| <= Kmax ||μ ν||_2 ||P||_2 <= min(|μ|_inf, |ν|_inf) ri_bound.
!> Cauchy-Schwarz with the positive erfc kernel, whose Fourier transform is at most
!> Kmax = pi/omega**2; the AOs and reference RI functions are normalized.
!> \param bs_env ...
!> \param U_atom ...
!> \param ao_max Bound of |φ_μ(r)| for every AO function, ao_max(μ, atom).
!> \param ri_bound Kmax max_p Σ_ref |U(ref, p)| of every atom.
! **************************************************************************************************
   SUBROUTINE integral_bounds(bs_env, U_atom, ao_max, ri_bound)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: U_atom(:, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: ao_max(:, :), ri_bound(:)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'integral_bounds'

      INTEGER                                            :: handle, iatom, ico, idir, ipgf, iset, &
                                                            isgf, n, power
      REAL(dp)                                           :: factor
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CALL timeset(routineN, handle)

      ALLOCATE (ao_max(MAXVAL(bs_env%i_ao_end_from_atom - bs_env%i_ao_start_from_atom + 1), bs_env%n_atom), &
                ri_bound(bs_env%n_atom))
      ao_max = 0.0_dp
      DO iatom = 1, bs_env%n_atom
         basis => bs_env%basis_set_AO(bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number)%gto_basis_set
         ! |x^n exp(-zeta r^2)| <= PRODUCT_d (n_d/(2 zeta e))^(n_d/2)
         DO iset = 1, basis%nset
            n = ncoset(basis%lmax(iset))
            DO ipgf = 1, basis%npgf(iset)
               DO ico = 1, n
                  factor = 1.0_dp
                  DO idir = 1, 3
                     power = indco(idir, ico)
                     IF (power > 0) factor = factor*(REAL(power, dp)/(2.0_dp*basis%zet(ipgf, iset)*EXP(1.0_dp)))** &
                                             (0.5_dp*REAL(power, dp))
                  END DO
                  DO isgf = basis%first_sgf(1, iset), basis%last_sgf(basis%nshell(iset), iset)
                     ao_max(isgf, iatom) = ao_max(isgf, iatom) + ABS(basis%sphi((ipgf - 1)*n + ico, isgf))*factor
                  END DO
               END DO
            END DO
         END DO
         ri_bound(iatom) = pi/bs_env%ri_metric%omega**2*MAXVAL(SUM(ABS(U_atom(:, :, iatom)), DIM=1))
      END DO

      CALL timestop(handle)

   END SUBROUTINE integral_bounds

! **************************************************************************************************
!> \brief Gaussian envelopes C exp(-e |r-A|^2) bounding every primitive of the AO and RI functions.
!> C is the largest sum of |sphi| times the polynomial bound over the functions of a primitive;
!> for the RI functions the AUTO_RI coefficients |U| are included.
!> \param bs_env ...
!> \param U_atom ...
!> \param ao_envelope (exponent e, prefactor C) of the AO primitives of all atoms.
!> \param ao_first First AO primitive of each atom, ao_first(n_atom + 1) = end + 1.
!> \param ri_envelope (exponent e, prefactor C) of the RI primitives of all atoms.
!> \param ri_first First RI primitive of each atom.
! **************************************************************************************************
   SUBROUTINE primitive_envelopes(bs_env, U_atom, ao_envelope, ao_first, ri_envelope, ri_first)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: U_atom(:, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: ao_envelope(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: ao_first(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: ri_envelope(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: ri_first(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'primitive_envelopes'

      INTEGER                                            :: first_sgf, handle, iatom, ikind, ilist, &
                                                            ipgf, iset, isgf, jatom, last_sgf, n, &
                                                            n_prim
      REAL(dp), ALLOCATABLE                              :: function_bound(:), work(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CALL timeset(routineN, handle)

      DO ilist = 1, 2
         ALLOCATE (work(2, 0))
         n_prim = 0
         IF (ilist == 1) ALLOCATE (ao_first(bs_env%n_atom + 1))
         IF (ilist == 2) ALLOCATE (ri_first(bs_env%n_atom + 1))
         DO iatom = 1, bs_env%n_atom
            ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
            IF (ilist == 1) THEN
               basis => bs_env%basis_set_AO(ikind)%gto_basis_set
               ao_first(iatom) = n_prim + 1
            ELSE
               basis => bs_env%basis_set_RI(ikind)%gto_basis_set
               ri_first(iatom) = n_prim + 1
            END IF
            CALL grow_envelopes(work, n_prim + SUM(basis%npgf(1:basis%nset)))
            ALLOCATE (function_bound(basis%nsgf))
            DO iset = 1, basis%nset
               n = ncoset(basis%lmax(iset))
               first_sgf = basis%first_sgf(1, iset)
               last_sgf = basis%last_sgf(basis%nshell(iset), iset)
               DO ipgf = 1, basis%npgf(iset)
                  function_bound = 0.0_dp
                  DO isgf = first_sgf, last_sgf
                     function_bound(isgf) = polynomial_bound(basis%sphi((ipgf - 1)*n + 1:ipgf*n, isgf), &
                                                             basis%zet(ipgf, iset))
                  END DO
                  n_prim = n_prim + 1
                  work(1, n_prim) = envelope_fraction*basis%zet(ipgf, iset)
                  IF (ilist == 1) THEN
                     work(2, n_prim) = MAXVAL(function_bound(first_sgf:last_sgf))
                  ELSE
                     work(2, n_prim) = MAXVAL(MATMUL(function_bound, ABS(U_atom(1:basis%nsgf, :, iatom))))
                  END IF
               END DO
            END DO
            DEALLOCATE (function_bound)
         END DO
         IF (ilist == 1) THEN
            ao_first(bs_env%n_atom + 1) = n_prim + 1
            ALLOCATE (ao_envelope(2, n_prim), SOURCE=work(:, 1:n_prim))
         ELSE
            ri_first(bs_env%n_atom + 1) = n_prim + 1
            ! Largest prefactor over the atoms of a kind, so that the envelope bound is the same
            ! for symmetry-equivalent blocks.
            DO iatom = 1, bs_env%n_atom
               DO jatom = 1, bs_env%n_atom
                  IF (bs_env%ri_rs%particle_set(jatom)%atomic_kind%kind_number /= &
                      bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number) CYCLE
                  work(2, ri_first(iatom):ri_first(iatom + 1) - 1) = &
                     MAX(work(2, ri_first(iatom):ri_first(iatom + 1) - 1), work(2, ri_first(jatom):ri_first(jatom + 1) - 1))
               END DO
            END DO
            ALLOCATE (ri_envelope(2, n_prim), SOURCE=work(:, 1:n_prim))
         END IF
         DEALLOCATE (work)
      END DO

      CALL timestop(handle)

   END SUBROUTINE primitive_envelopes

! **************************************************************************************************
!> \brief Enlarge the envelope work list to at least n entries.
!> \param work ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE grow_envelopes(work, n)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: work(:, :)
      INTEGER, INTENT(IN)                                :: n

      REAL(dp), ALLOCATABLE                              :: grown(:, :)

      IF (n <= SIZE(work, 2)) RETURN
      ALLOCATE (grown(2, n))
      grown(:, 1:SIZE(work, 2)) = work
      CALL MOVE_ALLOC(grown, work)

   END SUBROUTINE grow_envelopes

! **************************************************************************************************
!> \brief Σ_ico |coefficient(ico)| max_x |x^i y^j z^k exp(-(1-f) zeta r^2)| of one primitive.
!> \param coefficients sphi rows of the Cartesian components of the primitive.
!> \param zeta ...
!> \return ...
! **************************************************************************************************
   FUNCTION polynomial_bound(coefficients, zeta) RESULT(bound)
      REAL(dp), INTENT(IN)                               :: coefficients(:), zeta
      REAL(dp)                                           :: bound

      INTEGER                                            :: ico, idir, power
      REAL(dp)                                           :: factor

      bound = 0.0_dp
      DO ico = 1, SIZE(coefficients)
         factor = 1.0_dp
         DO idir = 1, 3
            power = indco(idir, ico)
            IF (power > 0) factor = factor*(REAL(power, dp)/(2.0_dp*(1.0_dp - envelope_fraction)*zeta* &
                                                             EXP(1.0_dp)))**(0.5_dp*REAL(power, dp))
         END DO
         bound = bound + ABS(coefficients(ico))*factor
      END DO

   END FUNCTION polynomial_bound

! **************************************************************************************************
!> \brief Upper bound of |(μ ν|P)| for all functions of a block from Gaussian envelopes.
!> Two s Gaussians exp(-p|r-Q|^2), exp(-q|r-C|^2) interact through erfc(omega r)/r as
!> (pi^2/(p q))^(3/2) [erf(sqrt(rho) D) - erf(sqrt(rho') D)]/D, rho = pq/(p+q),
!> 1/rho' = 1/p + 1/q + 1/omega^2, D = |Q - C|.
!> \param pair_envelope (p, Q, prefactor) of the AO primitive pairs.
!> \param ri_envelope (q, prefactor) of the RI primitives.
!> \param center_C RI center.
!> \param omega ...
!> \return ...
! **************************************************************************************************
   FUNCTION envelope_bound(pair_envelope, ri_envelope, center_C, omega) RESULT(bound)
      REAL(dp), INTENT(IN)                               :: pair_envelope(:, :), ri_envelope(:, :), &
                                                            center_C(3), omega
      REAL(dp)                                           :: bound

      INTEGER                                            :: ipair, iri
      REAL(dp)                                           :: distance, p, q, rho, rho_omega

      bound = 0.0_dp
      DO ipair = 1, SIZE(pair_envelope, 2)
         p = pair_envelope(1, ipair)
         distance = NORM2(pair_envelope(2:4, ipair) - center_C)
         DO iri = 1, SIZE(ri_envelope, 2)
            q = ri_envelope(1, iri)
            rho = SQRT(p*q/(p + q))
            rho_omega = SQRT(1.0_dp/(1.0_dp/p + 1.0_dp/q + 1.0_dp/omega**2))
            IF (distance*rho < eps_rho_distance) THEN
               bound = bound + pair_envelope(5, ipair)*ri_envelope(2, iri)*(pi*pi/(p*q))**1.5_dp* &
                       2.0_dp/rootpi*(rho - rho_omega)
            ELSE
               ! erf(x) - erf(y) = erfc(y) - erfc(x) avoids cancellation at large distances.
               bound = bound + pair_envelope(5, ipair)*ri_envelope(2, iri)*(pi*pi/(p*q))**1.5_dp* &
                       (ERFC(rho_omega*distance) - ERFC(rho*distance))/distance
            END IF
         END DO
      END DO

   END FUNCTION envelope_bound

! **************************************************************************************************
!> \brief List the unique AO image pairs (a 0, b S').
!> Only a <= b is kept; for a = b only S' >= 0 (lexicographic). A pair needs two images whose
!> AO cutoff spheres overlap.
!> \param bs_env ...
!> \param images ...
!> \param image_of ...
!> \param key_list key_list(1:5, pair) = (a, b, S').
! **************************************************************************************************
   SUBROUTINE list_unique_ao_pairs(bs_env, images, image_of, key_list)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      INTEGER, INTENT(IN)                                :: image_of(:, :, :, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: key_list(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'list_unique_ao_pairs'

      INTEGER                                            :: a, atom_a, atom_b, handle, i, &
                                                            index_3d(3), n_keys, s(3), s1, s2, s3, &
                                                            shift(3), span(3), upper(3)
      INTEGER, ALLOCATABLE                               :: grown(:, :)
      REAL(dp), ALLOCATABLE                              :: ao_radius(:), position(:, :)

      CALL timeset(routineN, handle)

      ! AO cutoff sphere of every atom; CUTOFF_RADIUS_RL_AO replaces the automatic radius
      ALLOCATE (ao_radius(bs_env%n_atom), position(3, bs_env%n_atom))
      DO a = 1, bs_env%n_atom
         ao_radius(a) = bs_env%ri_rs%radius_ao_per_atom(a)
         IF (bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) ao_radius(a) = bs_env%ri_rs%cutoff_radius_ri_ao
         position(:, a) = pbc(bs_env%ri_rs%particle_set(a)%r, bs_env%ri_rs%cell)
      END DO

      ! image_of is indexed from 1; shift maps it back to the cell vector.
      upper = [SIZE(image_of, 2), SIZE(image_of, 3), SIZE(image_of, 4)]
      shift = MINVAL(bs_env%index_to_cell_3c, DIM=2) - 1
      span = upper - 1
      n_keys = 0
      ALLOCATE (key_list(5, 1024))
      DO atom_a = 1, bs_env%n_atom
         DO atom_b = atom_a, bs_env%n_atom
            DO s3 = -span(3), span(3)
               DO s2 = -span(2), span(2)
                  DO s1 = -span(1), span(1)
                     s = [s1, s2, s3]
                     ! (a 0, a S') and (a 0, a -S') are the same pair: keep S' >= 0 (lexicographic)
                     IF (atom_a == atom_b .AND. (s1 < 0 .OR. (s1 == 0 .AND. (s2 < 0 .OR. (s2 == 0 .AND. s3 < 0))))) CYCLE
                     IF (NORM2(position(:, atom_b) + MATMUL(bs_env%hmat, REAL(s, dp)) - position(:, atom_a)) > &
                         ao_radius(atom_a) + ao_radius(atom_b)) CYCLE
                     DO i = 1, SIZE(images)
                        IF (images(i)%atom /= atom_a) CYCLE
                        index_3d = images(i)%cell + s - shift
                        IF (ANY(index_3d < 1) .OR. ANY(index_3d > upper)) CYCLE
                        IF (image_of(atom_b, index_3d(1), index_3d(2), index_3d(3)) == 0) CYCLE
                        n_keys = n_keys + 1
                        IF (n_keys > SIZE(key_list, 2)) THEN
                           ALLOCATE (grown(5, 2*SIZE(key_list, 2)))
                           grown(:, 1:n_keys - 1) = key_list(:, 1:n_keys - 1)
                           CALL MOVE_ALLOC(grown, key_list)
                        END IF
                        key_list(:, n_keys) = [atom_a, atom_b, s]
                        EXIT
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      ALLOCATE (grown(5, n_keys), SOURCE=key_list(:, 1:n_keys))
      CALL MOVE_ALLOC(grown, key_list)
      DEALLOCATE (ao_radius, position)

      CALL timestop(handle)

   END SUBROUTINE list_unique_ao_pairs

! **************************************************************************************************
!> \brief Makes sure that a work buffer holds at least n values.
!> \param buffer ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ensure_buffer_size(buffer, n)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: buffer(:)
      INTEGER, INTENT(IN)                                :: n

      IF (ALLOCATED(buffer)) THEN
         IF (SIZE(buffer) >= n) RETURN
         DEALLOCATE (buffer)
      END IF
      ALLOCATE (buffer(n))

   END SUBROUTINE ensure_buffer_size

! **************************************************************************************************
!> \brief Contracted overlap S_(μ0,νS') between the AOs of atom a and atom b in cell S'.
!> \param bs_env ...
!> \param a ...
!> \param b ...
!> \param s ...
!> \param overlap overlap(μ, ν).
! **************************************************************************************************
   SUBROUTINE contracted_overlap(bs_env, a, b, s, overlap)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: a, b, s(3)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: overlap(:, :)

      INTEGER                                            :: first_a, first_b, iset, jset, n_a, n_b, &
                                                            nco_a, nco_b
      REAL(dp)                                           :: rab(3)
      REAL(dp), ALLOCATABLE                              :: sab(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b

      cell => bs_env%ri_rs%cell
      basis_a => bs_env%basis_set_AO(bs_env%ri_rs%particle_set(a)%atomic_kind%kind_number)%gto_basis_set
      basis_b => bs_env%basis_set_AO(bs_env%ri_rs%particle_set(b)%atomic_kind%kind_number)%gto_basis_set
      rab = pbc(bs_env%ri_rs%particle_set(b)%r, cell) + MATMUL(bs_env%hmat, REAL(s, dp)) - &
            pbc(bs_env%ri_rs%particle_set(a)%r, cell)
      ALLOCATE (overlap(basis_a%nsgf, basis_b%nsgf))
      overlap = 0.0_dp
      DO iset = 1, basis_a%nset
         nco_a = basis_a%npgf(iset)*(ncoset(basis_a%lmax(iset)) - ncoset(basis_a%lmin(iset) - 1))
         first_a = basis_a%first_sgf(1, iset)
         n_a = basis_a%nsgf_set(iset)
         DO jset = 1, basis_b%nset
            IF (basis_a%set_radius(iset) + basis_b%set_radius(jset) < NORM2(rab)) CYCLE
            nco_b = basis_b%npgf(jset)*(ncoset(basis_b%lmax(jset)) - ncoset(basis_b%lmin(jset) - 1))
            first_b = basis_b%first_sgf(1, jset)
            n_b = basis_b%nsgf_set(jset)
            ALLOCATE (sab(nco_a, nco_b))
            sab = 0.0_dp
            CALL overlap_ab(basis_a%lmax(iset), basis_a%lmin(iset), basis_a%npgf(iset), &
                            basis_a%pgf_radius(:, iset), basis_a%zet(:, iset), &
                            basis_b%lmax(jset), basis_b%lmin(jset), basis_b%npgf(jset), &
                            basis_b%pgf_radius(:, jset), basis_b%zet(:, jset), rab, sab=sab)
            overlap(first_a:first_a + n_a - 1, first_b:first_b + n_b - 1) = &
               MATMUL(TRANSPOSE(basis_a%scon(1:nco_a, first_a:first_a + n_a - 1)), &
                      MATMUL(sab, basis_b%scon(1:nco_b, first_b:first_b + n_b - 1)))
            DEALLOCATE (sab)
         END DO
      END DO

   END SUBROUTINE contracted_overlap

! **************************************************************************************************
!> \brief Number of AO functions of an atom.
!> \param bs_env ...
!> \param iatom ...
!> \return ...
! **************************************************************************************************
   FUNCTION nsgf_atom(bs_env, iatom) RESULT(n)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom
      INTEGER                                            :: n

      n = bs_env%i_ao_end_from_atom(iatom) - bs_env%i_ao_start_from_atom(iatom) + 1

   END FUNCTION nsgf_atom

! **************************************************************************************************
!> \brief Number of reference RI functions of an atom.
!> \param bs_env ...
!> \param iatom ...
!> \return ...
! **************************************************************************************************
   FUNCTION nsgf_ri(bs_env, iatom) RESULT(n)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom
      INTEGER                                            :: n

      n = bs_env%basis_set_RI(bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number)%gto_basis_set%nsgf

   END FUNCTION nsgf_ri

END MODULE ri_rs_kp_compute_d_lP_analytic

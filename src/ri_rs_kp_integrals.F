!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Numerical integrals d_lP(R) for periodic RI-RS fitting.
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE ri_rs_kp_integrals
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_potential_short
   USE kinds,                           ONLY: dp,&
                                              dp_size,&
                                              int_8,&
                                              int_size
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE mathlib,                         ONLY: diamat_all,&
                                              gauss_legendre
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE omp_lib,                         ONLY: omp_get_max_threads,&
                                              omp_in_parallel
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset,&
                                              nso
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE ri_rs_kp_types,                  ONLY: ao_image_type,&
                                              block_class_cache_type,&
                                              l_block_partition_type,&
                                              radial_interpolation_type,&
                                              uniform_settings_type
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_kp_integrals'

   ! Uniform d_lP integration grid: spacing prefactor calibrated on the AO exponents
   ! (including the matched reference RI basis) and safety factor of the Gaussian padding.
   REAL(dp), PARAMETER :: spacing_prefactor = 0.4_dp*SQRT(2.693604_dp*LOG(1.0E8_dp)), &
                          padding_safety = 1.1_dp
   ! Largest Gaussian phase error alpha*|delta r**2| accepted for a tabulated radial moment.
   REAL(dp), PARAMETER :: eps_mesh_phase = 1.0E-8_dp
   ! Tolerance for integer lattice-metric ratios that allow tabulated radial moments.
   REAL(dp), PARAMETER :: eps_integer = 1.0E-10_dp
   ! Radial intervals of the tabulated bounds of the erfc RI potentials [bohr] and their number.
   REAL(dp), PARAMETER :: short_bound_spacing = 0.125_dp
   INTEGER, PARAMETER :: n_short_bounds = 1024
   ! Largest edge of the quadrature blocks of integration points.
   INTEGER, PARAMETER :: quadrature_block_edge = 8, quadrature_block_size = quadrature_block_edge**3
   ! Relative roundoff inflation of the tabulated potential bounds.
   REAL(dp), PARAMETER :: eps_bound_roundoff = 1.0E-10_dp
   ! Short-range radial grid in x = a*r**2: uniform spacing up to short_grid_switch, then
   ! geometric with ratio short_grid_ratio up to short_grid_xmax.
   REAL(dp), PARAMETER :: short_grid_spacing = 0.02_dp, short_grid_switch = 2.0_dp, &
                          short_grid_ratio = 1.01_dp, short_grid_xmax = 512.0_dp
   INTEGER, PARAMETER :: n_short_uniform = NINT(short_grid_switch/short_grid_spacing)
   ! Below this b = omega**2/(a+omega**2), far-field radial moments use the upward recurrence.
   REAL(dp), PARAMETER :: recurrence_b_max = 0.9_dp
   ! Radius [bohr] covered by the uniform-mesh radial moment tables.
   REAL(dp), PARAMETER :: mesh_table_radius = 40.0_dp
   ! Memory estimates: bytes per real, per MiB and per GB (mp_mem_avail_per_rank_GB).
   REAL(dp), PARAMETER :: real_bytes = REAL(dp_size, dp), bytes_per_MiB = 1024.0_dp**2, &
                          bytes_per_GB = 1.0E9_dp
   ! Fixed workspaces: library calls, OpenMP threads, MPI message batches, radial tables of all atoms.
   REAL(dp), PARAMETER :: library_workspace_bytes = 96.0_dp*bytes_per_MiB, &
                          thread_workspace_bytes = 128.0_dp*bytes_per_MiB, &
                          message_batch_bytes = 32.0_dp*bytes_per_MiB, &
                          radial_table_budget = 256.0_dp*bytes_per_MiB
   ! Safety factor on the estimated d_lP workspaces and the shares of the available memory
   ! for the cached RI potentials and AO values of the quadrature blocks.
   REAL(dp), PARAMETER :: workspace_safety = 1.3_dp, ri_cache_share = 0.4_dp, ao_cache_share = 0.1_dp

   PUBLIC :: compute_d_lP_R, &
             compute_d_lP_cutoffs, &
             get_local_l_range, &
             write_ri_rs_progress

CONTAINS

! **************************************************************************************************
!> \brief Allocate d_lP(R), prepare its integration grid, and evaluate the integrals.
!> \param bs_env ...
!> \param grid_points ...
!> \param cells_R_d_lP Retained integer cell vectors R, indexed (3,R).
!> \param d_lP_R Numerical integrals indexed by local point l, RI function P and cell R.
!> \param qpoints ...
!> \param charge_lq ...
! **************************************************************************************************
   SUBROUTINE compute_d_lP_R(bs_env, grid_points, cells_R_d_lP, d_lP_R, &
                             qpoints, charge_lq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: d_lP_R(:, :, :)
      REAL(dp), INTENT(IN)                               :: qpoints(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: charge_lq(:, :)

      INTEGER                                            :: first_l, handle, last_l
      TYPE(uniform_settings_type)                        :: settings

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,A9)') 'd_lP method', 'NUMERICAL'
      END IF
      CALL get_local_l_range(bs_env%ri_rs%n_grid_points, bs_env%para_env%num_pe, &
                             bs_env%para_env%mepos, first_l, last_l)
      ALLOCATE (d_lP_R(last_l - first_l + 1, bs_env%n_RI, SIZE(cells_R_d_lP, 2)))

      CALL timeset("d_lP_uniform_grid", handle)

      CALL prepare_d_lP_integration_grid(bs_env, grid_points, settings)
      ALLOCATE (charge_lq(bs_env%ri_rs%n_grid_points, SIZE(qpoints, 2)))
      CALL integrate_d_lP_R(bs_env, grid_points, cells_R_d_lP, settings, d_lP_R, first_l, &
                            qpoints, charge_lq)

      IF (.NOT. ALL(ieee_is_finite(d_lP_R))) CPABORT("Nonfinite DLP cell value")

      CALL timestop(handle)

   END SUBROUTINE compute_d_lP_R

! **************************************************************************************************
!> \brief Automatic l/P cutoff from the overlap of the diffuse AO-field square and RI potential.
!>
!> Outside the metric cutoff sphere, conservative Gaussian envelopes are
!> f_l(r)^2 ~ exp[-alpha_AO |r-r_l|^2] and
!> v_PT(r) ~ exp[-alpha_RI (|r-R_PT|-r_c)^2/2]. Maximizing their product over r gives
!> exp[-gamma (|r_l-R_PT|-r_c)^2], with 1/gamma=1/alpha_AO+2/alpha_RI.
!> Both methods use the same EPS_FILTER-dependent lattice domain.
!> \param bs_env ...
!> \param cutoffs ...
! **************************************************************************************************
   SUBROUTINE compute_d_lP_cutoffs(bs_env, cutoffs)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: cutoffs(:)

      INTEGER                                            :: iatom, ikind
      REAL(dp)                                           :: accuracy, ao_max, ao_min, ri_max, ri_min
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)

      particles => bs_env%ri_rs%particle_set
      accuracy = d_lP_integration_accuracy(bs_env)
      ALLOCATE (cutoffs(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         ikind = particles(iatom)%atomic_kind%kind_number
         basis => bs_env%basis_set_AO(ikind)%gto_basis_set
         CALL basis_exponent_range(basis, ao_min, ao_max)
         basis => bs_env%basis_set_RI(ikind)%gto_basis_set
         CALL basis_exponent_range(basis, ri_min, ri_max)
         cutoffs(iatom) = MAX(0.0_dp, bs_env%ri_metric%cutoff_radius) + &
                          SQRT(LOG(1.0_dp/accuracy)*(1.0_dp/ao_min + 2.0_dp/ri_min))
      END DO

   END SUBROUTINE compute_d_lP_cutoffs

! **************************************************************************************************
!> \brief Distribute individual l rows evenly across MPI ranks; empty ranges are valid.
!> \param n_grid Number of RI-RS grid points.
!> \param n_ranks Number of MPI ranks.
!> \param i_rank MPI rank.
!> \param first_l First locally owned grid row.
!> \param last_l Last locally owned grid row.
! **************************************************************************************************
   SUBROUTINE get_local_l_range(n_grid, n_ranks, i_rank, first_l, last_l)
      INTEGER, INTENT(IN)                                :: n_grid, n_ranks, i_rank
      INTEGER, INTENT(OUT)                               :: first_l, last_l

      CPASSERT(i_rank >= 0 .AND. i_rank < n_ranks)
      first_l = INT(INT(n_grid, int_8)*INT(i_rank, int_8)/INT(n_ranks, int_8)) + 1
      last_l = INT(INT(n_grid, int_8)*INT(i_rank + 1, int_8)/INT(n_ranks, int_8))

   END SUBROUTINE get_local_l_range

! **************************************************************************************************
!> \brief Print one rank-local RI-RS work milestone without synchronizing MPI ranks.
!> Percentages count completed local units; time is the elapsed wall time since the phase began.
!> \param unit_nr ...
!> \param phase ...
!> \param n_done ...
!> \param n_total ...
!> \param next_progress ...
!> \param phase_start ...
! **************************************************************************************************
   SUBROUTINE write_ri_rs_progress(unit_nr, phase, n_done, n_total, next_progress, phase_start)
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=*), INTENT(IN)                       :: phase
      INTEGER(KIND=int_8), INTENT(IN)                    :: n_done, n_total
      INTEGER, INTENT(INOUT)                             :: next_progress
      REAL(KIND=dp), INTENT(IN)                          :: phase_start

      REAL(KIND=dp), PARAMETER :: milestones(9) = [0.1_dp, 0.3_dp, 1.0_dp, 3.0_dp, 10.0_dp, 30.0_dp&
         , 60.0_dp, 90.0_dp, 100.0_dp]

      CHARACTER(LEN=100)                                 :: progress_line
      CHARACTER(LEN=20)                                  :: count_text
      CHARACTER(LEN=32)                                  :: phase_text
      CHARACTER(LEN=64)                                  :: progress_format
      CHARACTER(LEN=8)                                   :: time_text
      INTEGER                                            :: count_width, previous_progress
      REAL(KIND=dp)                                      :: percent_done, phase_time

      IF (unit_nr <= 0 .OR. n_total <= 0) RETURN
      percent_done = 100.0_dp*REAL(n_done, dp)/REAL(n_total, dp)
      previous_progress = next_progress
      DO WHILE (next_progress <= SIZE(milestones))
         IF (percent_done < milestones(next_progress)) EXIT
         next_progress = next_progress + 1
      END DO
      IF (next_progress == previous_progress) RETURN
      CPASSERT(LEN_TRIM(phase) <= 32)
      WRITE (count_text, '(I0)') n_total
      count_width = LEN_TRIM(count_text)
      WRITE (progress_format, '(A,I0,A,I0,A)') &
         '(A32,1X,I', count_width, ',A,I', count_width, ',1X,F5.1,A)'
      phase_text = phase
      WRITE (progress_line, progress_format) &
         phase_text, n_done, ' / ', n_total, percent_done, ' %'
      phase_time = m_walltime() - phase_start
      WRITE (time_text, '(F8.1)') phase_time
      IF (INDEX(time_text, '*') > 0) WRITE (time_text, '(ES8.1E3)') phase_time
      IF (LEN_TRIM(progress_line) > 56) THEN
         WRITE (unit_nr, '(T2,A)') TRIM(phase)
         WRITE (unit_nr, '(T2,A)') TRIM(ADJUSTL(progress_line(34:)))
         progress_line = ''
      END IF
      WRITE (unit_nr, '(T2,A,T60,A,T71,A,A)') TRIM(progress_line), 'Total time:', time_text, ' s'
      FLUSH (unit_nr)

   END SUBROUTINE write_ri_rs_progress

! **************************************************************************************************
!> \brief Choose the smallest fitting divisor of MPI ranks with 30 percent reserve.
!> A memory group has ngroups ranks and jointly holds all T. Corresponding T slices
!> across memory groups divide g. Equal integration ranks
!> share their uniquely evaluated fields across T groups. The caller's P result is already live.
!> Use the existing node-memory estimate to bound temporary workspaces.
!> Budget the maximum simultaneous integration, reduction and redistribution storage.
!> The final redistribution needs an extra receive buffer as large as the live result.
!> Selection masks remain live throughout; block masks have at most nl rows.
!> \param para_env ...
!> \param nl ...
!> \param bs_env ...
!> \param nt ...
!> \param nresult ...
!> \param integration_comm ...
!> \param field_comm ...
!> \param ngroups ...
! **************************************************************************************************
   SUBROUTINE prepare_d_lP_mpi_groups(para_env, nl, bs_env, nt, nresult, integration_comm, field_comm, ngroups)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nl
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: nt
      INTEGER(int_8), INTENT(IN)                         :: nresult
      TYPE(mp_comm_type), INTENT(OUT)                    :: integration_comm, field_comm
      INTEGER, INTENT(OUT)                               :: ngroups

      INTEGER                                            :: igroup, integration_rank, maxnp, &
                                                            nlocal_T, np
      REAL(dp) :: available_GB, cell_bytes, contraction_bytes, cyclic_bytes, local_bytes, &
         reduction_bytes, required_bytes, result_bytes, selection_bytes, workspace_bytes

      np = bs_env%n_RI
      CALL mp_mem_avail_per_rank_GB(para_env, available_GB)
      result_bytes = real_bytes*REAL(nresult, dp)
      CALL para_env%max(result_bytes)
      cell_bytes = real_bytes*REAL(nl, dp)*REAL(np, dp)
      cyclic_bytes = cell_bytes*REAL((nt + para_env%num_pe - 1)/para_env%num_pe, dp)
      maxnp = MAXVAL(bs_env%i_RI_end_from_atom - bs_env%i_RI_start_from_atom + 1)
      contraction_bytes = real_bytes*(REAL(nl, dp)*REAL(maxnp + quadrature_block_size, dp) + &
                                      REAL(quadrature_block_size, dp)*REAL(maxnp, dp))
      workspace_bytes = library_workspace_bytes + &
                        5.0_dp*real_bytes*REAL(nl + 3, dp)*REAL(quadrature_block_size, dp) + &
                        MAX(thread_workspace_bytes, contraction_bytes)
      DO ngroups = 1, para_env%num_pe
         IF (MODULO(para_env%num_pe, ngroups) /= 0) CYCLE
         nlocal_T = (nt + ngroups - 1)/ngroups
         ! The caller allocated d_lP_R before this call. Its bytes are already
         ! reflected in the measured free memory and must not be charged twice.
         local_bytes = cell_bytes*REAL(nlocal_T, dp)
         selection_bytes = 2.0_dp*REAL(STORAGE_SIZE(.TRUE.)/8, dp)* &
                           REAL(nl, dp)*REAL(bs_env%n_atom, dp)*REAL(nlocal_T, dp)
         selection_bytes = selection_bytes + 3.0_dp*REAL(int_size, dp)*REAL(nlocal_T, dp)
         reduction_bytes = MAX(message_batch_bytes, cell_bytes)
         required_bytes = workspace_safety*(selection_bytes + &
                                            MAX(local_bytes + workspace_bytes, &
                                                local_bytes + cyclic_bytes + reduction_bytes, &
                                                2.0_dp*cyclic_bytes + result_bytes + reduction_bytes))
         IF (required_bytes <= available_GB*bytes_per_GB) EXIT
      END DO
      IF (ngroups > para_env%num_pe) THEN
         ngroups = para_env%num_pe
         CPWARN('d_lP workspaces exceed memory with all MPI ranks')
      END IF
      igroup = MODULO(para_env%mepos, ngroups)
      integration_rank = para_env%mepos/ngroups
      CALL integration_comm%from_split(para_env, igroup, key=integration_rank)
      CALL field_comm%from_split(para_env, integration_rank, key=igroup)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T66,F12.1,A)') 'd_lP PSS-aware rank headroom:', available_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T2,A,T61,I20)') 'd_lP memory group size:', ngroups
         WRITE (bs_env%unit_nr, '(T2,A,T61,I20)') 'd_lP independent integration groups:', integration_comm%num_pe
         WRITE (bs_env%unit_nr, '(T2,A,T66,F12.1,A)') 'd_lP existing output storage:', result_bytes*1.0E-9_dp, ' GB'
         WRITE (bs_env%unit_nr, '(T2,A,T66,F12.1,A)') 'd_lP additional memory budget:', required_bytes*1.0E-9_dp, ' GB'
         FLUSH (bs_env%unit_nr)
      END IF

   END SUBROUTINE prepare_d_lP_mpi_groups

! **************************************************************************************************
!> \brief Sum Eq. (2) directly onto ranks owning cyclic T cells, T=r+1+j*N.
!> The first two indices are (L,P); contiguous T blocks are pointer views, with no packing copy.
!> Reduce-scatter retains only each rank's output cells instead of replicating the full sum.
!> \param comm ...
!> \param partial ...
!> \param cells ...
! **************************************************************************************************
   SUBROUTINE sum_d_lP_R_over_mpi_ranks(comm, partial, cells)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      REAL(dp), CONTIGUOUS, INTENT(IN), TARGET           :: partial(:, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT), TARGET         :: cells(:, :, :)

      INTEGER                                            :: block_cells, first_rank, first_T, &
                                                            icycle, last_rank, last_T, nl, nlocal, &
                                                            np, nranks, nt, rank, round_first, &
                                                            round_last
      INTEGER(int_8)                                     :: iround, plane_count
      REAL(dp), CONTIGUOUS, POINTER                      :: receive_block(:), send_block(:, :)
      REAL(dp), TARGET                                   :: empty(0)
      INTEGER                                            :: recv_counts(comm%num_pe)

      nl = SIZE(partial, 1)
      np = SIZE(partial, 2)
      nt = SIZE(partial, 3)
      nranks = comm%num_pe
      rank = comm%mepos
      nlocal = 0
      IF (nt > rank) nlocal = 1 + (nt - 1 - rank)/nranks
      ALLOCATE (cells(nl, np, nlocal))
      IF (nl == 0 .OR. np == 0 .OR. nt == 0) RETURN
      IF (nranks == 1) THEN
         cells(:, :, :) = partial(:, :, :)
         RETURN
      END IF
      plane_count = INT(nl, int_8)*INT(np, int_8)
      IF (plane_count > HUGE(0)) THEN
         CPABORT('DLP integration block exceeds the MPI count range')
      END IF
      block_cells = INT(MAX(1_int_8, INT(message_batch_bytes/real_bytes, int_8)/plane_count))
      block_cells = MIN(block_cells, INT(INT(HUGE(0), int_8)/plane_count))
      icycle = 0
      DO iround = 1_int_8, INT(nt, int_8), INT(nranks, int_8)
         icycle = icycle + 1
         round_first = INT(iround)
         round_last = INT(MIN(INT(nt, int_8), iround + INT(nranks, int_8) - 1_int_8))
         first_T = round_first
         DO WHILE (first_T <= round_last)
            last_T = INT(MIN(INT(round_last, int_8), INT(first_T, int_8) + INT(block_cells, int_8) - 1_int_8))
            first_rank = first_T - round_first
            last_rank = last_T - round_first
            recv_counts = 0
            recv_counts(first_rank + 1:last_rank + 1) = INT(plane_count)
            send_block(1:INT(plane_count), 1:last_T - first_T + 1) => partial(:, :, first_T:last_T)
            receive_block => empty
            IF (rank >= first_rank .AND. rank <= last_rank) THEN
               receive_block(1:INT(plane_count)) => cells(:, :, icycle)
            END IF
            CALL comm%sum_scatter(send_block, receive_block, recv_counts)
            IF (last_T == round_last) EXIT
            first_T = last_T + 1
         END DO
      END DO

   END SUBROUTINE sum_d_lP_R_over_mpi_ranks

! **************************************************************************************************
!> \brief Redistribute Eq. (2) from cyclic T cells to disjoint contiguous l rows.
!> Rank r initially owns T=r+1+j*N. The caller owns first_l:first_l+size(d_lP_R,1)-1.
!> Only the local T slice and local l slice are materialized, never a global DLP replica.
!> \param para_env ...
!> \param cell_values ...
!> \param first_l ...
!> \param d_lP_R Numerical integrals indexed by local point l, RI function P and cell R.
! **************************************************************************************************
   SUBROUTINE distribute_d_lP_R_by_l(para_env, cell_values, first_l, d_lP_R)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), INTENT(IN)                               :: cell_values(:, :, :)
      INTEGER, INTENT(IN)                                :: first_l
      REAL(dp), INTENT(OUT)                              :: d_lP_R(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'distribute_d_lP_R_by_l'

      INTEGER                                            :: handle, il, ip, irank, it, local_it, &
                                                            nlocal_l, nlocal_T, np, nranks, &
                                                            position
      INTEGER(int_8)                                     :: count_wide
      INTEGER, ALLOCATABLE :: cell_counts(:), first_rows(:), recv_counts(:), recv_offsets(:), &
         row_counts(:), send_counts(:), send_offsets(:)
      REAL(dp), ALLOCATABLE                              :: receive_values(:), send_values(:)

      CALL timeset(routineN, handle)

      nranks = para_env%num_pe
      np = SIZE(d_lP_R, 2)
      nlocal_l = SIZE(d_lP_R, 1)
      nlocal_T = SIZE(cell_values, 3)
      ALLOCATE (first_rows(nranks), row_counts(nranks), cell_counts(nranks), &
                send_counts(nranks), send_offsets(nranks), recv_counts(nranks), recv_offsets(nranks))
      CALL para_env%allgather(first_l, first_rows)
      CALL para_env%allgather(nlocal_l, row_counts)
      CALL para_env%allgather(nlocal_T, cell_counts)
      DO irank = 1, nranks
         count_wide = INT(np, int_8)*INT(row_counts(irank), int_8)*INT(nlocal_T, int_8)
         CPASSERT(count_wide <= HUGE(0))
         count_wide = INT(np, int_8)*INT(nlocal_l, int_8)*INT(cell_counts(irank), int_8)
         CPASSERT(count_wide <= HUGE(0))
         send_counts(irank) = np*row_counts(irank)*nlocal_T
         recv_counts(irank) = np*nlocal_l*cell_counts(irank)
      END DO
      CPASSERT(SUM(INT(send_counts, int_8)) <= HUGE(0))
      CPASSERT(SUM(INT(recv_counts, int_8)) <= HUGE(0))
      position = 0
      DO irank = 1, nranks
         send_offsets(irank) = position
         position = position + send_counts(irank)
      END DO
      ALLOCATE (send_values(position))
      position = 0
      DO irank = 1, nranks
         recv_offsets(irank) = position
         position = position + recv_counts(irank)
      END DO
      ALLOCATE (receive_values(position))
      DO irank = 1, nranks
         position = send_offsets(irank)
         DO local_it = 1, nlocal_T
            DO ip = 1, np
               DO il = first_rows(irank), first_rows(irank) + row_counts(irank) - 1
                  position = position + 1
                  send_values(position) = cell_values(il, ip, local_it)
               END DO
            END DO
         END DO
      END DO
      CALL para_env%alltoall(send_values, send_counts, send_offsets, &
                             receive_values, recv_counts, recv_offsets)
      DO irank = 1, nranks
         position = recv_offsets(irank)
         DO local_it = 1, cell_counts(irank)
            it = irank + (local_it - 1)*nranks
            DO ip = 1, np
               DO il = 1, nlocal_l
                  position = position + 1
                  d_lP_R(il, ip, it) = receive_values(position)
               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE distribute_d_lP_R_by_l

! **************************************************************************************************
!> \brief Accumulate Q_lq += w SUM_g f_l(r_g)^2 exp(i q.r_g) with real BLAS products.
!> The real field is the existing f_l(r_g)^2 on one quadrature integration-point block. Cosine and sine
!> columns share a DGEMM; OpenMP partitions disjoint l rows without a reduction.
!> The grid, phase convention, weights, and retained quadrature points are unchanged.
!> \param points ...
!> \param field ...
!> \param q_cart ...
!> \param weight ...
!> \param charge_lq ...
! **************************************************************************************************
   SUBROUTINE accumulate_fourier_charge(points, field, q_cart, weight, charge_lq)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: field(:, :)
      REAL(dp), INTENT(IN)                               :: q_cart(:, :), weight
      COMPLEX(dp), INTENT(INOUT)                         :: charge_lq(:, :)

      INTEGER, PARAMETER                                 :: row_block = 64

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: first_l, ig, il, iq, nblock, ng, nl, nq
      REAL(dp), ALLOCATABLE                              :: charge_parts(:, :), phase_parts(:, :)

      ng = SIZE(points, 2)
      nl = SIZE(field, 2)
      nq = SIZE(q_cart, 2)
      IF (ng == 0 .OR. nl == 0 .OR. nq == 0) RETURN
      ALLOCATE (phase_parts(ng, 2*nq))
      DO iq = 1, nq
         DO ig = 1, ng
            phase = EXP(CMPLX(0.0_dp, SUM(q_cart(:, iq)*points(:, ig)), dp))
            phase_parts(ig, iq) = REAL(phase, dp)
            phase_parts(ig, nq + iq) = AIMAG(phase)
         END DO
      END DO
      !$OMP PARALLEL DEFAULT(NONE) IF(nl > row_block) &
      !$OMP SHARED(field, phase_parts, charge_lq, ng, nl, nq, weight) &
      !$OMP PRIVATE(charge_parts, first_l, il, iq, nblock)
      ALLOCATE (charge_parts(row_block, 2*nq))
      !$OMP DO SCHEDULE(STATIC)
      DO first_l = 1, nl, row_block
         nblock = MIN(row_block, nl - first_l + 1)
         CALL dgemm('T', 'N', nblock, 2*nq, ng, weight, &
                    field(:, first_l:first_l + nblock - 1), SIZE(field, 1), &
                    phase_parts, ng, 0.0_dp, charge_parts, row_block)
         DO iq = 1, nq
            DO il = 1, nblock
               charge_lq(first_l + il - 1, iq) = charge_lq(first_l + il - 1, iq) + &
                                                 CMPLX(charge_parts(il, iq), charge_parts(il, nq + iq), dp)
            END DO
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (charge_parts)
      !$OMP END PARALLEL
      DEALLOCATE (phase_parts)

   END SUBROUTINE accumulate_fourier_charge

! **************************************************************************************************
!> \brief Validate the operator and define the integration box and mesh in Eq. (3).
!> \param bs_env ...
!> \param grid_points ...
!> \param settings ...
! **************************************************************************************************
   SUBROUTINE prepare_d_lP_integration_grid(bs_env, grid_points, settings)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      TYPE(uniform_settings_type), INTENT(OUT)           :: settings

      INTEGER(int_8)                                     :: npoints
      REAL(dp)                                           :: exponent_max, exponent_min

      settings%accuracy = d_lP_integration_accuracy(bs_env)
      ! The full-k RI-RS path always uses the erfc(omega r)/r RI metric with omega > 0.
      IF (bs_env%ri_metric%potential_type /= do_potential_short .OR. &
          .NOT. (bs_env%ri_metric%omega > 0.0_dp)) THEN
         CPABORT('Uniform DLP requires the short-range erfc RI metric with positive omega')
      END IF
      IF (ALL(bs_env%ri_rs%cell%perd == 0)) CPABORT('Uniform DLP requires a periodic cell')
      CALL ao_exponent_range(bs_env, exponent_min, exponent_max)
      settings%l_block_edge = 1.0_dp/angstrom
      settings%spacing = spacing_prefactor/ &
                         SQRT(exponent_max*LOG(1.0_dp/settings%accuracy))
      settings%padding = padding_safety*SQRT(LOG(1.0_dp/settings%accuracy)/exponent_min)
      CALL build_d_lP_integration_grid(bs_env, grid_points, settings)
      npoints = PRODUCT(INT(settings%counts, int_8))
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') 'd_lP uniform quadrature points', REAL(npoints, dp)
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1)') 'd_lP integration threshold', settings%accuracy
         WRITE (bs_env%unit_nr, '(T2,A,T69,F12.3)') 'd_lP largest AO exponent [bohr^-2]', exponent_max
         WRITE (bs_env%unit_nr, '(T2,A,T69,F12.3)') 'd_lP uniform spacing [bohr]', settings%spacing
      END IF

   END SUBROUTINE prepare_d_lP_integration_grid

! **************************************************************************************************
!> \brief Build Eq. (3) on a cell-commensurate lattice without coarsening its largest spacing.
!> The integration parallelepiped contains the original padded Cartesian box.
!> \param bs_env Band-structure environment with RI-RS cell and atoms.
!> \param grid_points ...
!> \param settings ...
! **************************************************************************************************
   SUBROUTINE build_d_lP_integration_grid(bs_env, grid_points, settings)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      TYPE(uniform_settings_type), INTENT(INOUT)         :: settings

      INTEGER                                            :: idir, subdivisions(3)
      REAL(dp)                                           :: box_center(3), box_extent(3), &
                                                            eigenvalues(3), inverse_step(3, 3), &
                                                            metric(3, 3), scale, upper(3)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particles(:)

      cell => bs_env%ri_rs%cell
      particles => bs_env%ri_rs%particle_set
      settings%origin = pbc(particles(1)%r, cell)
      DO idir = 1, 3
         subdivisions(idir) = MAX(1, CEILING(NORM2(cell%hmat(:, idir))/settings%spacing))
         settings%step(:, idir) = cell%hmat(:, idir)/REAL(subdivisions(idir), dp)
      END DO
      metric = MATMUL(TRANSPOSE(settings%step), settings%step)
      CALL diamat_all(metric, eigenvalues, dac=.TRUE.)
      scale = MAX(1.0_dp, SQRT(MAXVAL(eigenvalues))/settings%spacing)
      subdivisions = CEILING(scale*REAL(subdivisions, dp))
      settings%subdivisions = subdivisions
      DO idir = 1, 3
         settings%step(:, idir) = cell%hmat(:, idir)/REAL(subdivisions(idir), dp)
         inverse_step(idir, :) = REAL(subdivisions(idir), dp)*cell%h_inv(idir, :)
      END DO
      settings%weight = ABS(cell%deth)/PRODUCT(REAL(subdivisions, dp))
      box_center = 0.5_dp*(MINVAL(grid_points, DIM=2) + MAXVAL(grid_points, DIM=2)) - settings%origin
      box_extent = 0.5_dp*(MAXVAL(grid_points, DIM=2) - MINVAL(grid_points, DIM=2)) + settings%padding
      settings%lower = REAL(CEILING(MATMUL(inverse_step, box_center) - &
                                    MATMUL(ABS(inverse_step), box_extent)), dp)
      upper = REAL(FLOOR(MATMUL(inverse_step, box_center) + &
                         MATMUL(ABS(inverse_step), box_extent)), dp)
      settings%counts = INT(upper - settings%lower) + 1
      CPASSERT(ALL(settings%counts > 0))

   END SUBROUTINE build_d_lP_integration_grid

! **************************************************************************************************
!> \brief Partition the RI-RS points r_l into Cartesian cubes used by Eq. (4).
!> \param grid_points ...
!> \param edge ...
!> \param blocks ...
! **************************************************************************************************
   SUBROUTINE group_points_l_by_position(grid_points, edge, blocks)
      REAL(dp), INTENT(IN)                               :: grid_points(:, :), edge
      TYPE(l_block_partition_type), INTENT(OUT)          :: blocks

      INTEGER                                            :: cube(3), iblock, idir, il, insert
      INTEGER, ALLOCATABLE                               :: block_of_l(:), counts(:), cubes(:, :), &
                                                            next(:)
      REAL(dp)                                           :: origin(3)

      origin = MINVAL(grid_points, DIM=2)
      ALLOCATE (block_of_l(SIZE(grid_points, 2)), cubes(3, SIZE(grid_points, 2)))
      blocks%nblocks = 0
      DO il = 1, SIZE(grid_points, 2)
         cube = FLOOR((grid_points(:, il) - origin)/edge)
         iblock = 0
         DO insert = 1, blocks%nblocks
            IF (ALL(cubes(:, insert) == cube)) THEN
               iblock = insert
               EXIT
            END IF
         END DO
         IF (iblock == 0) THEN
            blocks%nblocks = blocks%nblocks + 1
            iblock = blocks%nblocks
            cubes(:, iblock) = cube
         END IF
         block_of_l(il) = iblock
      END DO
      ALLOCATE (counts(blocks%nblocks), next(blocks%nblocks), &
                blocks%first(blocks%nblocks + 1), blocks%indices(SIZE(grid_points, 2)), &
                blocks%lower(3, blocks%nblocks), blocks%upper(3, blocks%nblocks))
      counts = 0
      DO il = 1, SIZE(grid_points, 2)
         counts(block_of_l(il)) = counts(block_of_l(il)) + 1
      END DO
      blocks%first(1) = 1
      DO iblock = 1, blocks%nblocks
         blocks%first(iblock + 1) = blocks%first(iblock) + counts(iblock)
      END DO
      next(:) = blocks%first(1:blocks%nblocks)
      blocks%lower = HUGE(1.0_dp)
      blocks%upper = -HUGE(1.0_dp)
      DO il = 1, SIZE(grid_points, 2)
         iblock = block_of_l(il)
         blocks%indices(next(iblock)) = il
         next(iblock) = next(iblock) + 1
         DO idir = 1, 3
            blocks%lower(idir, iblock) = MIN(blocks%lower(idir, iblock), grid_points(idir, il))
            blocks%upper(idir, iblock) = MAX(blocks%upper(idir, iblock), grid_points(idir, il))
         END DO
      END DO

   END SUBROUTINE group_points_l_by_position

! **************************************************************************************************
!> \brief Return the smallest and largest positive AO primitive exponents for the grid heuristic.
!> \param bs_env ...
!> \param exponent_min ...
!> \param exponent_max ...
! **************************************************************************************************
   SUBROUTINE ao_exponent_range(bs_env, exponent_min, exponent_max)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(OUT)                              :: exponent_min, exponent_max

      INTEGER                                            :: ikind, ipgf, iset
      REAL(dp)                                           :: exponent
      TYPE(gto_basis_set_type), POINTER                  :: basis

      exponent_min = HUGE(1.0_dp)
      exponent_max = 0.0_dp
      DO ikind = 1, SIZE(bs_env%basis_set_AO)
         basis => bs_env%basis_set_AO(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis)) CYCLE
         DO iset = 1, basis%nset
            DO ipgf = 1, basis%npgf(iset)
               exponent = basis%zet(ipgf, iset)
               IF (exponent <= 0.0_dp) CYCLE
               exponent_min = MIN(exponent_min, exponent)
               exponent_max = MAX(exponent_max, exponent)
            END DO
         END DO
      END DO
      IF (exponent_max <= 0.0_dp .OR. exponent_min == HUGE(1.0_dp)) THEN
         CPABORT('Uniform DLP found no positive AO primitive exponent')
      END IF

   END SUBROUTINE ao_exponent_range

! **************************************************************************************************
!> \brief Use the sparse-filter threshold EPS_FILTER as the DLP threshold.
!> \param bs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION d_lP_integration_accuracy(bs_env) RESULT(accuracy)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp)                                           :: accuracy

      IF (.NOT. ieee_is_finite(bs_env%eps_filter) .OR. &
          bs_env%eps_filter <= 0.0_dp .OR. bs_env%eps_filter >= 1.0_dp) THEN
         CPABORT('EPS_FILTER must be between zero and one for periodic RI-RS DLP')
      END IF
      accuracy = bs_env%eps_filter
   END FUNCTION d_lP_integration_accuracy

! **************************************************************************************************
!> \brief Return the smallest and largest positive primitive exponent of one Gaussian basis.
!> \param basis ...
!> \param exponent_min ...
!> \param exponent_max ...
! **************************************************************************************************
   SUBROUTINE basis_exponent_range(basis, exponent_min, exponent_max)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(OUT)                              :: exponent_min, exponent_max

      INTEGER                                            :: ipgf, iset
      REAL(dp)                                           :: exponent

      CPASSERT(ASSOCIATED(basis))
      exponent_min = HUGE(1.0_dp)
      exponent_max = 0.0_dp
      DO iset = 1, basis%nset
         DO ipgf = 1, basis%npgf(iset)
            exponent = basis%zet(ipgf, iset)
            IF (exponent <= 0.0_dp) CYCLE
            exponent_min = MIN(exponent_min, exponent)
            exponent_max = MAX(exponent_max, exponent)
         END DO
      END DO
      IF (exponent_max <= 0.0_dp .OR. exponent_min == HUGE(1.0_dp)) THEN
         CPABORT('Automatic DLP cutoff found no positive primitive exponent')
      END IF

   END SUBROUTINE basis_exponent_range

! **************************************************************************************************
!> \brief Construct d_lP(T) by shared uniform quadrature with analytic auxiliary potentials.
!> (1) f_l(r)=sum_mu,R ϕ_mu,R(r_l) ϕ_mu,R(r), using the explicit AO cell set.
!> (2) d_lP(T)=integral f_l(r)^2 v_P,T(r) dr is approximated by |det S| sum_g f_l(r_g)^2 v_P,T(r_g).
!> (3) r_g=r_0+S n covers the bounding box of r_l enlarged by padding. The spacing follows
!> h=C/[sqrt(alpha_max log(1/epsilon))], calibrated to h=0.4 bohr for the validated Si basis
!> at epsilon=1e-8. The padding is 1.1 sqrt[log(1/epsilon)/alpha_min].
!> Periodic cells use r_0 at atom 1 and S=H diag(1/N_i), with integer N_i and ||S||_2<=h.
!> Translation and cell wrapping then only shift the quadrature lattice or relabel its points.
!> (4) Spatial blocks of r_l are screened separately for every P-atom/T pair immediately before
!> contraction. A triangle-inequality product bound controls discarded discrete terms, not
!> quadrature error.
!> The usual l/P/T cutoff removes individual l rows before contraction.
!> Eq. (1) retains all AO pairs in the explicit cell domain; cell-domain and screening errors
!> must be converged together with the numerical quadrature.
!> Integration-point blocks are MPI distributed and each AO field is evaluated once. Bounded batches
!> share these fields within memory groups whose ranks divide T; independent groups divide g.
!> Group sums assemble disjoint cyclic T cells. A final all-to-all
!> redistributes the completed d(T) from T ownership to the caller's l ownership.
!> Each memory group holds one full partial d(T), subject to the memory budget.
!> No global AO-pair matrix is formed.
!> EPS_FILTER controls the grid, finite box and screening threshold together.
!> \param bs_env ...
!> \param grid_points ...
!> \param cells_R_d_lP Retained integer cell vectors R, indexed (3,R).
!> \param settings ...
!> \param d_lP_R Numerical integrals indexed by local point l, RI function P and cell R.
!> \param first_l ...
!> \param qpoints ...
!> \param charge_lq ...
! **************************************************************************************************
   SUBROUTINE integrate_d_lP_R(bs_env, grid_points, cells_R_d_lP, settings, d_lP_R, first_l, &
                               qpoints, charge_lq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      TYPE(uniform_settings_type), INTENT(IN)            :: settings
      REAL(dp), INTENT(OUT)                              :: d_lP_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_l
      REAL(dp), INTENT(IN)                               :: qpoints(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: charge_lq(:, :)

      INTEGER :: block_class, block_point(3), copy(3), counts(3), first_block(3), handle, iatom, &
         ibatch, iblock, ig, igroup, il, indices(quadrature_block_size), irank, it, last_rank, &
         n_blocks(3), n_class(3), next_progress, ngroups, nkeep, nl, np, npoints, nranks, &
         nreceived, point_index(quadrature_block_size), position, quadrature_block_first(3), &
         quadrature_block_points(3), rank, ranks_per_batch
      INTEGER(int_8) :: contraction_products, iround, item, n_completed_quadrature_blocks, &
         n_items, n_local_quadrature_blocks, nrounds, ntotal, potential_points
      INTEGER(int_8), ALLOCATABLE                        :: chunk_first(:)
      INTEGER, ALLOCATABLE :: chunk_class(:), chunk_part(:), chunk_parts(:), kept_per_rank(:), &
         local_cells(:, :), recv_counts(:), recv_offsets(:)
      LOGICAL, ALLOCATABLE                               :: include_l_blocks_P_R(:, :, :), &
                                                            include_lP_R(:, :, :)
      REAL(dp) :: available_GB, bound, elapsed, field_bound(quadrature_block_size), lower(3), &
         mesh_index(3), points(3, quadrature_block_size), prep_start, RI_center_R(3), threshold, &
         weight
      REAL(dp), ALLOCATABLE :: cutoffs_lP(:), cyclic_d_lP_R(:, :, :), envelopes(:, :), &
         field(:, :), local_d_lP_R(:, :, :), q_cart(:, :), receive_fields(:), send_fields(:)
      TYPE(ao_image_type), ALLOCATABLE                   :: images(:)
      TYPE(block_class_cache_type)                       :: ao_cache, cache
      TYPE(cell_type), POINTER                           :: cell
      TYPE(l_block_partition_type)                       :: blocks
      TYPE(mp_comm_type)                                 :: field_comm, integration_comm
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(radial_interpolation_type), ALLOCATABLE       :: aux_potentials(:)

      nl = bs_env%ri_rs%n_grid_points
      cell => bs_env%ri_rs%cell
      particles => bs_env%ri_rs%particle_set
      ALLOCATE (field(quadrature_block_size, nl), q_cart(3, SIZE(qpoints, 2)))
      q_cart(:, :) = MATMUL(twopi*TRANSPOSE(cell%h_inv), qpoints)
      charge_lq = CMPLX(0.0_dp, 0.0_dp, dp)
      elapsed = m_walltime()
      prep_start = elapsed
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'd_lP quadrature: evaluating AO image values'
         FLUSH (bs_env%unit_nr)
      END IF
      CALL prepare_ao_images(bs_env, grid_points, images)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T67,F12.1,A)') &
            'd_lP AO image values, rank 0 time:', m_walltime() - prep_start, ' s'
         FLUSH (bs_env%unit_nr)
      END IF
      CALL prepare_aux_potentials(bs_env, settings, aux_potentials)
      CALL auxiliary_potential_bound(bs_env, aux_potentials, envelopes)
      CALL group_points_l_by_position(grid_points, settings%l_block_edge, blocks)
      bound = MAXVAL(envelopes(1, :))
      d_lP_R = 0.0_dp
      CALL prepare_d_lP_mpi_groups(bs_env%para_env, nl, bs_env, SIZE(cells_R_d_lP, 2), SIZE(d_lP_R, KIND=int_8), &
                                   integration_comm, field_comm, ngroups)
      igroup = field_comm%mepos
      nranks = field_comm%num_pe
      rank = field_comm%mepos
      np = bs_env%n_RI
      ALLOCATE (local_cells(3, SIZE(cells_R_d_lP(:, igroup + 1::ngroups), 2)))
      local_cells(:, :) = cells_R_d_lP(:, igroup + 1::ngroups)

      CALL timeset('d_lP_prepare_RI_geometry', handle)

      CALL compute_d_lP_cutoffs(bs_env, cutoffs_lP)
      ALLOCATE (include_lP_R(nl, bs_env%n_atom, SIZE(local_cells, 2)))
      ALLOCATE (include_l_blocks_P_R(blocks%nblocks, bs_env%n_atom, SIZE(local_cells, 2)), source=.FALSE.)
      DO it = 1, SIZE(local_cells, 2)
         DO iatom = 1, bs_env%n_atom
            RI_center_R = pbc(particles(iatom)%r, cell) + MATMUL(bs_env%hmat, REAL(local_cells(:, it), dp))
            DO il = 1, nl
               include_lP_R(il, iatom, it) = SUM((grid_points(:, il) - RI_center_R)**2) <= cutoffs_lP(iatom)**2
            END DO
            DO iblock = 1, blocks%nblocks
               IF (distance_to_box_squared(RI_center_R, blocks%lower(:, iblock), &
                                           blocks%upper(:, iblock)) > cutoffs_lP(iatom)**2) CYCLE
               include_l_blocks_P_R(iblock, iatom, it) = &
                  ANY(include_lP_R(blocks%indices(blocks%first(iblock):blocks%first(iblock + 1) - 1), iatom, it))
            END DO
         END DO
      END DO
      DEALLOCATE (cutoffs_lP)
      IF (bs_env%unit_nr > 0) THEN
         FLUSH (bs_env%unit_nr)
      END IF

      CALL timestop(handle)

      ALLOCATE (local_d_lP_R(nl, np, SIZE(local_cells, 2)), source=0.0_dp)
      ALLOCATE (kept_per_rank(nranks), recv_counts(nranks), recv_offsets(nranks))
      ! At most 32 MiB of shared fields per rank, independent of the MPI size.
      ranks_per_batch = MAX(1, INT(message_batch_bytes/ &
                                   (real_bytes*REAL(nl + 4, dp)*REAL(quadrature_block_size, dp))))
      ranks_per_batch = MIN(nranks, ranks_per_batch)
      ALLOCATE (send_fields((nl + 4)*quadrature_block_size), &
                receive_fields((nl + 4)*quadrature_block_size*ranks_per_batch))
      contraction_products = 0_int_8
      potential_points = 0_int_8
      weight = settings%weight
      lower = settings%lower
      counts = settings%counts
      ntotal = PRODUCT(INT(counts, int_8))
      ! Each integration group handles parts of classes of cell-aligned blocks, its ranks share the
      ! fields of each block, and every rank keeps the potentials of the current class.
      CALL plan_quadrature_blocks(settings, cell%perd, integration_comm%num_pe, integration_comm%mepos, &
                                  n_class, first_block, n_blocks, chunk_class, chunk_part, chunk_parts, &
                                  chunk_first)
      n_items = chunk_first(SIZE(chunk_class))
      nrounds = (n_items + INT(nranks, int_8) - 1_int_8)/INT(nranks, int_8)
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, available_GB)
      CALL prepare_block_class_cache(cache, local_cells, bs_env%n_atom, quadrature_block_size, &
                                     MAXVAL(bs_env%i_RI_end_from_atom - bs_env%i_RI_start_from_atom + 1), &
                                     n_class, first_block, n_blocks, ri_cache_share*bytes_per_GB*available_GB, .TRUE.)
      CALL prepare_block_class_cache(ao_cache, bs_env%index_to_cell_3c, bs_env%n_atom, quadrature_block_size, &
                                     MAXVAL(bs_env%i_ao_end_from_atom - bs_env%i_ao_start_from_atom + 1), &
                                     n_class, first_block, n_blocks, ao_cache_share*bytes_per_GB*available_GB, .FALSE.)
      next_progress = 1
      n_local_quadrature_blocks = nrounds
      threshold = 0.5_dp*settings%accuracy/(REAL(ntotal, dp)*weight)
      DO iround = 0_int_8, nrounds - 1_int_8
         item = iround*INT(nranks, int_8) + INT(rank, int_8)
         nkeep = 0
         IF (item < n_items) THEN
            CALL locate_quadrature_block(item, settings, cell%perd, n_class, first_block, n_blocks, &
                                         chunk_class, chunk_part, chunk_parts, chunk_first, block_class, &
                                         copy, quadrature_block_first, quadrature_block_points)
            ! Points are indexed in the complete block, so all copies address the same potentials.
            npoints = 0
            DO ig = 1, quadrature_block_size
               block_point(1) = MODULO(ig - 1, quadrature_block_edge)
               block_point(2) = MODULO((ig - 1)/quadrature_block_edge, quadrature_block_edge)
               block_point(3) = (ig - 1)/quadrature_block_edge**2
               IF (ANY(block_point >= quadrature_block_points)) CYCLE
               IF (ANY(quadrature_block_first + block_point < 0 .OR. &
                       quadrature_block_first + block_point >= counts)) CYCLE
               npoints = npoints + 1
               point_index(npoints) = ig
               mesh_index = lower + REAL(quadrature_block_first + block_point, dp)
               points(:, npoints) = settings%origin + MATMUL(settings%step, mesh_index)
            END DO
            CALL timeset('d_lP_f_l_squared', handle)
            IF (block_class /= ao_cache%block_class) CALL reset_block_class_cache(ao_cache, block_class)
            CALL evaluate_rho_l_squared(images, points(:, 1:npoints), point_index(1:npoints), copy, ao_cache, &
                                        field(1:npoints, :), &
                                        settings%accuracy/(8.0_dp*weight*REAL(ntotal, dp)*MAX(1.0_dp, bound)))
            CALL timestop(handle)
            CALL timeset('d_lP_charge_fourier', handle)
            CALL accumulate_fourier_charge(points(:, 1:npoints), field, q_cart, weight, charge_lq)
            CALL timestop(handle)
            CALL timeset('d_lP_screen_pack_fields', handle)
            nkeep = 0
            field_bound(1:npoints) = 0.0_dp
            DO il = 1, nl
               !$OMP SIMD
               DO ig = 1, npoints
                  field_bound(ig) = MAX(field_bound(ig), ABS(field(ig, il)))
               END DO
            END DO
            DO ig = 1, npoints
               IF (field_bound(ig)*bound > threshold) THEN
                  nkeep = nkeep + 1
                  indices(nkeep) = ig
               END IF
            END DO
            DO ig = 1, nkeep
               send_fields(4*ig - 3:4*ig - 1) = points(:, indices(ig))
               send_fields(4*ig) = REAL(point_index(indices(ig)), dp)
            END DO
            DO il = 1, nl
               DO ig = 1, nkeep
                  send_fields(4*nkeep + (il - 1)*nkeep + ig) = field(indices(ig), il)
               END DO
            END DO
            CALL timestop(handle)
         END IF
         CALL field_comm%allgather(nkeep, kept_per_rank)
         DO ibatch = 1, nranks, ranks_per_batch
            last_rank = MIN(nranks, ibatch + ranks_per_batch - 1)
            recv_counts = 0
            recv_counts(ibatch:last_rank) = (nl + 4)*kept_per_rank(ibatch:last_rank)
            nreceived = 0
            DO irank = 1, nranks
               recv_offsets(irank) = nreceived
               nreceived = nreceived + recv_counts(irank)
            END DO
            IF (nreceived == 0) CYCLE
            CPASSERT(nreceived <= SIZE(receive_fields))
            CALL timeset('d_lP_uniform_exchange', handle)
            CALL field_comm%allgatherv(send_fields(1:recv_counts(rank + 1)), &
                                       receive_fields(1:nreceived), recv_counts, recv_offsets)
            CALL timestop(handle)
            DO irank = ibatch, last_rank
               nkeep = kept_per_rank(irank)
               IF (nkeep == 0 .OR. SIZE(local_cells, 2) == 0) CYCLE
               CALL timeset('d_lP_unpack_fields', handle)
               position = recv_offsets(irank)
               DO ig = 1, nkeep
                  points(:, ig) = receive_fields(position + 4*ig - 3:position + 4*ig - 1)
                  point_index(ig) = NINT(receive_fields(position + 4*ig))
               END DO
               DO il = 1, nl
                  DO ig = 1, nkeep
                     field(ig, il) = receive_fields(position + 4*nkeep + (il - 1)*nkeep + ig)
                  END DO
               END DO
               CALL locate_quadrature_block(iround*INT(nranks, int_8) + INT(irank - 1, int_8), settings, &
                                            cell%perd, n_class, first_block, n_blocks, chunk_class, &
                                            chunk_part, chunk_parts, chunk_first, block_class, copy, &
                                            quadrature_block_first, quadrature_block_points)
               IF (block_class /= cache%block_class) CALL reset_block_class_cache(cache, block_class)
               CALL timestop(handle)
               CALL timeset('d_lP_integrate_quadrature_blocks', handle)
               CALL integrate_potentials(bs_env, local_cells, include_lP_R, include_l_blocks_P_R, &
                                         points(:, 1:nkeep), point_index(1:nkeep), field(1:nkeep, :), &
                                         weight, envelopes, blocks, threshold, copy, cache, potential_points, &
                                         contraction_products, local_d_lP_R, aux_potentials)
               CALL timestop(handle)
            END DO
         END DO
         n_completed_quadrature_blocks = iround + 1_int_8
         CALL write_ri_rs_progress(bs_env%unit_nr, 'd_lP Integration batch done', &
                                   n_completed_quadrature_blocks, n_local_quadrature_blocks, next_progress, elapsed)
      END DO

      CALL timeset('d_lP_uniform_reduce', handle)

      DEALLOCATE (field, send_fields, receive_fields, kept_per_rank, recv_counts, recv_offsets, &
                  chunk_class, chunk_part, chunk_parts, chunk_first)
      CALL sum_d_lP_R_over_mpi_ranks(integration_comm, local_d_lP_R, cyclic_d_lP_R)
      DEALLOCATE (local_d_lP_R, local_cells)
      CALL distribute_d_lP_R_by_l(bs_env%para_env, cyclic_d_lP_R, first_l, d_lP_R)
      DEALLOCATE (cyclic_d_lP_R)
      CALL integration_comm%free()
      CALL field_comm%free()
      CALL bs_env%para_env%sum(charge_lq)
      CALL bs_env%para_env%sum(potential_points)
      CALL bs_env%para_env%sum(contraction_products)

      CALL timestop(handle)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') &
            'Potential evaluations V_P^R(r_g) (point, atom image):', REAL(potential_points, dp)
         WRITE (bs_env%unit_nr, '(T2,A,T72,ES9.1E2)') &
            'd_lP products f_l(r_g)^2 V_P^R(r_g):', REAL(contraction_products, dp)
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'f_l(r) = SUM_(μ,R) ϕ_μ^R(r_l) ϕ_μ^R(r)'
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'd_lP^R = SUM_g w_g f_l(r_g)^2 V_P^R(r_g)'
      END IF
      elapsed = m_walltime() - elapsed
      CALL bs_env%para_env%max(elapsed)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') &
         'd_lP completed numerical integration, Exec. time:', elapsed, ' s'
      IF (bs_env%unit_nr > 0) FLUSH (bs_env%unit_nr)

   END SUBROUTINE integrate_d_lP_R

! **************************************************************************************************
!> \brief Distribute the cell-aligned quadrature blocks over the integration groups by class.
!> In a periodic direction, blocks of up to quadrature_block_edge points tile each cell period of
!> N subdivisions. Block j = c + M t is block c of the M blocks of period t, so all lattice copies t
!> of the class c see the same RI potentials. In a non-periodic direction every block is a class.
!> A group takes whole classes; a class with more blocks than a group share is split into strided
!> parts, each of which spans the whole integration box.
!> \param settings ...
!> \param perd periodicity of the cell
!> \param n_groups ...
!> \param my_group ...
!> \param n_class classes M per direction
!> \param first_block index j of the first box block per direction
!> \param n_blocks box blocks per direction
!> \param chunk_class class of each part of this group (from 0)
!> \param chunk_part ...
!> \param chunk_parts ...
!> \param chunk_first number of blocks of this group before each part, indexed from 0
! **************************************************************************************************
   SUBROUTINE plan_quadrature_blocks(settings, perd, n_groups, my_group, n_class, first_block, &
                                     n_blocks, chunk_class, chunk_part, chunk_parts, chunk_first)
      TYPE(uniform_settings_type), INTENT(IN)            :: settings
      INTEGER, INTENT(IN)                                :: perd(3), n_groups, my_group
      INTEGER, INTENT(OUT)                               :: n_class(3), first_block(3), n_blocks(3)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: chunk_class(:), chunk_part(:), &
                                                            chunk_parts(:)
      INTEGER(int_8), ALLOCATABLE, INTENT(OUT)           :: chunk_first(:)

      INTEGER :: block_class, copy_lower(3), group, ichunk, idir, ipart, last_point, n_chunks, &
         n_copies(3), n_mine, n_parts, period, position(3)
      INTEGER(int_8)                                     :: class_blocks, part_length
      INTEGER(int_8), ALLOCATABLE                        :: load(:), part_blocks(:)
      INTEGER, ALLOCATABLE                               :: owner(:), part_class(:), part_count(:), &
                                                            part_index(:)

      DO idir = 1, 3
         IF (perd(idir) == 1) THEN
            n_class(idir) = (settings%subdivisions(idir) + quadrature_block_edge - 1)/quadrature_block_edge
            period = FLOOR(settings%lower(idir)/REAL(settings%subdivisions(idir), dp))
            first_block(idir) = n_class(idir)*period + &
                                (NINT(settings%lower(idir)) - settings%subdivisions(idir)*period)/quadrature_block_edge
            last_point = NINT(settings%lower(idir)) + settings%counts(idir) - 1
            period = FLOOR(REAL(last_point, dp)/REAL(settings%subdivisions(idir), dp))
            n_blocks(idir) = n_class(idir)*period + &
                             (last_point - settings%subdivisions(idir)*period)/quadrature_block_edge - &
                             first_block(idir) + 1
         ELSE
            n_blocks(idir) = (settings%counts(idir) + quadrature_block_edge - 1)/quadrature_block_edge
            n_class(idir) = n_blocks(idir)
            first_block(idir) = 0
         END IF
      END DO
      ! About 64 parts per group balance the load and still reuse each potential many times.
      part_length = MAX(1_int_8, PRODUCT(INT(n_blocks, int_8))/(64_int_8*INT(n_groups, int_8)))
      n_chunks = 0
      DO block_class = 0, PRODUCT(n_class) - 1
         CALL class_copy_range(block_class, n_class, first_block, n_blocks, position, copy_lower, n_copies)
         class_blocks = PRODUCT(INT(n_copies, int_8))
         n_chunks = n_chunks + INT((class_blocks + part_length - 1_int_8)/part_length)
      END DO
      ALLOCATE (part_class(n_chunks), part_index(n_chunks), part_count(n_chunks), &
                part_blocks(n_chunks), owner(n_chunks), load(0:n_groups - 1))
      load = 0_int_8
      ichunk = 0
      DO block_class = 0, PRODUCT(n_class) - 1
         CALL class_copy_range(block_class, n_class, first_block, n_blocks, position, copy_lower, n_copies)
         class_blocks = PRODUCT(INT(n_copies, int_8))
         n_parts = INT((class_blocks + part_length - 1_int_8)/part_length)
         DO ipart = 0, n_parts - 1
            ichunk = ichunk + 1
            part_class(ichunk) = block_class
            part_index(ichunk) = ipart
            part_count(ichunk) = n_parts
            part_blocks(ichunk) = (class_blocks - INT(ipart, int_8) + INT(n_parts - 1, int_8))/INT(n_parts, int_8)
            group = MINLOC(load, 1) - 1
            owner(ichunk) = group
            load(group) = load(group) + part_blocks(ichunk)
         END DO
      END DO
      n_mine = COUNT(owner == my_group)
      ALLOCATE (chunk_class(n_mine), chunk_part(n_mine), chunk_parts(n_mine), chunk_first(0:n_mine))
      chunk_first(0) = 0_int_8
      n_mine = 0
      DO ichunk = 1, n_chunks
         IF (owner(ichunk) /= my_group) CYCLE
         n_mine = n_mine + 1
         chunk_class(n_mine) = part_class(ichunk)
         chunk_part(n_mine) = part_index(ichunk)
         chunk_parts(n_mine) = part_count(ichunk)
         chunk_first(n_mine) = chunk_first(n_mine - 1) + part_blocks(ichunk)
      END DO
      DEALLOCATE (part_class, part_index, part_count, part_blocks, owner, load)

   END SUBROUTINE plan_quadrature_blocks

! **************************************************************************************************
!> \brief Position of a block class and the range of its lattice copies in the box.
!> \param block_class class index from 0
!> \param n_class ...
!> \param first_block ...
!> \param n_blocks ...
!> \param position c per direction
!> \param copy_lower first copy t per direction
!> \param n_copies copies per direction
! **************************************************************************************************
   PURE SUBROUTINE class_copy_range(block_class, n_class, first_block, n_blocks, position, &
                                    copy_lower, n_copies)
      INTEGER, INTENT(IN)                                :: block_class, n_class(3), first_block(3), &
                                                            n_blocks(3)
      INTEGER, INTENT(OUT)                               :: position(3), copy_lower(3), n_copies(3)

      INTEGER                                            :: idir

      position(1) = MODULO(block_class, n_class(1))
      position(2) = MODULO(block_class/n_class(1), n_class(2))
      position(3) = block_class/(n_class(1)*n_class(2))
      DO idir = 1, 3
         copy_lower(idir) = CEILING(REAL(first_block(idir) - position(idir), dp)/REAL(n_class(idir), dp))
         n_copies(idir) = MAX(0, FLOOR(REAL(first_block(idir) + n_blocks(idir) - 1 - position(idir), dp)/ &
                                       REAL(n_class(idir), dp)) - copy_lower(idir) + 1)
      END DO

   END SUBROUTINE class_copy_range

! **************************************************************************************************
!> \brief Class, lattice copy and box points of the item-th block of this integration group.
!> \param item block of the group, from 0
!> \param settings ...
!> \param perd ...
!> \param n_class ...
!> \param first_block ...
!> \param n_blocks ...
!> \param chunk_class ...
!> \param chunk_part ...
!> \param chunk_parts ...
!> \param chunk_first ...
!> \param block_class ...
!> \param copy lattice vector t of the copy
!> \param block_first first point of the block in the box per direction, may lie before the box
!> \param block_points points of the block per direction
! **************************************************************************************************
   SUBROUTINE locate_quadrature_block(item, settings, perd, n_class, first_block, n_blocks, &
                                      chunk_class, chunk_part, chunk_parts, chunk_first, &
                                      block_class, copy, block_first, block_points)
      INTEGER(int_8), INTENT(IN)                         :: item
      TYPE(uniform_settings_type), INTENT(IN)            :: settings
      INTEGER, INTENT(IN)                                :: perd(3), n_class(3), first_block(3), &
                                                            n_blocks(3), chunk_class(:), &
                                                            chunk_part(:), chunk_parts(:)
      INTEGER(int_8), INTENT(IN)                         :: chunk_first(0:)
      INTEGER, INTENT(OUT)                               :: block_class, copy(3), block_first(3), &
                                                            block_points(3)

      INTEGER                                            :: copy_lower(3), ichunk, idir, &
                                                            lower_chunk, n_copies(3), position(3), &
                                                            upper_chunk
      INTEGER(int_8)                                     :: k

      lower_chunk = 1
      upper_chunk = SIZE(chunk_class)
      DO WHILE (lower_chunk < upper_chunk)
         ichunk = (lower_chunk + upper_chunk)/2
         IF (chunk_first(ichunk) <= item) THEN
            lower_chunk = ichunk + 1
         ELSE
            upper_chunk = ichunk
         END IF
      END DO
      ichunk = lower_chunk
      block_class = chunk_class(ichunk)
      CALL class_copy_range(block_class, n_class, first_block, n_blocks, position, copy_lower, n_copies)
      k = INT(chunk_part(ichunk), int_8) + INT(chunk_parts(ichunk), int_8)*(item - chunk_first(ichunk - 1))
      copy(1) = copy_lower(1) + INT(MODULO(k, INT(n_copies(1), int_8)))
      copy(2) = copy_lower(2) + INT(MODULO(k/INT(n_copies(1), int_8), INT(n_copies(2), int_8)))
      copy(3) = copy_lower(3) + INT(k/INT(n_copies(1)*n_copies(2), int_8))
      DO idir = 1, 3
         IF (perd(idir) == 1) THEN
            block_first(idir) = settings%subdivisions(idir)*copy(idir) + &
                                quadrature_block_edge*position(idir) - NINT(settings%lower(idir))
            block_points(idir) = MIN(quadrature_block_edge, &
                                     settings%subdivisions(idir) - quadrature_block_edge*position(idir))
         ELSE
            block_first(idir) = quadrature_block_edge*position(idir)
            block_points(idir) = quadrature_block_edge
         END IF
      END DO

   END SUBROUTINE locate_quadrature_block

! **************************************************************************************************
!> \brief Allocate the slot map over all U = T - t of this rank's cells T and the box copies t.
!> \param cache ...
!> \param cells ...
!> \param n_atom ...
!> \param n_points points of a complete block
!> \param maxnp functions per atom
!> \param n_class ...
!> \param first_block ...
!> \param n_blocks ...
!> \param budget bytes the values may occupy
!> \param point_bounds whether the slots also keep point bounds
! **************************************************************************************************
   SUBROUTINE prepare_block_class_cache(cache, cells, n_atom, n_points, maxnp, n_class, first_block, &
                                        n_blocks, budget, point_bounds)
      TYPE(block_class_cache_type), INTENT(OUT)          :: cache
      INTEGER, INTENT(IN)                                :: cells(:, :), n_atom, n_points, maxnp, &
                                                            n_class(3), first_block(3), n_blocks(3)
      REAL(dp), INTENT(IN)                               :: budget
      LOGICAL, INTENT(IN)                                :: point_bounds

      INTEGER                                            :: copy_lower(3), copy_upper(3), &
                                                            u_lower(3), u_upper(3)

      copy_lower = FLOOR(REAL(first_block, dp)/REAL(n_class, dp))
      copy_upper = FLOOR(REAL(first_block + n_blocks - 1, dp)/REAL(n_class, dp))
      u_lower = 0
      u_upper = 0
      IF (SIZE(cells, 2) > 0) THEN
         u_lower = MINVAL(cells, DIM=2) - copy_upper
         u_upper = MAXVAL(cells, DIM=2) - copy_lower
      END IF
      ALLOCATE (cache%slot(u_lower(1):u_upper(1), u_lower(2):u_upper(2), u_lower(3):u_upper(3), n_atom), &
                SOURCE=0)
      ! Only the pages of used slots are touched, so the budget is reserved, not occupied.
      cache%n_slots = INT(MIN(REAL(SIZE(cache%slot), dp), &
                              budget/(real_bytes*REAL(n_points, dp)*REAL(maxnp + 2, dp))))
      ALLOCATE (cache%keys(4, cache%n_slots), cache%evaluated(n_points, cache%n_slots), &
                cache%values(maxnp, n_points, cache%n_slots))
      IF (point_bounds) ALLOCATE (cache%point_bound(n_points, cache%n_slots))

   END SUBROUTINE prepare_block_class_cache

! **************************************************************************************************
!> \brief Start the values of a new block class.
!> \param cache ...
!> \param block_class ...
! **************************************************************************************************
   SUBROUTINE reset_block_class_cache(cache, block_class)
      TYPE(block_class_cache_type), INTENT(INOUT)        :: cache
      INTEGER, INTENT(IN)                                :: block_class

      INTEGER                                            :: islot

      DO islot = 1, cache%n_used
         cache%slot(cache%keys(1, islot), cache%keys(2, islot), cache%keys(3, islot), cache%keys(4, islot)) = 0
      END DO
      cache%n_used = 0
      cache%block_class = block_class

   END SUBROUTINE reset_block_class_cache

! **************************************************************************************************
!> \brief Form each atom-local auxiliary potential's spherical coefficients just once.
!> Without AUTO_RI, the contraction U_Pp is the identity.
!> \param bs_env ...
!> \param settings ...
!> \param aux_potentials ...
! **************************************************************************************************
   SUBROUTINE prepare_aux_potentials(bs_env, settings, aux_potentials)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(uniform_settings_type), INTENT(IN)            :: settings
      TYPE(radial_interpolation_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: aux_potentials(:)

      CHARACTER(LEN=12)                                  :: spacing_text
      INTEGER                                            :: first_U, iatom, ikind, ip, nmesh, np, &
                                                            nref
      REAL(dp)                                           :: mesh_bytes, mesh_displacement(3), &
                                                            mesh_spacing
      REAL(dp), ALLOCATABLE                              :: U_Pp(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)

      particles => bs_env%ri_rs%particle_set
      cell => bs_env%ri_rs%cell
      ALLOCATE (aux_potentials(bs_env%n_atom))
      nmesh = 0
      mesh_bytes = 0.0_dp
      mesh_spacing = 0.0_dp
      DO iatom = 1, bs_env%n_atom
         ikind = particles(iatom)%atomic_kind%kind_number
         basis => bs_env%basis_set_RI(ikind)%gto_basis_set
         IF (bs_env%auto_ri%enabled) THEN
            nref = bs_env%auto_ri%sizes_ref_RI(iatom)
            np = bs_env%auto_ri%sizes_opt_RI(iatom)
            first_U = bs_env%auto_ri%U_Pp_atom_offset(iatom)
            ALLOCATE (U_Pp(nref, np))
            U_Pp(:, :) = RESHAPE(bs_env%auto_ri%U_Pp_atom(first_U:first_U + nref*np - 1), [nref, np])
         ELSE
            ALLOCATE (U_Pp(basis%nsgf, basis%nsgf), source=0.0_dp)
            DO ip = 1, basis%nsgf
               U_Pp(ip, ip) = 1.0_dp
            END DO
         END IF
         mesh_displacement = settings%origin - pbc(particles(iatom)%r, cell)
         CALL prepare_contracted_potential(basis, U_Pp, aux_potentials(iatom), &
                                           MAX(0.0_dp, bs_env%ri_metric%omega), &
                                           radial_table_bytes=radial_table_budget/REAL(MAX(1, bs_env%n_atom), dp), &
                                           mesh_step=settings%step, &
                                           mesh_displacement=mesh_displacement)
         IF (ALLOCATED(aux_potentials(iatom)%mesh_moments)) THEN
            nmesh = nmesh + 1
            mesh_bytes = mesh_bytes + real_bytes*REAL(SIZE(aux_potentials(iatom)%mesh_moments, KIND=int_8), dp)
            mesh_spacing = aux_potentials(iatom)%mesh_spacing
         END IF
         DEALLOCATE (U_Pp)
      END DO
      IF (bs_env%unit_nr > 0) THEN
         IF (nmesh > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T69,F12.1)') 'd_lP radial mesh storage [MiB/rank]', mesh_bytes/bytes_per_MiB
            WRITE (spacing_text, '(ES12.1E1)') mesh_spacing
            IF (INDEX(spacing_text, '*') > 0) WRITE (spacing_text, '(ES12.1E2)') mesh_spacing
            IF (INDEX(spacing_text, '*') > 0) WRITE (spacing_text, '(ES12.1E3)') mesh_spacing
            WRITE (bs_env%unit_nr, '(T2,A,T69,A12)') 'd_lP radial mesh squared spacing [bohr^2]', spacing_text
         END IF
      END IF

   END SUBROUTINE prepare_aux_potentials

! **************************************************************************************************
!> \brief Store nonzero atom-image coefficients ϕ_mu,R(r_l) in Eq. (1) once.
!> \param bs_env ...
!> \param grid_points ...
!> \param images ...
! **************************************************************************************************
   SUBROUTINE prepare_ao_images(bs_env, grid_points, images)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(IN)                               :: grid_points(:, :)
      TYPE(ao_image_type), ALLOCATABLE, INTENT(OUT)      :: images(:)

      INTEGER                                            :: iatom, icell, iimage, ikind, maxao, nao
      REAL(dp)                                           :: center(3), cutoff2
      REAL(dp), ALLOCATABLE                              :: primary(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)

      cell => bs_env%ri_rs%cell
      particles => bs_env%ri_rs%particle_set
      maxao = MAXVAL(bs_env%i_ao_end_from_atom - bs_env%i_ao_start_from_atom + 1)
      ALLOCATE (primary(bs_env%ri_rs%n_grid_points, maxao))
      ALLOCATE (images(bs_env%n_atom*SIZE(bs_env%index_to_cell_3c, 2)))
      iimage = 0
      DO icell = 1, SIZE(bs_env%index_to_cell_3c, 2)
         DO iatom = 1, bs_env%n_atom
            ikind = particles(iatom)%atomic_kind%kind_number
            basis => bs_env%basis_set_AO(ikind)%gto_basis_set
            nao = bs_env%i_ao_end_from_atom(iatom) - bs_env%i_ao_start_from_atom(iatom) + 1
            center = pbc(particles(iatom)%r, cell) + MATMUL(bs_env%hmat, REAL(bs_env%index_to_cell_3c(:, icell), dp))
            cutoff2 = bs_env%ri_rs%radius_ao_per_atom(iatom)**2
            IF (bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) cutoff2 = bs_env%ri_rs%cutoff_radius_ri_ao**2
            primary = 0.0_dp
            CALL evaluate_ao_basis_on_points(primary(:, 1:nao), grid_points, basis, center, cell, &
                                             cutoff_squared=cutoff2, minimum_image=.FALSE.)
            IF (.NOT. ANY(primary(:, 1:nao) /= 0.0_dp)) CYCLE
            iimage = iimage + 1
            images(iimage)%basis => basis
            images(iimage)%atom = iatom
            images(iimage)%cell = bs_env%index_to_cell_3c(:, icell)
            images(iimage)%center = center
            ALLOCATE (images(iimage)%primary(bs_env%ri_rs%n_grid_points, nao))
            images(iimage)%primary(:, :) = primary(:, 1:nao)
            CALL pack_nonzero_ao_rows(images(iimage))
            CALL bound_ao_image(images(iimage))
         END DO
      END DO

   END SUBROUTINE prepare_ao_images

! **************************************************************************************************
!> \brief Store only exactly nonzero l rows of Eq. (1) when at least half the rows are zero.
!> primary_l maps packed rows back to the original r_l; no finite AO coefficient is discarded.
!> \param image ...
! **************************************************************************************************
   SUBROUTINE pack_nonzero_ao_rows(image)
      TYPE(ao_image_type), INTENT(INOUT)                 :: image

      INTEGER                                            :: il, ipacked, nkeep, nl
      LOGICAL, ALLOCATABLE                               :: nonzero(:)
      REAL(dp), ALLOCATABLE                              :: primary(:, :)

      nl = SIZE(image%primary, 1)
      ALLOCATE (nonzero(nl))
      DO il = 1, nl
         nonzero(il) = ANY(ABS(image%primary(il, :)) > 0.0_dp)
      END DO
      nkeep = COUNT(nonzero)
      IF (nkeep == 0 .OR. nkeep > nl/2) RETURN
      ALLOCATE (primary(nkeep, SIZE(image%primary, 2)), image%primary_l(nkeep))
      ipacked = 0
      DO il = 1, nl
         IF (.NOT. nonzero(il)) CYCLE
         ipacked = ipacked + 1
         primary(ipacked, :) = image%primary(il, :)
         image%primary_l(ipacked) = il
      END DO
      CALL MOVE_ALLOC(primary, image%primary)

   END SUBROUTINE pack_nonzero_ao_rows

! **************************************************************************************************
!> \brief Bound each image contribution to f_l(r) in Eq. (1) by B exp(-a_min r^2/2).
!> |x^n exp(-a x^2/2)| <= (n/(a e))^(n/2) bounds the remaining polynomial and Gaussian.
!> \param image ...
! **************************************************************************************************
   SUBROUTINE bound_ao_image(image)
      TYPE(ao_image_type), INTENT(INOUT)                 :: image

      INTEGER                                            :: ico, idir, il, ipgf, irow, iset, isgf, &
                                                            ncart, power
      REAL(dp)                                           :: amplitude, exponent, factor
      REAL(dp), ALLOCATABLE                              :: ao_bound(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      basis => image%basis
      image%decay = HUGE(0.0_dp)
      ALLOCATE (ao_bound(SIZE(image%primary, 2)), source=0.0_dp)
      DO iset = 1, basis%nset
         ncart = ncoset(basis%lmax(iset))
         DO ipgf = 1, basis%npgf(iset)
            exponent = basis%zet(ipgf, iset)
            CPASSERT(exponent > 0.0_dp)
            image%decay = MIN(image%decay, 0.5_dp*exponent)
            DO ico = 1, ncart
               factor = 1.0_dp
               DO idir = 1, 3
                  power = indco(idir, ico)
                  IF (power > 0) factor = factor*(REAL(power, dp)/(exponent*EXP(1.0_dp)))** &
                                          (0.5_dp*REAL(power, dp))
               END DO
               irow = (ipgf - 1)*ncart + ico
               DO isgf = basis%first_sgf(1, iset), basis%last_sgf(basis%nshell(iset), iset)
                  ao_bound(isgf) = ao_bound(isgf) + ABS(basis%sphi(irow, isgf))*factor
               END DO
            END DO
         END DO
      END DO
      image%bound = 0.0_dp
      DO il = 1, SIZE(image%primary, 1)
         amplitude = 0.0_dp
         DO isgf = 1, SIZE(ao_bound)
            amplitude = amplitude + ABS(image%primary(il, isgf))*ao_bound(isgf)
         END DO
         image%bound = MAX(image%bound, amplitude)
      END DO
      CALL prepare_ao_radial_envelope(image)

   END SUBROUTINE bound_ao_image

! **************************************************************************************************
!> \brief Bound each Cartesian primitive in Eq. (1) by its exact angular maximum at fixed radius.
!> For degree n, |x^i y^j z^k| <= r^n product_d (n_d/n)^(n_d/2).
!> Positive coefficient sums bound every l row without changing the field error budget.
!> \param image ...
! **************************************************************************************************
   SUBROUTINE prepare_ao_radial_envelope(image)
      TYPE(ao_image_type), INTENT(INOUT)                 :: image

      INTEGER                                            :: ico, idegree, idir, iexponent, il, ipgf, &
                                                            irow, iset, isgf, ncart, power
      REAL(dp)                                           :: amplitude, angular_factor
      REAL(dp), ALLOCATABLE                              :: ao_weight(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis

      basis => image%basis
      ALLOCATE (image%envelope_exponents(SUM(basis%npgf(:basis%nset))), &
                image%radial_envelope(0:MAXVAL(basis%lmax(:basis%nset)), SUM(basis%npgf(:basis%nset))), &
                ao_weight(SIZE(image%primary, 2)))
      ao_weight = 0.0_dp
      DO isgf = 1, SIZE(ao_weight)
         DO il = 1, SIZE(image%primary, 1)
            ao_weight(isgf) = MAX(ao_weight(isgf), ABS(image%primary(il, isgf)))
         END DO
      END DO
      image%radial_envelope = 0.0_dp
      iexponent = 0
      DO iset = 1, basis%nset
         ncart = ncoset(basis%lmax(iset))
         DO ipgf = 1, basis%npgf(iset)
            iexponent = iexponent + 1
            image%envelope_exponents(iexponent) = basis%zet(ipgf, iset)
            DO ico = 1, ncart
               idegree = SUM(indco(:, ico))
               angular_factor = 1.0_dp
               DO idir = 1, 3
                  power = indco(idir, ico)
                  IF (power > 0) angular_factor = angular_factor* &
                                                  (REAL(power, dp)/REAL(idegree, dp))**(0.5_dp*REAL(power, dp))
               END DO
               irow = (ipgf - 1)*ncart + ico
               amplitude = 0.0_dp
               DO isgf = basis%first_sgf(1, iset), basis%last_sgf(basis%nshell(iset), iset)
                  amplitude = amplitude + ao_weight(isgf)*ABS(basis%sphi(irow, isgf))
               END DO
               image%radial_envelope(idegree, iexponent) = &
                  image%radial_envelope(idegree, iexponent) + angular_factor*amplitude
            END DO
         END DO
      END DO
      image%radial_envelope(:, :) = image%radial_envelope* &
                                    (1.0_dp + EPSILON(1.0_dp)*(64.0_dp + &
                                                            4.0_dp*REAL(SIZE(image%radial_envelope), dp)*REAL(SIZE(ao_weight), dp)))

   END SUBROUTINE prepare_ao_radial_envelope

! **************************************************************************************************
!> \brief Bound Eq. (1) outside a sphere using max_{r>=r_min} r^n exp(-alpha*r²).
!> The maximum is at r²=max(r_min²,n/(2*alpha)); far-field polynomials share one exponential.
!> \param image ...
!> \param distance_squared ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION bound_ao_image_distance(image, distance_squared) RESULT(bound)
      TYPE(ao_image_type), INTENT(IN)                    :: image
      REAL(dp), INTENT(IN)                               :: distance_squared
      REAL(dp)                                           :: bound

      INTEGER                                            :: idegree, iexponent, lmax
      REAL(dp)                                           :: exponent, polynomial, radius, &
                                                            radius_squared

      bound = 0.0_dp
      lmax = UBOUND(image%radial_envelope, 1)
      radius = SQRT(distance_squared)
      DO iexponent = 1, SIZE(image%envelope_exponents)
         exponent = image%envelope_exponents(iexponent)
         IF (exponent*distance_squared >= 0.5_dp*REAL(lmax, dp)) THEN
            polynomial = image%radial_envelope(lmax, iexponent)
            DO idegree = lmax - 1, 0, -1
               polynomial = image%radial_envelope(idegree, iexponent) + radius*polynomial
            END DO
            bound = bound + polynomial*EXP(-exponent*distance_squared)
         ELSE
            DO idegree = 0, lmax
               radius_squared = MAX(distance_squared, 0.5_dp*REAL(idegree, dp)/exponent)
               bound = bound + image%radial_envelope(idegree, iexponent)* &
                       radius_squared**(0.5_dp*REAL(idegree, dp))*EXP(-exponent*radius_squared)
            END DO
         END IF
      END DO
   END FUNCTION bound_ao_image_distance

! **************************************************************************************************
!> \brief Evaluate Eq. (1), then square, on one quadrature block.
!> The AO image (A, R) on the lattice copy t of a block equals the image (A, R - t) on the class
!> block, so its values come from the class cache. The threads first fill the cache, one image
!> each, then sum Eq. (1) on separate point slices.
!> \param images ...
!> \param points ...
!> \param point_index index of each point in its complete quadrature block
!> \param copy lattice copy t of the quadrature block
!> \param cache AO values of the block class, shared by its copies
!> \param field_squared ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE evaluate_rho_l_squared(images, points, point_index, copy, cache, field_squared, tolerance)
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: point_index(:), copy(3)
      TYPE(block_class_cache_type), INTENT(INOUT)        :: cache
      REAL(dp), INTENT(OUT)                              :: field_squared(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance

      INTEGER                                            :: first, handle, ig, iimage, islot, last, &
                                                            nao, ng, nimages, nmissing, nslice, &
                                                            nthreads, row, u(3)
      INTEGER, ALLOCATABLE                               :: missing_indices(:)
      REAL(dp), ALLOCATABLE                              :: values(:, :)
      LOGICAL                                            :: image_needed(SIZE(images))
      INTEGER                                            :: image_slot(SIZE(images))

      ng = SIZE(points, 2)
      nimages = 0
      DO iimage = 1, SIZE(images)
         IF (.NOT. ALLOCATED(images(iimage)%primary)) EXIT
         nimages = nimages + 1
      END DO
      CALL screen_ao_images(images(1:nimages), points, tolerance, image_needed(1:nimages))
      IF (.NOT. ANY(image_needed(1:nimages))) THEN
         field_squared = 0.0_dp
         RETURN
      END IF
      image_slot = 0
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(images, points, point_index, copy, cache, ng, nimages, image_needed, image_slot) &
      !$OMP PRIVATE(handle, ig, iimage, islot, nao, nmissing, row, u, missing_indices, values)
      ALLOCATE (missing_indices(ng))
      !$OMP DO SCHEDULE(DYNAMIC)
      DO iimage = 1, nimages
         IF (.NOT. image_needed(iimage)) CYCLE
         u = images(iimage)%cell - copy
         ! Only this thread handles the image, so only the slot counter is shared.
         islot = cache%slot(u(1), u(2), u(3), images(iimage)%atom)
         IF (islot == 0) THEN
            !$OMP CRITICAL(d_lP_ao_cache)
            IF (cache%n_used < cache%n_slots) THEN
               cache%n_used = cache%n_used + 1
               islot = cache%n_used
               cache%keys(:, islot) = [u, images(iimage)%atom]
            END IF
            !$OMP END CRITICAL(d_lP_ao_cache)
            IF (islot == 0) CYCLE
            cache%slot(u(1), u(2), u(3), images(iimage)%atom) = islot
            cache%evaluated(:, islot) = .FALSE.
         END IF
         image_slot(iimage) = islot
         nmissing = 0
         DO ig = 1, ng
            IF (.NOT. cache%evaluated(point_index(ig), islot)) THEN
               nmissing = nmissing + 1
               missing_indices(nmissing) = ig
            END IF
         END DO
         IF (nmissing == 0) CYCLE
         CALL timeset('d_lP_AO_values', handle)
         nao = SIZE(images(iimage)%primary, 2)
         ALLOCATE (values(nmissing, nao))
         CALL evaluate_mu_r_values(values, points(:, missing_indices(1:nmissing)), &
                                   images(iimage)%basis, images(iimage)%center)
         DO row = 1, nmissing
            cache%values(1:nao, point_index(missing_indices(row)), islot) = values(row, :)
            cache%evaluated(point_index(missing_indices(row)), islot) = .TRUE.
         END DO
         DEALLOCATE (values)
         CALL timestop(handle)
      END DO
      !$OMP END DO
      DEALLOCATE (missing_indices)
      !$OMP END PARALLEL
      nthreads = 1
!$    nthreads = omp_get_max_threads()
      nslice = MAX(1, (ng + nthreads - 1)/nthreads)
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(images, points, point_index, cache, field_squared, ng, nslice, nimages, image_slot, tolerance) &
      !$OMP PRIVATE(first, last)
      DO first = 1, ng, nslice
         last = MIN(ng, first + nslice - 1)
         CALL evaluate_rho_l_squared_slice(images(1:nimages), points(:, first:last), point_index(first:last), &
                                           image_slot(1:nimages), cache, field_squared(first:last, :), tolerance)
      END DO
      !$OMP END PARALLEL DO

   END SUBROUTINE evaluate_rho_l_squared

! **************************************************************************************************
!> \brief Select the AO images a set of points needs for the squared-field error bound.
!> For image bounds b_i and A=sum_i b_i, omit at most delta=sqrt(A²+epsilon)-A in amplitude;
!> then 2A*delta+delta²<=epsilon. Points with A²<=epsilon need no image. A nonpositive tolerance
!> keeps every image. A subset of the points needs a subset of the images.
!> \param images ...
!> \param points ...
!> \param tolerance ...
!> \param image_needed ...
! **************************************************************************************************
   SUBROUTINE screen_ao_images(images, points, tolerance, image_needed)
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      REAL(dp), INTENT(IN)                               :: points(:, :), tolerance
      LOGICAL, INTENT(OUT)                               :: image_needed(:)

      INTEGER                                            :: iimage, nimages
      REAL(dp)                                           :: amplitude, delta, lower(3), &
                                                            root_tolerance, upper(3), &
                                                            image_bounds(SIZE(images))

      image_needed = .TRUE.
      IF (tolerance <= 0.0_dp) RETURN
      nimages = SIZE(images)
      lower = MINVAL(points, DIM=2)
      upper = MAXVAL(points, DIM=2)
      DO iimage = 1, nimages
         image_bounds(iimage) = images(iimage)%bound* &
                                EXP(-images(iimage)%decay*distance_to_box_squared(images(iimage)%center, lower, upper))
      END DO
      amplitude = SUM(image_bounds)*(1.0_dp + REAL(nimages, dp)*EPSILON(1.0_dp))
      root_tolerance = SQRT(tolerance)
      IF (amplitude <= root_tolerance) THEN
         image_needed = .FALSE.
         RETURN
      END IF
      delta = (tolerance/amplitude)/(SQRT(1.0_dp + (root_tolerance/amplitude)**2) + 1.0_dp)
      delta = delta/REAL(MAX(1, nimages), dp)
      DO iimage = 1, nimages
         IF (image_bounds(iimage) <= delta) THEN
            image_needed(iimage) = .FALSE.
         ELSE
            image_needed(iimage) = bound_ao_image_distance(images(iimage), &
                                                           distance_to_box_squared(images(iimage)%center, lower, upper)) > delta
         END IF
      END DO

   END SUBROUTINE screen_ao_images

! **************************************************************************************************
!> \brief Evaluate Eq. (1), then square, on a slice of a quadrature block.
!> Images with a cache slot take their values from the class cache, the others are evaluated.
!> \param images ...
!> \param points ...
!> \param point_index index of each point in its complete quadrature block
!> \param image_slot cache slot of each image, or 0
!> \param cache AO values of the block class
!> \param field_squared ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE evaluate_rho_l_squared_slice(images, points, point_index, image_slot, cache, field_squared, &
                                           tolerance)
      TYPE(ao_image_type), INTENT(IN)                    :: images(:)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: point_index(:), image_slot(:)
      TYPE(block_class_cache_type), INTENT(IN)           :: cache
      REAL(dp), INTENT(OUT)                              :: field_squared(:, :)
      REAL(dp), INTENT(IN)                               :: tolerance

      INTEGER                                            :: handle, ig, iimage, il, nao, ng, nl, &
                                                            nprimary
      REAL(dp), ALLOCATABLE                              :: field(:, :), image_field(:, :), &
                                                            values(:, :)
      LOGICAL                                            :: image_needed(SIZE(images))

      ng = SIZE(points, 2)
      nl = SIZE(field_squared, 2)
      CALL screen_ao_images(images, points, tolerance, image_needed)
      ALLOCATE (field(ng, nl), SOURCE=0.0_dp)
      DO iimage = 1, SIZE(images)
         IF (.NOT. image_needed(iimage)) CYCLE
         nao = SIZE(images(iimage)%primary, 2)
         ALLOCATE (values(ng, nao))
         CALL timeset('d_lP_AO_values', handle)
         IF (image_slot(iimage) > 0) THEN
            DO ig = 1, ng
               values(ig, :) = cache%values(1:nao, point_index(ig), image_slot(iimage))
            END DO
         ELSE
            CALL evaluate_mu_r_values(values, points, images(iimage)%basis, images(iimage)%center)
         END IF
         CALL timestop(handle)
         CALL timeset('d_lP_f_l_sum_AO', handle)
         nprimary = SIZE(images(iimage)%primary, 1)
         IF (ALLOCATED(images(iimage)%primary_l)) THEN
            IF (.NOT. ALLOCATED(image_field)) ALLOCATE (image_field(ng, nl))
            CALL dgemm('N', 'T', ng, nprimary, nao, 1.0_dp, values, ng, &
                       images(iimage)%primary, nprimary, 0.0_dp, image_field, ng)
            DO il = 1, nprimary
               field(:, images(iimage)%primary_l(il)) = &
                  field(:, images(iimage)%primary_l(il)) + image_field(:, il)
            END DO
         ELSE
            CALL dgemm('N', 'T', ng, nl, nao, 1.0_dp, values, ng, &
                       images(iimage)%primary, nl, 1.0_dp, field, ng)
         END IF
         CALL timestop(handle)
         DEALLOCATE (values)
      END DO
      field_squared = field**2

   END SUBROUTINE evaluate_rho_l_squared_slice

! **************************************************************************************************
!> \brief Evaluate the AO values in Eq. (1), sharing primitive factors before contraction.
!> Columns of primitive are Cartesian Gaussian primitives, grouped by basis set and exponent.
!> The original spherical contraction coefficients define exactly the same unscreened AO values.
!> \param values ...
!> \param points ...
!> \param basis ...
!> \param center ...
! **************************************************************************************************
   SUBROUTINE evaluate_mu_r_values(values, points, basis, center)
      REAL(dp), INTENT(OUT)                              :: values(:, :)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: center(3)

      INTEGER                                            :: first_sgf, ico, ig, ipgf, ipower, iset, &
                                                            last_sgf, ncart, ng, nprimitive
      INTEGER, ALLOCATABLE                               :: offsets(:)
      LOGICAL                                            :: parallel_points
      REAL(dp)                                           :: relative(3)
      REAL(dp), ALLOCATABLE                              :: coordinate_power(:, :), exponential(:), &
                                                            monomial(:, :), primitive(:, :), &
                                                            radius2(:)

      ng = SIZE(points, 2)
      IF (ng == 0) RETURN
      ALLOCATE (offsets(basis%nset), monomial(ng, ncoset(MAXVAL(basis%lmax))), &
                coordinate_power(3, 0:MAXVAL(basis%lmax)), radius2(ng), exponential(ng))
      nprimitive = 0
      DO iset = 1, basis%nset
         offsets(iset) = nprimitive
         nprimitive = nprimitive + basis%npgf(iset)*ncoset(basis%lmax(iset))
      END DO
      ALLOCATE (primitive(ng, nprimitive))
      parallel_points = .TRUE.
!$    parallel_points = .NOT. omp_in_parallel()
      !$OMP PARALLEL DEFAULT(NONE) IF(parallel_points) &
      !$OMP SHARED(ng, points, center, basis, primitive, offsets, ncoset, indco, radius2, monomial) &
      !$OMP PRIVATE(ig, relative, ico, coordinate_power, ipower, iset, ipgf, ncart, exponential)
      !$OMP DO SCHEDULE(STATIC)
      DO ig = 1, ng
         relative = points(:, ig) - center
         radius2(ig) = DOT_PRODUCT(relative, relative)
         coordinate_power(:, 0) = 1.0_dp
         DO ipower = 1, UBOUND(coordinate_power, 2)
            coordinate_power(:, ipower) = coordinate_power(:, ipower - 1)*relative
         END DO
         DO ico = 1, SIZE(monomial, 2)
            monomial(ig, ico) = coordinate_power(1, indco(1, ico))* &
                                coordinate_power(2, indco(2, ico))*coordinate_power(3, indco(3, ico))
         END DO
      END DO
      !$OMP END DO
      ! Fill primitive columns contiguously in g, reusing each Cartesian monomial.
      DO iset = 1, basis%nset
         ncart = ncoset(basis%lmax(iset))
         !$OMP DO SCHEDULE(STATIC)
         DO ipgf = 1, basis%npgf(iset)
            !$OMP SIMD
            DO ig = 1, ng
               exponential(ig) = EXP(-basis%zet(ipgf, iset)*radius2(ig))
            END DO
            DO ico = 1, ncart
               !$OMP SIMD
               DO ig = 1, ng
                  primitive(ig, offsets(iset) + (ipgf - 1)*ncart + ico) = monomial(ig, ico)*exponential(ig)
               END DO
            END DO
         END DO
         !$OMP END DO
      END DO
      !$OMP END PARALLEL
      values = 0.0_dp
      DO iset = 1, basis%nset
         nprimitive = basis%npgf(iset)*ncoset(basis%lmax(iset))
         first_sgf = basis%first_sgf(1, iset)
         last_sgf = basis%last_sgf(basis%nshell(iset), iset)
         CALL dgemm('N', 'N', ng, last_sgf - first_sgf + 1, nprimitive, &
                    1.0_dp, primitive(:, offsets(iset) + 1:), ng, &
                    basis%sphi(:, first_sgf:last_sgf), SIZE(basis%sphi, 1), &
                    1.0_dp, values(:, first_sgf:last_sgf), ng)
      END DO

   END SUBROUTINE evaluate_mu_r_values

! **************************************************************************************************
!> \brief Bound every auxiliary potential in Eq. (2), including contracted angular functions.
!> |x^n exp(-a x²/2)| <= (n/(a e))^(n/2); the remaining Gaussian potential is <=4 pi/a.
!> The triangle inequality uses the contracted primitive coefficients of each function.
!> Truncating the positive kernel preserves the bound; envelopes are indexed (bound/exponents,atom).
!> \param bs_env ...
!> \param aux_potentials ...
!> \param envelopes ...
! **************************************************************************************************
   SUBROUTINE auxiliary_potential_bound(bs_env, aux_potentials, envelopes)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(radial_interpolation_type), ALLOCATABLE, &
         INTENT(IN)                                      :: aux_potentials(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: envelopes(:, :)

      INTEGER                                            :: iatom, ico, idir, ikind, ipgf, &
                                                            irow_packed, iset, ishell, l, power
      LOGICAL, ALLOCATABLE                               :: angular_used(:)
      REAL(dp)                                           :: exponent, factor
      REAL(dp), ALLOCATABLE                              :: bounds(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)

      particles => bs_env%ri_rs%particle_set
      ALLOCATE (envelopes(3, bs_env%n_atom))
      envelopes = 0.0_dp
      DO iatom = 1, bs_env%n_atom
         ikind = particles(iatom)%atomic_kind%kind_number
         basis => bs_env%basis_set_RI(ikind)%gto_basis_set
         CPASSERT(ASSOCIATED(basis))
         ALLOCATE (bounds(SIZE(aux_potentials(iatom)%coefficients, 2)))
         ALLOCATE (angular_used(0:MAXVAL(basis%lmax)))
         bounds = 0.0_dp
         irow_packed = 0
         envelopes(2, iatom) = HUGE(1.0_dp)
         DO iset = 1, basis%nset
            angular_used = .FALSE.
            DO ishell = 1, basis%nshell(iset)
               angular_used(basis%l(ishell, iset)) = .TRUE.
            END DO
            DO ipgf = 1, basis%npgf(iset)
               exponent = basis%zet(ipgf, iset)
               envelopes(2, iatom) = MIN(envelopes(2, iatom), exponent)
               envelopes(3, iatom) = MAX(envelopes(3, iatom), exponent)
               DO l = 0, basis%lmax(iset)
                  IF (.NOT. angular_used(l)) CYCLE
                  DO ico = ncoset(l - 1) + 1, ncoset(l)
                     irow_packed = irow_packed + 1
                     factor = 4.0_dp*pi/exponent
                     DO idir = 1, 3
                        power = indco(idir, ico)
                        IF (power > 0) factor = factor*(REAL(power, dp)/(exponent*EXP(1.0_dp)))** &
                                                (0.5_dp*REAL(power, dp))
                     END DO
                     bounds(:) = bounds(:) + ABS(aux_potentials(iatom)%coefficients(irow_packed, :))*factor
                  END DO
               END DO
            END DO
         END DO
         envelopes(1, iatom) = MAXVAL(bounds)
         DEALLOCATE (bounds, angular_used)
      END DO

   END SUBROUTINE auxiliary_potential_bound

! **************************************************************************************************
!> \brief Contract Eqs. (2),(4), screening every spatial l/P-atom/T block separately.
!> B bounds every selected potential on this atom; a_min/max bound its primitive exponents.
!> For erfc(ωr)/r, |v_P(r)|<=B exp[-a_min ω² |r-center|²/(a_min+2ω²)].
!> This follows from erfc(x)<=exp(-x²), completing the square, and the Gaussian Coulomb bound.
!> The tabulated radial bounds of prepare_short_potential_bounds tighten it where smaller.
!> Half the absolute budget screens common points. The other half is applied independently to
!> every d_lP(T): a point is contracted for one l block only when its local product bound exceeds
!> the per-point threshold. The exact l/P/T distance criterion is applied before each GEMM.
!> \param bs_env ...
!> \param cells_R_d_lP Retained integer cell vectors R, indexed (3,R).
!> \param include_lP_R ...
!> \param include_l_blocks_P_R ...
!> \param points ...
!> \param point_index index of each point in its complete quadrature block
!> \param field_squared ...
!> \param weight ...
!> \param envelopes ...
!> \param blocks ...
!> \param threshold ...
!> \param copy lattice copy t of the quadrature block
!> \param cache potentials of the block class, shared by its copies
!> \param potential_points ...
!> \param contraction_products ...
!> \param d_lP_R Numerical integrals indexed by local point l, RI function P and cell R.
!> \param aux_potentials ...
! **************************************************************************************************
   SUBROUTINE integrate_potentials(bs_env, cells_R_d_lP, include_lP_R, include_l_blocks_P_R, &
                                   points, point_index, field_squared, weight, envelopes, blocks, &
                                   threshold, copy, cache, potential_points, contraction_products, &
                                   d_lP_R, aux_potentials)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells_R_d_lP(:, :)
      LOGICAL, INTENT(IN)                                :: include_lP_R(:, :, :), &
                                                            include_l_blocks_P_R(:, :, :)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      INTEGER, INTENT(IN)                                :: point_index(:)
      REAL(dp), INTENT(IN)                               :: field_squared(:, :), weight, &
                                                            envelopes(:, :)
      TYPE(l_block_partition_type), INTENT(IN)           :: blocks
      REAL(dp), INTENT(IN)                               :: threshold
      INTEGER, INTENT(IN)                                :: copy(3)
      TYPE(block_class_cache_type), INTENT(INOUT)        :: cache
      INTEGER(int_8), INTENT(INOUT)                      :: potential_points, contraction_products
      REAL(dp), INTENT(INOUT)                            :: d_lP_R(:, :, :)
      TYPE(radial_interpolation_type), ALLOCATABLE, &
         INTENT(IN)                                      :: aux_potentials(:)

      INTEGER :: firstp, handle, iatom, iblock, ig, ikind, ilocal, iradial, islot, it, lastp, &
         maxnp, nactive, nmissing, np, nthreads, nunion, row, u(3)
      REAL(dp)                                           :: center(3), decay, field_max, &
                                                            points_lower(3), points_upper(3)
      REAL(dp), ALLOCATABLE :: block_field_bound(:, :), block_RI(:, :), field_work(:, :), &
         maximum_field(:), point_bound(:), potential(:, :), potential_t(:, :)
      LOGICAL                                            :: point_needed(SIZE(points, 2))
      INTEGER :: union_indices(SIZE(points, 2)), missing_indices(SIZE(points, 2))
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      LOGICAL                                            :: block_needed(blocks%nblocks)
      INTEGER :: block_active_last(blocks%nblocks), block_active_first(blocks%nblocks)
      TYPE(particle_type), POINTER                       :: particles(:)
      INTEGER :: active_l(bs_env%ri_rs%n_grid_points)

      cell => bs_env%ri_rs%cell
      particles => bs_env%ri_rs%particle_set
      maxnp = MAXVAL(bs_env%i_RI_end_from_atom - bs_env%i_RI_start_from_atom + 1)

      CALL timeset('d_lP_allocate_field_bounds', handle)

      ALLOCATE (block_field_bound(SIZE(points, 2), blocks%nblocks))

      CALL timestop(handle)

      CALL timeset('d_lP_prepare_field_bounds', handle)

      DO iblock = 1, blocks%nblocks
         block_field_bound(:, iblock) = 0.0_dp
         DO ilocal = blocks%first(iblock), blocks%first(iblock + 1) - 1
            DO ig = 1, SIZE(points, 2)
               block_field_bound(ig, iblock) = MAX(block_field_bound(ig, iblock), &
                                                   ABS(field_squared(ig, blocks%indices(ilocal))))
            END DO
         END DO
      END DO

      field_max = MAXVAL(block_field_bound)
      points_lower = MINVAL(points, DIM=2)
      points_upper = MAXVAL(points, DIM=2)

      CALL timestop(handle)

      ! Bound the additional thread workspaces to 128 MiB per MPI rank.
      ! A single thread retains the original memory bound for large global l grids.
      nthreads = 1
!$    nthreads = omp_get_max_threads()
      nthreads = MIN(nthreads, MAX(1, INT(thread_workspace_bytes/ &
                                          (real_bytes*(REAL(bs_env%ri_rs%n_grid_points, dp)*(maxnp + SIZE(points, 2)) + &
                                                       2.0_dp*REAL(SIZE(points, 2), dp)*maxnp)))))
      ! Each atom/image owns disjoint P/T columns of Eq. (2).

      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP SHARED(bs_env, cells_R_d_lP, points, point_index, field_squared, weight, envelopes, &
      !$OMP blocks, threshold, copy, cache, d_lP_R, aux_potentials, cell, particles, &
      !$OMP maxnp, block_field_bound, field_max, points_lower, points_upper, include_lP_R, &
      !$OMP include_l_blocks_P_R) &
      !$OMP PRIVATE(firstp, handle, iatom, iblock, ig, ikind, ilocal, iradial, islot, it, lastp, &
      !$OMP nactive, nmissing, np, nunion, row, u, center, decay, basis, &
      !$OMP block_RI, field_work, maximum_field, point_bound, potential, potential_t, point_needed, &
      !$OMP union_indices, missing_indices, &
      !$OMP active_l, block_needed, block_active_last, block_active_first) &
      !$OMP REDUCTION(+:potential_points, contraction_products)

      CALL timeset('d_lP_allocate_thread_buffers', handle)

      ALLOCATE (block_RI(bs_env%ri_rs%n_grid_points, maxnp), &
                field_work(SIZE(points, 2), bs_env%ri_rs%n_grid_points), point_bound(SIZE(points, 2)), &
                maximum_field(SIZE(points, 2)), potential(SIZE(points, 2), maxnp), &
                potential_t(maxnp, SIZE(points, 2)))

      CALL timestop(handle)

      !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
      DO it = 1, SIZE(cells_R_d_lP, 2)
         DO iatom = 1, bs_env%n_atom
            ikind = particles(iatom)%atomic_kind%kind_number
            basis => bs_env%basis_set_RI(ikind)%gto_basis_set
            firstp = bs_env%i_RI_start_from_atom(iatom)
            lastp = bs_env%i_RI_end_from_atom(iatom)
            np = lastp - firstp + 1
            center = pbc(particles(iatom)%r, cell) + MATMUL(bs_env%hmat, REAL(cells_R_d_lP(:, it), dp))
            IF (.NOT. ANY(include_l_blocks_P_R(:, iatom, it))) CYCLE
            ! The block of the lattice copy T' sees V_P around R_A + hT at U = T - T'.
            u = cells_R_d_lP(:, it) - copy
            islot = cache%slot(u(1), u(2), u(3), iatom)
            CALL timeset('d_lP_screen_atom_image', handle)
            decay = envelopes(2, iatom)*bs_env%ri_metric%omega**2/ &
                    (envelopes(2, iatom) + 2.0_dp*bs_env%ri_metric%omega**2)
            ! No point bound exceeds the envelope at the block point nearest to the atom image.
            IF (field_max*envelopes(1, iatom)* &
                EXP(-decay*distance_to_box_squared(center, points_lower, points_upper)) <= threshold) THEN
               CALL timestop(handle)
               CYCLE
            END IF
            DO ig = 1, SIZE(points, 2)
               IF (islot > 0) THEN
                  IF (cache%point_bound(point_index(ig), islot) >= 0.0_dp) THEN
                     point_bound(ig) = cache%point_bound(point_index(ig), islot)
                     CYCLE
                  END IF
               END IF
               point_bound(ig) = envelopes(1, iatom)
               IF (decay > 0.0_dp) point_bound(ig) = point_bound(ig)* &
                                                     EXP(-decay*SUM((points(:, ig) - center)**2))
               IF (ALLOCATED(aux_potentials(iatom)%short_bounds)) THEN
                  iradial = INT(NORM2(points(:, ig) - center)/short_bound_spacing)
                  IF (iradial < n_short_bounds) THEN
                     point_bound(ig) = MIN(point_bound(ig), aux_potentials(iatom)%short_bounds(iradial))
                  END IF
               END IF
               IF (islot > 0) cache%point_bound(point_index(ig), islot) = point_bound(ig)
            END DO
            ! A point is needed if it passes for one needed l block, i.e. for their largest field.
            block_needed = .FALSE.
            maximum_field = 0.0_dp
            DO iblock = 1, blocks%nblocks
               IF (.NOT. include_l_blocks_P_R(iblock, iatom, it)) CYCLE
               block_needed(iblock) = ANY(block_field_bound(:, iblock)*point_bound > threshold)
               IF (block_needed(iblock)) maximum_field(:) = MAX(maximum_field, block_field_bound(:, iblock))
            END DO
            point_needed = maximum_field*point_bound > threshold
            nunion = 0
            DO ig = 1, SIZE(points, 2)
               IF (point_needed(ig)) THEN
                  nunion = nunion + 1
                  union_indices(nunion) = ig
               END IF
            END DO
            CALL timestop(handle)
            IF (nunion == 0) CYCLE
            IF (islot == 0) THEN
               !$OMP CRITICAL(d_lP_potential_cache)
               IF (cache%n_used < cache%n_slots) THEN
                  cache%n_used = cache%n_used + 1
                  islot = cache%n_used
                  cache%keys(:, islot) = [u, iatom]
                  cache%slot(u(1), u(2), u(3), iatom) = islot
               END IF
               !$OMP END CRITICAL(d_lP_potential_cache)
               IF (islot > 0) THEN
                  cache%evaluated(:, islot) = .FALSE.
                  cache%point_bound(:, islot) = -1.0_dp
                  cache%point_bound(point_index, islot) = point_bound(1:SIZE(points, 2))
               END IF
            END IF
            IF (islot > 0) THEN
               ! Evaluate only the points that no earlier copy of the block needed.
               nmissing = 0
               DO row = 1, nunion
                  IF (.NOT. cache%evaluated(point_index(union_indices(row)), islot)) THEN
                     nmissing = nmissing + 1
                     missing_indices(nmissing) = union_indices(row)
                  END IF
               END DO
               IF (nmissing > 0) THEN
                  CALL timeset('d_lP_potential_eval', handle)
                  CALL contract_selected_potentials(potential(:, 1:np), &
                                                    points(:, missing_indices(1:nmissing)), &
                                                    basis, center, aux_potentials(iatom))
                  CALL timestop(handle)
                  DO row = 1, nmissing
                     cache%values(1:np, point_index(missing_indices(row)), islot) = potential(row, 1:np)
                     cache%evaluated(point_index(missing_indices(row)), islot) = .TRUE.
                  END DO
               END IF
               DO row = 1, nunion
                  potential_t(1:np, row) = cache%values(1:np, point_index(union_indices(row)), islot)
               END DO
               potential_points = potential_points + INT(nmissing, int_8)
            ELSE
               CALL timeset('d_lP_potential_eval', handle)
               CALL contract_selected_potentials(potential(:, 1:np), &
                                                 points(:, union_indices(1:nunion)), &
                                                 basis, center, aux_potentials(iatom))
               CALL timestop(handle)
               potential_t(1:np, 1:nunion) = TRANSPOSE(potential(1:nunion, 1:np))
               potential_points = potential_points + INT(nunion, int_8)
            END IF
            CALL timeset('d_lP_prepare_quadrature', handle)
            nactive = 0
            block_active_first = 0
            block_active_last = -1
            DO iblock = 1, blocks%nblocks
               IF (.NOT. block_needed(iblock)) CYCLE
               block_active_first(iblock) = nactive + 1
               DO ilocal = blocks%first(iblock), blocks%first(iblock + 1) - 1
                  IF (include_lP_R(blocks%indices(ilocal), iatom, it)) THEN
                     nactive = nactive + 1
                     active_l(nactive) = blocks%indices(ilocal)
                  END IF
               END DO
               block_active_last(iblock) = nactive
            END DO
            DO iblock = 1, blocks%nblocks
               IF (.NOT. block_needed(iblock)) CYCLE
               DO row = 1, nunion
                  ig = union_indices(row)
                  point_needed(row) = block_field_bound(ig, iblock)*point_bound(ig) > threshold
               END DO
               DO ilocal = block_active_first(iblock), block_active_last(iblock)
                  DO row = 1, nunion
                     field_work(row, ilocal) = MERGE( &
                                               field_squared(union_indices(row), active_l(ilocal)), 0.0_dp, point_needed(row))
                  END DO
               END DO
            END DO
            CALL timestop(handle)
            CALL timeset('d_lP_sum_g_f2_v', handle)
            CALL dgemm('T', 'T', nactive, np, nunion, weight, field_work, SIZE(field_work, 1), &
                       potential_t, SIZE(potential_t, 1), 0.0_dp, block_RI, SIZE(block_RI, 1))
            DO ilocal = 1, nactive
               d_lP_R(active_l(ilocal), firstp:lastp, it) = &
                  d_lP_R(active_l(ilocal), firstp:lastp, it) + block_RI(ilocal, 1:np)
            END DO
            CALL timestop(handle)
            contraction_products = contraction_products + &
                                   INT(nactive, int_8)*INT(np, int_8)*INT(nunion, int_8)
         END DO
      END DO
      !$OMP END DO NOWAIT

      CALL timeset('d_lP_wait_for_threads', handle)

      !$OMP BARRIER

      CALL timestop(handle)

      CALL timeset('d_lP_release_thread_buffers', handle)

      DEALLOCATE (block_RI, field_work, maximum_field, point_bound, potential, potential_t)

      CALL timestop(handle)

      !$OMP END PARALLEL

   END SUBROUTINE integrate_potentials

! **************************************************************************************************
!> \brief Squared distance from a point to an axis-aligned spatial l block.
!> \param point ...
!> \param lower ...
!> \param upper ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION distance_to_box_squared(point, lower, upper) RESULT(distance_squared)
      REAL(dp), INTENT(IN)                               :: point(3), lower(3), upper(3)
      REAL(dp)                                           :: distance_squared

      INTEGER                                            :: idir
      REAL(dp)                                           :: delta

      distance_squared = 0.0_dp
      DO idir = 1, 3
         delta = MAX(lower(idir) - point(idir), point(idir) - upper(idir), 0.0_dp)
         distance_squared = distance_squared + delta**2
      END DO
   END FUNCTION distance_to_box_squared

! **************************************************************************************************
!> \brief Transform spherical coefficients once; bound the radial tables by radial_table_bytes.
!> \param basis ...
!> \param U_Pp ...
!> \param contracted ...
!> \param attenuation ...
!> \param radial_table_bytes ...
!> \param mesh_step ...
!> \param mesh_displacement ...
! **************************************************************************************************
   SUBROUTINE prepare_contracted_potential(basis, U_Pp, contracted, attenuation, radial_table_bytes, &
                                           mesh_step, mesh_displacement)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: U_Pp(:, :)
      TYPE(radial_interpolation_type), INTENT(OUT)       :: contracted
      REAL(dp), INTENT(IN)                               :: attenuation, radial_table_bytes, &
                                                            mesh_step(3, 3), mesh_displacement(3)

      INTEGER                                            :: ico, iopt, ipgf, irow, irow_packed, &
                                                            iset, isgf, ishell, l, ncart, &
                                                            ninterval, np, nrows
      LOGICAL, ALLOCATABLE                               :: angular_used(:)
      REAL(dp)                                           :: bytes_per_node, table_bytes
      REAL(dp), ALLOCATABLE                              :: full_abscissae(:), moments(:, :)

      ALLOCATE (angular_used(0:MAXVAL(basis%lmax)))
      CPASSERT(SIZE(U_Pp, 1) == basis%nsgf)
      np = SIZE(U_Pp, 2)
      nrows = 0
      DO iset = 1, basis%nset
         angular_used = .FALSE.
         DO ishell = 1, basis%nshell(iset)
            angular_used(basis%l(ishell, iset)) = .TRUE.
         END DO
         DO l = 0, basis%lmax(iset)
            IF (angular_used(l)) nrows = nrows + basis%npgf(iset)*(ncoset(l) - ncoset(l - 1))
         END DO
      END DO
      ! W_g,ia is independent of the contracted shell: sum its coefficients before evaluating W.
      ! Rows enumerate (set, primitive, distinct angular momentum, Cartesian monomial).
      ALLOCATE (contracted%coefficients(nrows, np), source=0.0_dp)
      irow_packed = 0
      DO iset = 1, basis%nset
         ncart = ncoset(basis%lmax(iset))
         angular_used = .FALSE.
         DO ishell = 1, basis%nshell(iset)
            angular_used(basis%l(ishell, iset)) = .TRUE.
         END DO
         DO ipgf = 1, basis%npgf(iset)
            DO l = 0, basis%lmax(iset)
               IF (.NOT. angular_used(l)) CYCLE
               DO ico = ncoset(l - 1) + 1, ncoset(l)
                  irow = (ipgf - 1)*ncart + ico
                  irow_packed = irow_packed + 1
                  DO ishell = 1, basis%nshell(iset)
                     IF (basis%l(ishell, iset) /= l) CYCLE
                     DO isgf = basis%first_sgf(ishell, iset), basis%last_sgf(ishell, iset)
                        DO iopt = 1, np
                           contracted%coefficients(irow_packed, iopt) = &
                              contracted%coefficients(irow_packed, iopt) + basis%sphi(irow, isgf)*U_Pp(isgf, iopt)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (angular_used)
      IF (basis%nsgf < nrows) CALL prepare_spherical_potential_factorization(basis, U_Pp, contracted)
      contracted%attenuation = MAX(0.0_dp, attenuation)
      IF (contracted%attenuation > 0.0_dp) THEN
         CALL prepare_short_potential_bounds(basis, contracted%coefficients, contracted%attenuation, &
                                             contracted%short_bounds)
      END IF
      CALL prepare_radial_series(basis, contracted%attenuation, contracted%series, contracted%short_series)
      table_bytes = MAX(0.0_dp, radial_table_bytes)
      ! Budget the actual logarithmic nodes, including one primitive preparation buffer.
      bytes_per_node = real_bytes*REAL(MAXVAL(basis%lmax) + 7, dp)* &
                       REAL(MAXVAL(basis%npgf)*basis%nset + 1, dp)
      CALL prepare_short_radial_grid(short_grid_xmax, contracted%short_abscissae)
      table_bytes = MAX(0.0_dp, table_bytes - real_bytes*REAL(SIZE(contracted%short_abscissae), dp))
      ninterval = INT(MIN(REAL(SIZE(contracted%short_abscissae), dp), &
                          table_bytes/bytes_per_node)) - 1
      IF (ninterval < 1) THEN
         DEALLOCATE (contracted%short_abscissae)
         RETURN
      END IF
      IF (ninterval < UBOUND(contracted%short_abscissae, 1)) THEN
         CALL MOVE_ALLOC(contracted%short_abscissae, full_abscissae)
         ALLOCATE (contracted%short_abscissae(0:ninterval))
         contracted%short_abscissae(:) = full_abscissae(0:ninterval)
         DEALLOCATE (full_abscissae)
      END IF
      ALLOCATE (contracted%short_moments(0:MAXVAL(basis%lmax) + 6, &
                                         0:ninterval, MAXVAL(basis%npgf), basis%nset))
      DO iset = 1, basis%nset
         DO ipgf = 1, basis%npgf(iset)
            CALL prepare_short_radial_table(basis%zet(ipgf, iset), contracted%attenuation, &
                                            MAXVAL(basis%lmax), moments, ninterval, abscissae=contracted%short_abscissae)
            contracted%short_moments(:, :, ipgf, iset) = moments
            DEALLOCATE (moments)
         END DO
      END DO
      ! Keep the generic radial lookup at all radii before using remaining storage for mesh nodes.
      table_bytes = MAX(0.0_dp, table_bytes - real_bytes*REAL(SIZE(contracted%short_moments, KIND=int_8), dp))
      CALL prepare_mesh_radial_moments(basis, contracted, mesh_step, mesh_displacement, table_bytes)

   END SUBROUTINE prepare_contracted_potential

! **************************************************************************************************
!> \brief Bound every contracted erfc potential over radial intervals without cancellation.
!> H_l(r) = (2 pi/a) int_sqrt(b)^1 t^(2l) exp(-a r^2 t^2) dt, b = omega^2/(a+omega^2), decreases with r.
!> On [r_lo, r_hi], |x^i y^j z^k| <= r_hi^l prod_d (n_d/l)^(n_d/2), l = i+j+k. Upper rectangle sums
!> at the maxima of the integrand bound H_l(r_lo); summing absolute contracted coefficients then
!> bounds every P and every direction. For each homogeneous polynomial, Cauchy also gives
!> |sum_alpha c_alpha r^alpha| <= r^l sqrt(sum_alpha c_alpha^2/multinomial(l, alpha)); the smaller
!> angular bound is taken. Points beyond the last interval keep the Gaussian envelope bound.
!> \param basis ...
!> \param coefficients contracted Cartesian primitive rows of the RI functions P
!> \param omega attenuation of the erfc kernel
!> \param radial_bounds bound on |V_P(r)| for all P and r in each interval
! **************************************************************************************************
   SUBROUTINE prepare_short_potential_bounds(basis, coefficients, omega, radial_bounds)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: coefficients(:, :), omega
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: radial_bounds(:)

      INTEGER, PARAMETER                                 :: n_rectangles = 32

      INTEGER                                            :: ico, idir, ifactor, il, ipgf, ipower, &
                                                            iq, ir, irow, iset, ishell, istart
      REAL(dp)                                           :: a, angular, factor, integral, lower, &
                                                            multinomial, peak, rhi, rlo, thi, tlo, &
                                                            width, x
      REAL(dp), DIMENSION(SIZE(coefficients, 2))         :: bounds, cscale, sumsq, triangle
      LOGICAL :: angular_used(0:MAXVAL(basis%lmax))

      ALLOCATE (radial_bounds(0:n_short_bounds - 1))
      DO ir = 0, n_short_bounds - 1
         rlo = REAL(ir, dp)*short_bound_spacing
         rhi = REAL(ir + 1, dp)*short_bound_spacing
         bounds = 0.0_dp
         irow = 0
         DO iset = 1, basis%nset
            angular_used = .FALSE.
            DO ishell = 1, basis%nshell(iset)
               angular_used(basis%l(ishell, iset)) = .TRUE.
            END DO
            DO ipgf = 1, basis%npgf(iset)
               a = basis%zet(ipgf, iset)
               lower = SQRT(omega**2/(a + omega**2))
               width = 1.0_dp - lower
               x = a*rlo**2
               DO il = 0, basis%lmax(iset)
                  IF (.NOT. angular_used(il)) CYCLE
                  integral = 0.0_dp
                  DO iq = 0, n_rectangles - 1
                     tlo = lower + width*REAL(iq, dp)/REAL(n_rectangles, dp)
                     thi = lower + width*REAL(iq + 1, dp)/REAL(n_rectangles, dp)
                     peak = thi
                     IF (x > 0.0_dp) peak = MAX(tlo, MIN(thi, SQRT(REAL(il, dp)/x)))
                     integral = integral + width/REAL(n_rectangles, dp)*peak**(2*il)*EXP(-x*peak**2)
                  END DO
                  factor = 2.0_dp*pi/a*integral*rhi**il
                  istart = irow
                  cscale = 0.0_dp
                  DO ico = 1, ncoset(il) - ncoset(il - 1)
                     cscale = MAX(cscale, ABS(coefficients(istart + ico, :)))
                  END DO
                  triangle = 0.0_dp
                  sumsq = 0.0_dp
                  DO ico = ncoset(il - 1) + 1, ncoset(il)
                     angular = 1.0_dp
                     multinomial = 1.0_dp
                     DO ifactor = 2, il
                        multinomial = multinomial*REAL(ifactor, dp)
                     END DO
                     DO idir = 1, 3
                        ipower = indco(idir, ico)
                        IF (ipower > 0) angular = angular*(REAL(ipower, dp)/REAL(il, dp))**(0.5_dp*REAL(ipower, dp))
                        DO ifactor = 2, ipower
                           multinomial = multinomial/REAL(ifactor, dp)
                        END DO
                     END DO
                     irow = irow + 1
                     triangle = triangle + ABS(coefficients(irow, :))*angular
                     sumsq = sumsq + (coefficients(irow, :)/MAX(TINY(1.0_dp), cscale))**2/multinomial
                  END DO
                  WHERE (cscale >= TINY(1.0_dp)) triangle = MIN(triangle, cscale*SQRT(sumsq))
                  bounds = bounds + triangle*MAX(TINY(1.0_dp), factor)
               END DO
            END DO
         END DO
         CPASSERT(irow == SIZE(coefficients, 1))
         ! Inflate by roundoff.
         radial_bounds(ir) = MAX(TINY(1.0_dp), MAXVAL(bounds)*(1.0_dp + eps_bound_roundoff))
      END DO

   END SUBROUTINE prepare_short_potential_bounds

! **************************************************************************************************
!> \brief Factor sphi into a common angular map and shell/primitive radial scales when exact.
!> This is checked against every spherical coefficient; arbitrary contracted or modified bases
!> retain the original Cartesian path if the factorization does not hold.
!> \param basis ...
!> \param U_Pp ...
!> \param contracted ...
! **************************************************************************************************
   SUBROUTINE prepare_spherical_potential_factorization(basis, U_Pp, contracted)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: U_Pp(:, :)
      TYPE(radial_interpolation_type), INTENT(INOUT)     :: contracted

      INTEGER                                            :: c, first, ico, ipgf, iset, ishell, l, &
                                                            lmax, m, ncart, ref_c, ref_m, ref_p, &
                                                            ref_shell, row
      LOGICAL, ALLOCATABLE                               :: used(:)
      REAL(dp)                                           :: actual, expected, ref_abs, ref_value, &
                                                            scale
      REAL(dp), ALLOCATABLE                              :: angular_coeff(:, :, :, :), &
                                                            primitive_scale(:, :, :)

      IF (SIZE(U_Pp, 1) /= basis%nsgf) RETURN
      ALLOCATE (used(basis%nsgf), source=.FALSE.)
      ALLOCATE (angular_coeff(MAXVAL(ncoset(0:MAXVAL(basis%lmax)) - &
                                     ncoset(-1:MAXVAL(basis%lmax) - 1)), &
                              MAXVAL(nso(0:MAXVAL(basis%lmax))), &
                              0:MAXVAL(basis%lmax), basis%nset), source=0.0_dp)
      ALLOCATE (primitive_scale(MAXVAL(basis%npgf), MAXVAL(basis%nshell), basis%nset), &
                source=0.0_dp)
      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            first = basis%first_sgf(ishell, iset)
            IF (first < 1 .OR. basis%last_sgf(ishell, iset) /= first + nso(l) - 1) RETURN
            IF (basis%last_sgf(ishell, iset) > basis%nsgf) RETURN
            IF (ANY(used(first:first + nso(l) - 1))) RETURN
            used(first:first + nso(l) - 1) = .TRUE.
         END DO
      END DO
      IF (.NOT. ALL(used)) RETURN
      DO iset = 1, basis%nset
         lmax = basis%lmax(iset)
         DO l = 0, lmax
            ref_abs = 0.0_dp
            ref_shell = 0
            ncart = ncoset(l) - ncoset(l - 1)
            DO ishell = 1, basis%nshell(iset)
               IF (basis%l(ishell, iset) /= l) CYCLE
               first = basis%first_sgf(ishell, iset)
               DO ipgf = 1, basis%npgf(iset)
                  DO m = 1, nso(l)
                     DO c = 1, ncart
                        ico = ncoset(l - 1) + c
                        row = (ipgf - 1)*ncoset(lmax) + ico
                        actual = basis%sphi(row, first + m - 1)
                        IF (ABS(actual) > ref_abs) THEN
                           ref_abs = ABS(actual)
                           ref_value = actual
                           ref_shell = ishell
                           ref_p = ipgf
                           ref_c = c
                           ref_m = m
                        END IF
                     END DO
                  END DO
               END DO
            END DO
            IF (ref_shell == 0) CYCLE
            first = basis%first_sgf(ref_shell, iset)
            DO m = 1, nso(l)
               DO c = 1, ncart
                  ico = ncoset(l - 1) + c
                  row = (ref_p - 1)*ncoset(lmax) + ico
                  angular_coeff(c, m, l, iset) = basis%sphi(row, first + m - 1)
               END DO
            END DO
            DO ishell = 1, basis%nshell(iset)
               IF (basis%l(ishell, iset) /= l) CYCLE
               first = basis%first_sgf(ishell, iset)
               DO ipgf = 1, basis%npgf(iset)
                  row = (ipgf - 1)*ncoset(lmax) + ncoset(l - 1) + ref_c
                  scale = basis%sphi(row, first + ref_m - 1)/ref_value
                  primitive_scale(ipgf, ishell, iset) = scale
                  DO m = 1, nso(l)
                     DO c = 1, ncart
                        ico = ncoset(l - 1) + c
                        row = (ipgf - 1)*ncoset(lmax) + ico
                        actual = basis%sphi(row, first + m - 1)
                        expected = scale*angular_coeff(c, m, l, iset)
                        IF (ABS(actual - expected) > 64.0_dp*EPSILON(1.0_dp)* &
                            MAX(ref_abs, ABS(actual), ABS(expected))) RETURN
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL MOVE_ALLOC(angular_coeff, contracted%spherical_angular)
      CALL MOVE_ALLOC(primitive_scale, contracted%spherical_primitive)
      ALLOCATE (contracted%spherical_selection(SIZE(U_Pp, 1), SIZE(U_Pp, 2)))
      contracted%spherical_selection(:, :) = U_Pp
      contracted%spherical_factorized = .TRUE.

   END SUBROUTINE prepare_spherical_potential_factorization

! **************************************************************************************************
!> \brief Prepare the origin series in Eq. (2) once for a fixed basis and erfc(omega*r)/r kernel.
!> The coefficients depend on neither the integration point nor the translated atomic center.
!> \param basis ...
!> \param omega ...
!> \param series ...
!> \param short_series ...
! **************************************************************************************************
   SUBROUTINE prepare_radial_series(basis, omega, series, short_series)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: omega
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: series(:, :, :), short_series(:, :, :)

      INTEGER                                            :: in, ipgf, iset, lmax
      REAL(dp)                                           :: beta_ratio, exponent

      lmax = MAXVAL(basis%lmax)
      ALLOCATE (series(0:lmax + 24, MAXVAL(basis%npgf), basis%nset), &
                short_series(0:lmax + 24, MAXVAL(basis%npgf), basis%nset))
      DO iset = 1, basis%nset
         DO ipgf = 1, basis%npgf(iset)
            exponent = basis%zet(ipgf, iset)
            CALL origin_series(series(:, ipgf, iset))
            beta_ratio = omega**2/(exponent + omega**2)
            DO in = 0, UBOUND(series, 1)
               short_series(in, ipgf, iset) = series(in, ipgf, iset)* &
                                              (1.0_dp - beta_ratio**(REAL(in, dp) + 0.5_dp))
            END DO
         END DO
      END DO

   END SUBROUTINE prepare_radial_series

! **************************************************************************************************
!> \brief Evaluate the auxiliary potentials v_P(r) for the erfc(omega*r)/r kernel.
!> For a primitive H_l(r) exp(-a r²), H_l is a homogeneous harmonic polynomial.
!> (1) Its potential is H_l(r) (-1/a)^l d^l F(s)/ds^l, s=r².
!> (2) For the full Coulomb kernel, F(r²)=pi^(3/2) erf(sqrt(a)r)/(a^(3/2) r);
!> radial_potential_short subtracts the long-range part of erfc(omega*r)/r.
!> The potential columns are the contracted functions P of contracted%coefficients
!> (AUTO_RI functions or the plain RI basis); a final GEMM applies the contraction.
!> The erfc kernel and its origin series are those prepared in contracted.
!> \param values ...
!> \param points ...
!> \param basis ...
!> \param center ...
!> \param contracted ...
! **************************************************************************************************
   SUBROUTINE contract_selected_potentials(values, points, basis, center, contracted)
      REAL(dp), CONTIGUOUS, INTENT(OUT)                  :: values(:, :)
      REAL(dp), INTENT(IN)                               :: points(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), INTENT(IN)                               :: center(3)
      TYPE(radial_interpolation_type), INTENT(IN), &
         TARGET                                          :: contracted

      INTEGER                                            :: handle, iangular, ico, ig, imesh, ipgf, &
                                                            irow, iset, isgf, ishell, l, lmax, m, &
                                                            nchannels, ng, np, nrows, nwork
      INTEGER, ALLOCATABLE                               :: angular_count(:, :, :), &
                                                            angular_index(:), &
                                                            angular_term(:, :, :, :), &
                                                            mesh_channel(:), minimum_l(:)
      LOGICAL                                            :: factorized, in_table, parallel_points, &
                                                            tabulated, use_mesh
      LOGICAL, ALLOCATABLE                               :: angular_used(:, :)
      REAL(dp)                                           :: mesh_delta, mesh_position, omega, &
                                                            radial_sum, radius, radius2, &
                                                            relative(3), spherical_value, x
      REAL(dp), ALLOCATABLE :: angular_values(:, :, :), coordinate_power(:, :), long_radial(:), &
         monomial(:), radial(:), radial_channels(:, :, :), radial_mesh(:), work(:, :)
      REAL(dp), CONTIGUOUS, POINTER                      :: coefficient_matrix(:, :)
      REAL(dp), POINTER                                  :: series(:, :, :), short_series(:, :, :)

      CALL timeset('d_lP_potential_prepare', handle)

      ng = SIZE(points, 2)
      np = SIZE(values, 2)
      lmax = MAXVAL(basis%lmax)
      nrows = SIZE(contracted%coefficients, 1)
      factorized = contracted%spherical_factorized
      nwork = nrows
      IF (factorized) nwork = basis%nsgf
      ALLOCATE (angular_used(0:lmax, basis%nset), source=.FALSE.)
      ALLOCATE (minimum_l(basis%nset))
      DO iset = 1, basis%nset
         minimum_l(iset) = MINVAL(basis%l(1:basis%nshell(iset), iset))
         DO ishell = 1, basis%nshell(iset)
            angular_used(basis%l(ishell, iset), iset) = .TRUE.
         END DO
      END DO
      omega = contracted%attenuation
      ALLOCATE (work(nwork, ng))
      series => contracted%series
      short_series => contracted%short_series
      ! Primitive/angular row labels are invariant over the integration points.
      ! Evaluate each tabulated radial channel once, then apply its Cartesian monomial.
      nchannels = 0
      ALLOCATE (mesh_channel(nrows), angular_index(nrows))
      IF (ALLOCATED(contracted%mesh_moments)) THEN
         nchannels = SIZE(contracted%mesh_moments, 2)
         irow = 0
         DO iset = 1, basis%nset
            DO ipgf = 1, basis%npgf(iset)
               DO l = 0, basis%lmax(iset)
                  IF (.NOT. angular_used(l, iset)) CYCLE
                  DO ico = ncoset(l - 1) + 1, ncoset(l)
                     irow = irow + 1
                     mesh_channel(irow) = contracted%mesh_first(ipgf, iset) + l - minimum_l(iset)
                     angular_index(irow) = ico
                  END DO
               END DO
            END DO
         END DO
      END IF
      ! Preserve the order of the nonzero polynomial terms; omit only exact zeros.
      IF (factorized) THEN
         ALLOCATE (angular_count(MAXVAL(nso(0:lmax)), 0:lmax, basis%nset), source=0)
         ALLOCATE (angular_term(SIZE(contracted%spherical_angular, 1), MAXVAL(nso(0:lmax)), 0:lmax, basis%nset))
         DO iset = 1, basis%nset
            DO l = 0, basis%lmax(iset)
               IF (.NOT. angular_used(l, iset)) CYCLE
               DO m = 1, nso(l)
                  DO ico = ncoset(l - 1) + 1, ncoset(l)
                     IF (ABS(contracted%spherical_angular(ico - ncoset(l - 1), m, l, iset)) <= 0.0_dp) CYCLE
                     angular_count(m, l, iset) = angular_count(m, l, iset) + 1
                     angular_term(angular_count(m, l, iset), m, l, iset) = ico
                  END DO
               END DO
            END DO
         END DO
      END IF
      tabulated = ALLOCATED(contracted%short_moments)

      CALL timestop(handle)

      CALL timeset('d_lP_potential_radial', handle)

      parallel_points = .TRUE.
!$    parallel_points = .NOT. omp_in_parallel()
      !$OMP PARALLEL DEFAULT(NONE) IF(parallel_points) &
      !$OMP SHARED(ng, points, center, basis, series, short_series, omega, &
      !$OMP work, indco, ncoset, lmax, angular_used, minimum_l, contracted, tabulated, &
      !$OMP nchannels, nrows, nwork, factorized, mesh_channel, angular_index, nso, angular_count, angular_term) &
      !$OMP PRIVATE(ig, ico, relative, radius, radius2, iset, ipgf, radial, monomial, l, irow, long_radial, &
      !$OMP x, in_table, coordinate_power, imesh, mesh_delta, mesh_position, use_mesh, radial_mesh, &
      !$OMP radial_channels, angular_values, isgf, m, radial_sum, spherical_value, iangular)
      ALLOCATE (radial(0:lmax), long_radial(0:lmax))
      ALLOCATE (monomial(ncoset(lmax)), coordinate_power(3, 0:lmax), radial_mesh(nchannels))
      IF (factorized) ALLOCATE (radial_channels(0:lmax, MAXVAL(basis%npgf), basis%nset), &
                                angular_values(MAXVAL(nso(0:lmax)), 0:lmax, basis%nset))
      !$OMP DO SCHEDULE(STATIC)
      DO ig = 1, ng
         relative = points(:, ig) - center
         radius2 = SUM(relative**2)
         coordinate_power(:, 0) = 1.0_dp
         DO l = 1, lmax
            coordinate_power(:, l) = coordinate_power(:, l - 1)*relative
         END DO
         DO ico = 1, SIZE(monomial)
            monomial(ico) = coordinate_power(1, indco(1, ico))* &
                            coordinate_power(2, indco(2, ico))*coordinate_power(3, indco(3, ico))
         END DO
         use_mesh = .FALSE.
         mesh_delta = 0.0_dp
         imesh = 0
         IF (ALLOCATED(contracted%mesh_moments)) THEN
            mesh_position = radius2/contracted%mesh_spacing - contracted%mesh_offset
            IF (mesh_position >= -0.5_dp .AND. mesh_position <= REAL(UBOUND(contracted%mesh_moments, 3), dp) + 0.5_dp) THEN
               imesh = NINT(mesh_position)
               IF (imesh >= 0 .AND. imesh <= UBOUND(contracted%mesh_moments, 3)) THEN
                  mesh_delta = radius2 - contracted%mesh_spacing*(REAL(imesh, dp) + contracted%mesh_offset)
                  use_mesh = ABS(mesh_delta) <= 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, radius2)
                  use_mesh = use_mesh .AND. contracted%mesh_exponent*ABS(mesh_delta) <= eps_mesh_phase
               END IF
            END IF
         END IF
         IF (use_mesh) THEN
            IF (factorized .AND. ALLOCATED(contracted%spherical_mesh_first)) THEN
               ! The radial shell contraction was formed once for each squared-radius node.
            ELSE
               radial_mesh(:) = contracted%mesh_moments(1, :, imesh) + &
                                mesh_delta*contracted%mesh_moments(2, :, imesh)
            END IF
            IF (factorized .AND. .NOT. ALLOCATED(contracted%spherical_mesh_first)) THEN
               DO iset = 1, basis%nset
                  DO ipgf = 1, basis%npgf(iset)
                     DO l = 0, basis%lmax(iset)
                        IF (.NOT. angular_used(l, iset)) CYCLE
                        radial_channels(l, ipgf, iset) = &
                           radial_mesh(contracted%mesh_first(ipgf, iset) + l - minimum_l(iset))
                     END DO
                  END DO
               END DO
            ELSE IF (.NOT. factorized) THEN
               DO irow = 1, nrows
                  work(irow, ig) = monomial(angular_index(irow))*radial_mesh(mesh_channel(irow))
               END DO
               CYCLE
            END IF
         ELSE
            radius = SQRT(radius2)
            irow = 0
            DO iset = 1, basis%nset
               DO ipgf = 1, basis%npgf(iset)
                  x = basis%zet(ipgf, iset)*radius2
                  in_table = .FALSE.
                  IF (tabulated) in_table = x <= contracted%short_abscissae(UBOUND(contracted%short_abscissae, 1))
                  IF (in_table) THEN
        CALL evaluate_short_radial_table(x, contracted%short_abscissae, contracted%short_moments(minimum_l(iset):, :, ipgf, iset), &
                                                      radial(minimum_l(iset):basis%lmax(iset)))
                  ELSE
                     CALL radial_potential_short(radius, basis%zet(ipgf, iset), omega, &
                                                 series(:, ipgf, iset), short_series(:, ipgf, iset), &
                                                 radial(0:basis%lmax(iset)), long_radial)
                  END IF
                  IF (factorized) THEN
                     DO l = 0, basis%lmax(iset)
                        IF (angular_used(l, iset)) radial_channels(l, ipgf, iset) = radial(l)
                     END DO
                  ELSE
                     DO l = 0, basis%lmax(iset)
                        IF (.NOT. angular_used(l, iset)) CYCLE
                        DO ico = ncoset(l - 1) + 1, ncoset(l)
                           irow = irow + 1
                           work(irow, ig) = monomial(ico)*radial(l)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END IF
         IF (factorized) THEN
            work(:, ig) = 0.0_dp
            DO iset = 1, basis%nset
               DO l = 0, basis%lmax(iset)
                  IF (.NOT. angular_used(l, iset)) CYCLE
                  DO m = 1, nso(l)
                     spherical_value = 0.0_dp
                     DO iangular = 1, angular_count(m, l, iset)
                        ico = angular_term(iangular, m, l, iset)
                        spherical_value = spherical_value + monomial(ico)* &
                                          contracted%spherical_angular(ico - ncoset(l - 1), m, l, iset)
                     END DO
                     angular_values(m, l, iset) = spherical_value
                  END DO
               END DO
               DO ishell = 1, basis%nshell(iset)
                  l = basis%l(ishell, iset)
                  IF (use_mesh .AND. ALLOCATED(contracted%spherical_mesh_first)) THEN
                     irow = contracted%spherical_mesh_first(ishell, iset)
                     radial_sum = contracted%mesh_moments(1, irow, imesh) + &
                                  mesh_delta*contracted%mesh_moments(2, irow, imesh)
                  ELSE
                     radial_sum = SUM(radial_channels(l, 1:basis%npgf(iset), iset)* &
                                      contracted%spherical_primitive(1:basis%npgf(iset), ishell, iset))
                  END IF
                  DO m = 1, nso(l)
                     isgf = basis%first_sgf(ishell, iset) + m - 1
                     work(isgf, ig) = radial_sum*angular_values(m, l, iset)
                  END DO
               END DO
            END DO
         END IF
      END DO
      !$OMP END DO
      DEALLOCATE (radial, long_radial)
      DEALLOCATE (monomial, coordinate_power, radial_mesh)
      IF (factorized) DEALLOCATE (radial_channels, angular_values)
      !$OMP END PARALLEL

      CALL timestop(handle)

      CALL timeset('d_lP_potential_RI_transform', handle)

      ! The caller supplies the full leading dimension; write only the active point rows.
      IF (factorized) THEN
         coefficient_matrix => contracted%spherical_selection
      ELSE
         coefficient_matrix => contracted%coefficients
      END IF
      CALL dgemm('T', 'N', ng, np, nwork, 1.0_dp, work, nwork, &
                 coefficient_matrix, nwork, 0.0_dp, values, SIZE(values, 1))
      IF (factorized) DEALLOCATE (angular_count, angular_term)
      DEALLOCATE (angular_used, minimum_l, mesh_channel, angular_index)

      CALL timestop(handle)

   END SUBROUTINE contract_selected_potentials

! **************************************************************************************************
!> \brief Gaussian potential for erfc(omega*r)/r via two full-Coulomb potentials.
!> With beta=a*omega²/(a+omega²), the long-range part equals the beta potential
!> multiplied in angular channel l by (beta/a)^(l+3/2). The origin coefficients
!> combine both terms before evaluation to avoid subtraction at r=0.
!> \param radius ...
!> \param exponent ...
!> \param omega ...
!> \param coefficients ...
!> \param short_coefficients ...
!> \param radial ...
!> \param long_radial ...
! **************************************************************************************************
   SUBROUTINE radial_potential_short(radius, exponent, omega, coefficients, short_coefficients, radial, long_radial)
      REAL(dp), INTENT(IN)                               :: radius, exponent, omega, &
                                                            coefficients(0:), &
                                                            short_coefficients(0:)
      REAL(dp), INTENT(OUT)                              :: radial(0:)
      REAL(dp), INTENT(INOUT)                            :: long_radial(0:)

      INTEGER                                            :: il
      REAL(dp)                                           :: beta, beta_ratio, factor, lower_exp, &
                                                            ratio_power, root_x, upper_exp, x

      IF (exponent*radius**2 < 0.5_dp) THEN
         CALL radial_potential(radius, exponent, short_coefficients, radial)
         RETURN
      END IF
      beta_ratio = omega**2/(exponent + omega**2)
      beta = exponent*beta_ratio
      x = exponent*radius**2
      IF (x >= REAL(2*UBOUND(radial, 1) + 8, dp)) THEN
         ! H_l(x,b)=integral_{sqrt(b)}^1 t^(2l) exp(-x t^2) dt; radial_l=(2 pi/a) H_l.
         ! Its upward recurrence is stable here and avoids subtracting two Coulomb derivatives.
         root_x = SQRT(x)
         factor = 2.0_dp*pi/exponent
         radial(0) = factor*SQRT(pi)/(2.0_dp*root_x)* &
                     (ERFC(SQRT(beta_ratio)*root_x) - ERFC(root_x))
         upper_exp = EXP(-x)
         lower_exp = EXP(-beta_ratio*x)
         ratio_power = SQRT(beta_ratio)
         DO il = 1, UBOUND(radial, 1)
            radial(il) = (REAL(2*il - 1, dp)*radial(il - 1) + &
                          factor*(ratio_power*lower_exp - upper_exp))/(2.0_dp*x)
            ratio_power = ratio_power*beta_ratio
         END DO
         RETURN
      END IF
      CALL radial_potential(radius, exponent, coefficients, radial)
      CALL radial_potential(radius, beta, coefficients, long_radial(0:UBOUND(radial, 1)))
      DO il = 0, UBOUND(radial, 1)
         radial(il) = radial(il) - beta_ratio**(REAL(il, dp) + 1.5_dp)*long_radial(il)
      END DO

   END SUBROUTINE radial_potential_short

! **************************************************************************************************
!> \brief Expand the full-Coulomb F of Eq. (2) in x=a r²: F=(2 pi/a) sum_n c_n x^n near the origin.
!> The coefficients c_n are independent of the exponent a.
!> \param coefficients ...
! **************************************************************************************************
   SUBROUTINE origin_series(coefficients)
      REAL(dp), INTENT(OUT)                              :: coefficients(0:)

      INTEGER                                            :: ik, in
      REAL(dp) :: factorial, weight, terms(0:UBOUND(coefficients, 1))

      weight = 1.0_dp
      DO ik = 0, UBOUND(terms, 1)
         terms(ik) = weight
         weight = weight*2.0_dp/REAL(2*ik + 3, dp)
      END DO
      DO in = 0, UBOUND(coefficients, 1)
         coefficients(in) = 0.0_dp
         factorial = 1.0_dp
         DO ik = in, 0, -1
            coefficients(in) = coefficients(in) + factorial*terms(ik)
            factorial = -factorial/REAL(in - ik + 1, dp)
         END DO
      END DO

   END SUBROUTINE origin_series

! **************************************************************************************************
!> \brief Evaluate the radial derivatives in Eq. (1) by series or Taylor differentiation of (2).
!> The full Coulomb kernel is used; radial_potential_short forms erfc(omega*r)/r from it.
!> \param radius ...
!> \param exponent ...
!> \param coefficients ...
!> \param radial ...
! **************************************************************************************************
   SUBROUTINE radial_potential(radius, exponent, coefficients, radial)
      REAL(dp), INTENT(IN)                               :: radius, exponent, coefficients(0:)
      REAL(dp), INTENT(OUT)                              :: radial(0:)

      INTEGER                                            :: ij, il, in, lmax
      REAL(dp) :: factor, power, prefactor, value, x, jet(0:UBOUND(radial, 1)), &
         gaussian(0:UBOUND(radial, 1)), derivative(0:UBOUND(radial, 1)), &
         bracket(0:UBOUND(radial, 1))

      lmax = UBOUND(radial, 1)
      x = exponent*radius**2
      IF (x < 0.5_dp) THEN
         DO il = 0, lmax
            value = 0.0_dp
            power = 1.0_dp
            DO in = il, UBOUND(coefficients, 1)
               factor = 1.0_dp
               DO ij = in - il + 1, in
                  factor = factor*REAL(ij, dp)
               END DO
               value = value + coefficients(in)*factor*power
               power = power*x
            END DO
            radial(il) = 2.0_dp*pi/exponent*(-1.0_dp)**il*value
         END DO
         RETURN
      END IF
      bracket(0:lmax) = 0.0_dp
      gaussian(0) = EXP(-exponent*radius**2)
      IF (lmax > 0) gaussian(1) = -2.0_dp*exponent*radius*gaussian(0)
      DO in = 2, lmax
         gaussian(in) = -2.0_dp*exponent* &
                        (radius*gaussian(in - 1) + gaussian(in - 2))/REAL(in, dp)
      END DO
      DO in = 1, lmax
         bracket(in) = bracket(in) + &
                       4.0_dp*SQRT(exponent/pi)*gaussian(in - 1)/REAL(in, dp)
      END DO
      bracket(0) = 2.0_dp*ERF(SQRT(exponent)*radius)
      prefactor = pi**1.5_dp/(2.0_dp*exponent**1.5_dp)
      jet(0) = bracket(0)/radius
      DO in = 1, lmax
         jet(in) = (bracket(in) - jet(in - 1))/radius
      END DO
      radial(0) = prefactor*jet(0)
      DO il = 1, lmax
         derivative(0:lmax) = 0.0_dp
         derivative(0) = jet(1)/(2.0_dp*radius)
         DO in = 1, lmax - il
            derivative(in) = (REAL(in + 1, dp)*jet(in + 1)/2.0_dp - derivative(in - 1))/radius
         END DO
         jet(0:lmax) = derivative(0:lmax)
         radial(il) = prefactor*(-1.0_dp/exponent)**il*jet(0)
      END DO

   END SUBROUTINE radial_potential

! **************************************************************************************************
!> \brief Place radial nodes at x=0.02*i up to x=2, then at x=2*1.01**(i-100).
!> The final node retains the original table extent. For a nearest-node displacement δ,
!> the sixth-order Taylor remainder of H_j is bounded by
!> |δ|**7/7! * integral_0^1 t**(2*j+14)*exp(-(x-|δ|)*t**2) dt.
!> On the geometric part |δ| <= 0.005*x, so this bound is below 1.1e-17 for j>=0.
!> This is an absolute bound before the Gaussian prefactor 2*pi/a; roundoff is separate.
!> \param xmax ...
!> \param abscissae ...
! **************************************************************************************************
   SUBROUTINE prepare_short_radial_grid(xmax, abscissae)
      REAL(dp), INTENT(IN)                               :: xmax
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: abscissae(:)

      INTEGER                                            :: ix, last

      CPASSERT(xmax > 0.0_dp)
      IF (xmax <= short_grid_switch) THEN
         last = CEILING(xmax/short_grid_spacing)
      ELSE
         last = n_short_uniform + CEILING(LOG(xmax/short_grid_switch)/LOG(short_grid_ratio))
      END IF
      ALLOCATE (abscissae(0:last))
      DO ix = 0, last - 1
         IF (ix <= n_short_uniform) THEN
            abscissae(ix) = short_grid_spacing*REAL(ix, dp)
         ELSE
            abscissae(ix) = short_grid_switch*EXP(REAL(ix - n_short_uniform, dp)*LOG(short_grid_ratio))
         END IF
      END DO
      abscissae(last) = xmax

   END SUBROUTINE prepare_short_radial_grid

! **************************************************************************************************
!> \brief Tabulate H_j(x,b)=integral_sqrt(b)^1 t^(2j)*exp(-x*t²)dt for erfc potentials.
!> A 64-point positive quadrature supplies near-field moments; a stable upward recurrence
!> supplies far-field moments for b<0.9. Narrow intervals retain quadrature to avoid cancellation.
!> Uniform mesh tables select spacing and offset; nonuniform lookup tables supply their abscissae.
!> The caller bounds storage, including one primitive preparation buffer.
!> \param a ...
!> \param omega ...
!> \param lmax ...
!> \param table ...
!> \param ninterval ...
!> \param point_spacing ...
!> \param point_offset ...
!> \param last_moment ...
!> \param abscissae ...
! **************************************************************************************************
   SUBROUTINE prepare_short_radial_table(a, omega, lmax, table, ninterval, point_spacing, point_offset, &
                                         last_moment, abscissae)
      REAL(dp), INTENT(IN)                               :: a, omega
      INTEGER, INTENT(IN)                                :: lmax
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: table(:, :)
      INTEGER, INTENT(IN)                                :: ninterval
      REAL(dp), INTENT(IN), OPTIONAL                     :: point_spacing, point_offset
      INTEGER, INTENT(IN), OPTIONAL                      :: last_moment
      REAL(dp), INTENT(IN), OPTIONAL                     :: abscissae(0:)

      INTEGER, PARAMETER                                 :: n_legendre = 64

      INTEGER                                            :: ig, il, ix, last, order, quadrature_last
      REAL(dp) :: b, factor, legendre_nodes(n_legendre), legendre_weights(n_legendre), lower, &
         lower_exp, offset, ratio_power, root_x, spacing, t, t2, term, upper_exp, width, x

      CPASSERT(PRESENT(point_spacing) .OR. PRESENT(abscissae))
      last = ninterval
      spacing = 0.0_dp
      offset = 0.0_dp
      order = lmax + 6
      IF (PRESENT(point_spacing)) spacing = point_spacing
      IF (PRESENT(point_offset)) offset = point_offset
      IF (PRESENT(last_moment)) order = last_moment
      ALLOCATE (table(0:order, 0:last), source=0.0_dp)
      b = omega**2/(a + omega**2)
      lower = SQRT(b)
      width = (a/(a + omega**2))/(1.0_dp + lower)
      quadrature_last = last
      IF (b < recurrence_b_max) THEN
         IF (PRESENT(abscissae)) THEN
            quadrature_last = COUNT(abscissae <= REAL(2*order + 8, dp)) - 1
         ELSE
            quadrature_last = INT(MIN(REAL(last, dp), MAX(0.0_dp, REAL(2*order + 8, dp)/spacing - offset)))
         END IF
      END IF
      CALL gauss_legendre(n_legendre, legendre_nodes, legendre_weights)
      DO ig = 1, n_legendre
         t = lower + 0.5_dp*width*(1.0_dp + legendre_nodes(ig))
         t2 = t*t
         DO ix = 0, quadrature_last
            x = spacing*(REAL(ix, dp) + offset)
            IF (PRESENT(abscissae)) x = abscissae(ix)
            term = 0.5_dp*width*legendre_weights(ig)*EXP(-x*t2)*2.0_dp*pi/a
            DO il = 0, UBOUND(table, 1)
               table(il, ix) = table(il, ix) + term
               term = term*t2
            END DO
         END DO
      END DO
      factor = 2.0_dp*pi/a
      DO ix = quadrature_last + 1, last
         x = spacing*(REAL(ix, dp) + offset)
         IF (PRESENT(abscissae)) x = abscissae(ix)
         root_x = SQRT(x)
         table(0, ix) = factor*SQRT(pi)/(2.0_dp*root_x)*(ERFC(lower*root_x) - ERFC(root_x))
         lower_exp = EXP(-b*x)
         upper_exp = EXP(-x)
         ratio_power = lower
         DO il = 1, UBOUND(table, 1)
            table(il, ix) = (REAL(2*il - 1, dp)*table(il - 1, ix) + &
                             factor*(ratio_power*lower_exp - upper_exp))/(2.0_dp*x)
            ratio_power = ratio_power*b
         END DO
      END DO

   END SUBROUTINE prepare_short_radial_table

! **************************************************************************************************
!> \brief Evaluate H_l at the nearest table point with its first six derivatives.
!> Since d^k H_l/dx^k=(-1)^k H_(l+k), the stored higher moments provide the derivatives.
!> The nonuniform-node error bound is given in prepare_short_radial_grid.
!> \param x ...
!> \param abscissae ...
!> \param table ...
!> \param radial ...
! **************************************************************************************************
   SUBROUTINE evaluate_short_radial_table(x, abscissae, table, radial)
      REAL(dp), INTENT(IN)                               :: x, abscissae(0:), table(0:, 0:)
      REAL(dp), INTENT(OUT)                              :: radial(0:)

      REAL(dp), PARAMETER :: &
         inverse(6) = [1.0_dp, 0.5_dp, 1.0_dp/3.0_dp, 0.25_dp, 0.2_dp, 1.0_dp/6.0_dp]

      INTEGER                                            :: ik, il, ileft, imiddle, iright, ix
      REAL(dp)                                           :: delta, value

      ileft = 0
      iright = UBOUND(abscissae, 1)
      DO WHILE (iright - ileft > 1)
         imiddle = (ileft + iright)/2
         IF (x < abscissae(imiddle)) THEN
            iright = imiddle
         ELSE
            ileft = imiddle
         END IF
      END DO
      ix = ileft
      IF (abscissae(iright) - x < x - abscissae(ileft)) ix = iright
      delta = abscissae(ix) - x
      DO il = 0, UBOUND(radial, 1)
         value = table(il + 6, ix)
         DO ik = 5, 0, -1
            value = table(il + ik, ix) + delta*value*inverse(ik + 1)
         END DO
         radial(il) = value
      END DO

   END SUBROUTINE evaluate_short_radial_table

! **************************************************************************************************
!> \brief Tabulate H_l(r²) only when squared mesh distances share a one-dimensional lattice.
!> With M=SᵀS and d=r_origin-r_atom, r²=nᵀMn+2nᵀSᵀd+d². Integer coefficients
!> in units q imply r²=q(i+offset). Whole-cell translations preserve this mesh.
!> H_l(s+δ)=H_l(s)-a δ H_(l+1)(s)+O((a δ)²); each lookup checks the residual.
!> Tables cover up to 40 bohr within the storage budget; other points use the generic lookup.
!> Account for all stored channels, their index map and the largest temporary moment table.
!> \param basis ...
!> \param contracted ...
!> \param step ...
!> \param displacement ...
!> \param table_bytes ...
! **************************************************************************************************
   SUBROUTINE prepare_mesh_radial_moments(basis, contracted, step, displacement, table_bytes)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(radial_interpolation_type), INTENT(INOUT)     :: contracted
      REAL(dp), INTENT(IN)                               :: step(3, 3), displacement(3)
      REAL(dp), INTENT(INOUT)                            :: table_bytes

      INTEGER :: idir, ifirst, il, ipgf, iset, ishell, jdir, last, max_l, max_primitive, &
         minimum_l, nchannels, nshell_total, nstored, shell_index
      LOGICAL                                            :: direct_shells
      REAL(dp)                                           :: bytes_per_node, index_bytes, linear(3), &
                                                            metric(3, 3), offset, ratio, &
                                                            required_bytes, spacing
      REAL(dp), ALLOCATABLE                              :: moments(:, :)

      IF (basis%nset < 1) RETURN
      IF (ANY(basis%nshell(1:basis%nset) < 1)) RETURN
      IF (ANY(basis%npgf(1:basis%nset) < 1)) RETURN
      IF (contracted%attenuation <= 0.0_dp .OR. table_bytes <= 0.0_dp) RETURN
      metric = MATMUL(TRANSPOSE(step), step)
      spacing = MIN(metric(1, 1), metric(2, 2), metric(3, 3))
      IF (spacing <= TINY(1.0_dp)) RETURN
      DO idir = 1, 3
         DO jdir = idir, 3
            ratio = metric(idir, jdir)/spacing
            IF (idir /= jdir) ratio = 2.0_dp*ratio
            IF (ABS(ratio - ANINT(ratio)) > eps_integer) RETURN
         END DO
      END DO
      linear = 2.0_dp*MATMUL(TRANSPOSE(step), displacement)/spacing
      IF (ANY(ABS(linear - ANINT(linear)) > eps_integer)) RETURN
      offset = MODULO(SUM(displacement**2)/spacing, 1.0_dp)
      IF (mesh_table_radius**2/spacing > REAL(HUGE(last) - 1, dp)) RETURN
      last = CEILING(mesh_table_radius**2/spacing)
      IF (last < 1) RETURN
      nchannels = 0
      DO iset = 1, basis%nset
         minimum_l = MINVAL(basis%l(1:basis%nshell(iset), iset))
         nchannels = nchannels + basis%npgf(iset)*(basis%lmax(iset) - minimum_l + 1)
      END DO
      nshell_total = SUM(basis%nshell(1:basis%nset))
      max_l = MAXVAL(basis%lmax(1:basis%nset))
      max_primitive = MAXVAL(basis%npgf(1:basis%nset))
      direct_shells = contracted%spherical_factorized .AND. nshell_total <= nchannels
      nstored = nchannels
      IF (direct_shells) nstored = nshell_total
      bytes_per_node = real_bytes*(2.0_dp*REAL(nstored, dp) + REAL(max_l + 2, dp))
      index_bytes = REAL(int_size, dp)*REAL(max_primitive, dp)*REAL(basis%nset, dp)
      IF (direct_shells) THEN
         index_bytes = index_bytes + &
                       REAL(int_size, dp)*REAL(MAXVAL(basis%nshell)*basis%nset, dp)
      END IF
      ! A smaller table still accelerates nearby points; its upper bound selects the generic fallback.
      last = MIN(last, INT(MIN(REAL(HUGE(last), dp), &
                               MAX(0.0_dp, (table_bytes - index_bytes)/bytes_per_node))) - 1)
      IF (last < 1) RETURN
      required_bytes = bytes_per_node*REAL(last + 1, dp) + index_bytes
      IF (required_bytes > table_bytes) RETURN
      ALLOCATE (contracted%mesh_moments(2, nstored, 0:last))
      ALLOCATE (contracted%mesh_first(max_primitive, basis%nset), source=0)
      contracted%mesh_spacing = spacing
      contracted%mesh_offset = offset
      IF (direct_shells) THEN
         ALLOCATE (contracted%spherical_mesh_first(MAXVAL(basis%nshell), basis%nset), source=0)
         shell_index = 0
         DO iset = 1, basis%nset
            DO ishell = 1, basis%nshell(iset)
               shell_index = shell_index + 1
               contracted%spherical_mesh_first(ishell, iset) = shell_index
            END DO
         END DO
         contracted%mesh_moments = 0.0_dp
      END IF
      ifirst = 1
      DO iset = 1, basis%nset
         minimum_l = MINVAL(basis%l(1:basis%nshell(iset), iset))
         DO ipgf = 1, basis%npgf(iset)
            contracted%mesh_exponent = MAX(contracted%mesh_exponent, basis%zet(ipgf, iset))
            contracted%mesh_first(ipgf, iset) = ifirst
            CALL prepare_short_radial_table(basis%zet(ipgf, iset), contracted%attenuation, &
                                            basis%lmax(iset), moments, last, &
                                            basis%zet(ipgf, iset)*spacing, offset, basis%lmax(iset) + 1)
            IF (direct_shells) THEN
               ! Accumulate primitive contributions in their original order for each shell.
               ! The larger primitive-channel table is never materialized.
               DO ishell = 1, basis%nshell(iset)
                  il = basis%l(ishell, iset)
                  shell_index = contracted%spherical_mesh_first(ishell, iset)
                  contracted%mesh_moments(1, shell_index, :) = &
                     contracted%mesh_moments(1, shell_index, :) + &
                     contracted%spherical_primitive(ipgf, ishell, iset)*moments(il, :)
                  contracted%mesh_moments(2, shell_index, :) = &
                     contracted%mesh_moments(2, shell_index, :) + &
                     contracted%spherical_primitive(ipgf, ishell, iset)*(-basis%zet(ipgf, iset)*moments(il + 1, :))
               END DO
               ifirst = ifirst + basis%lmax(iset) - minimum_l + 1
            ELSE
               DO il = minimum_l, basis%lmax(iset)
                  contracted%mesh_moments(1, ifirst, :) = moments(il, :)
                  contracted%mesh_moments(2, ifirst, :) = -basis%zet(ipgf, iset)*moments(il + 1, :)
                  ifirst = ifirst + 1
               END DO
            END IF
            DEALLOCATE (moments)
         END DO
      END DO
      table_bytes = table_bytes - required_bytes

   END SUBROUTINE prepare_mesh_radial_moments

END MODULE ri_rs_kp_integrals

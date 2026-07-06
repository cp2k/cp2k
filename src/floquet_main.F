!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Floquet stuff
!> \par History
!> \author Shridhar Shanbhag (27.01.2026)
! **************************************************************************************************
MODULE floquet_main
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE floquet_utils,                   ONLY: build_floquet_matrix_local,&
                                              calculate_floquet_spectral_function,&
                                              check_floquet_conversion,&
                                              floquet_determine_subgroup_size,&
                                              floquet_precompute_kp_data,&
                                              get_quasi_energies,&
                                              write_floquet_spectral_function,&
                                              write_quasi_energy
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_any_source,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_request_type,&
                                              mp_testall,&
                                              mp_testany,&
                                              mp_waitall
   USE physcon,                         ONLY: a_bohr,&
                                              evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'floquet_main'

   ! Public subroutines
   PUBLIC :: floquet

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE floquet(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'floquet'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: dipole_kp_spin
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: dipole
      INTEGER :: group_size, handle, ikp, ispin, loc_idx, local_src, max_f_index, max_printed_ikp, &
                 my_group, n_E, n_f_size, n_fbands, n_kp_dos, n_recv, n_send, n_spin, n_streams, nao, &
                 ngroups, nkp, nkp_start, owner, stride_kp, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: group_distribution, remaining, req_ikp
      LOGICAL                                            :: all_sent, am_collector
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: received
      REAL(KIND=dp)                                      :: energy_step, energy_window
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a_k, e_k_kp_spin, eigenvalues, &
                                                            quasi_energies
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: de_dk_kp_spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: all_a_k, all_quasi_energies, e_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: de_dk
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors, floquet_matrix
      TYPE(cp_fm_struct_type), POINTER                   :: floquet_struct
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_comm_type)                                 :: comm_stream
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_reqs, send_reqs

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      unit_nr = bs_env%unit_nr
      nkp = bs_env%nkp_bs_and_DOS
      nkp_start = bs_env%nkp_only_DOS
      n_spin = bs_env%n_spin
      max_f_index = bs_env%max_floquet_index
      n_fbands = 1 + 2*max_f_index
      n_f_size = nao*n_fbands

      energy_step = bs_env%energy_step_floquet
      energy_window = 2*bs_env%energy_window_floquet

      n_E = INT(energy_window/energy_step)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(T2,A)') ' '
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A,A51,A27)') '-', 'FLOQUET BANDSTRUCTURE CALCULATION', '-'
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A)') ' '

         WRITE (unit_nr, '(T2,A,T56,ES12.2E2)') &
            "FLOQUET PARAMETERS      Amplitude [V/m]:", bs_env%floquet_amplitude*evolt/a_bohr
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Frequency [eV]:", bs_env%floquet_omega*evolt
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T56,3F8.4)') &
            "                        Polarisation:", bs_env%floquet_polarisation(1:3)
         WRITE (unit_nr, '(T2,A,T56,3F8.4)') &
            "                        Phase offsets:", pi*bs_env%floquet_phi(1:3)
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T56,I12)') &
            "                        Max Floquet index:", bs_env%max_floquet_index
         WRITE (unit_nr, '(T2,A,T56,I12)') &
            "                        Floquet Hamiltonian Size:", n_f_size
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Energy window [eV]:", bs_env%energy_window_floquet*evolt
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Energy step [eV]:", bs_env%energy_step_floquet*evolt
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Broadening [eV]:", bs_env%broadening_floquet*evolt
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(A)') ""

         WRITE (unit_nr, '(T2,A)') &
            "We construct the Floquet-Bloch Hamiltonian and diagonalise it to obtain the"
         WRITE (unit_nr, '(T2,A)') &
            "eigenvalues which give us the quasi-energies stored in QUASI_ENERGIES.bs"
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(T2,A)') &
            "The k-resolved density of states is obtained by computing the trace of"
         WRITE (unit_nr, '(T2,A)') &
            "the retarded Green's function and stored in FLOQUET_DOS.out"
         WRITE (unit_nr, '(T2,A)') &
            "DOS(ω,k) = -1/π*Im[Tr_KS(G^R(ω,k))]"
         WRITE (unit_nr, '(A)') ""
      END IF

      n_kp_dos = nkp - nkp_start

      ! Phase 1 (subgroup sizing + split): choose the number of ranks per subgroup as the
      ! larger of two floors. The memory floor requires the subgroup to hold the Floquet
      ! working set (n_fm_work_copies copies of the n_f_size^2 complex Hamiltonian) distributed
      ! across its ranks; for large systems (e.g. Bi2Te3 / TZVP) the matrix does not fit on one
      ! rank, forcing a larger subgroup. The k-point floor (floor(num_pe / n_kp_dos)) makes the
      ! subgroups as large as possible while still keeping at least n_kp_dos of them, so spare
      ! ranks speed up each diagonalisation rather than sitting idle.
      group_size = floquet_determine_subgroup_size(qs_env, bs_env)

      ! Split the global communicator into subgroups, mirroring qs_tddfpt2_subgroups.
      ! group_distribution(r) gives the subgroup of global rank r. The band data of k-point ikp
      ! lives only on global rank MOD(ikp-1,num_pe), i.e. on ranks 0..n_kp_dos-1 when ranks
      ! outnumber k-points. A contiguous split would pack those owners into the first few
      ! subgroups (the rest idle), so when num_pe > n_kp_dos we stride-interleave the split
      ! (stride = number of subgroups) to scatter the owners across all subgroups; otherwise a
      ! contiguous split (stride 1) keeps better intra-subgroup communication locality.
      IF (n_kp_dos < para_env%num_pe) THEN
         stride_kp = para_env%num_pe/group_size
      ELSE
         stride_kp = 1
      END IF
      ALLOCATE (group_distribution(0:para_env%num_pe - 1))
      ALLOCATE (para_env_sub)
      CALL para_env_sub%from_split(comm=para_env, ngroups=ngroups, &
                                   group_distribution=group_distribution, &
                                   subgroup_min_size=group_size, stride=stride_kp)
      my_group = group_distribution(para_env%mepos)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Subgroup rank stride:", stride_kp
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Number of subgroups:", ngroups
         WRITE (unit_nr, '(T2,A)') ""
         WRITE (unit_nr, '(T2,A)') "Calculations completed for:"
         WRITE (unit_nr, '(2X,A,T12,A,T22,A)') "KPoint", "Spin", "Coordinate"
      END IF

      NULLIFY (blacs_env_sub)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      ! The struct (and its BLACS context) is k-point independent and kept for the whole run;
      ! the heavy distributed cfm matrices are (de)allocated per k-point (see below).
      CALL cp_fm_struct_create(floquet_struct, &
                               nrow_global=n_f_size, &
                               ncol_global=n_f_size, &
                               context=blacs_env_sub, &
                               para_env=para_env_sub)

      ALLOCATE (eigenvalues(n_f_size))
      ALLOCATE (quasi_energies(nao))
      ALLOCATE (a_k(n_E))

      ! Per-single-(k-point,spin) broadcast buffers: the owner rank fills its slice and
      ! broadcasts them within the owning subgroup so every subgroup rank can build its piece
      ! of the Floquet matrix. These same buffers are passed to the Floquet matrix builder.
      ALLOCATE (e_k_kp_spin(nao))
      ALLOCATE (de_dk_kp_spin(3, nao))
      ALLOCATE (dipole_kp_spin(3, nao, nao))

      ! Per-(spin,k-point) result slots. Each is filled once by the owning subgroup's source
      ! rank and never overwritten, so a slot doubles as a persistent isend/irecv buffer.
      ALLOCATE (all_quasi_energies(nao, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_a_k(n_E, n_spin, n_kp_dos), source=0.0_dp)

      ! Streaming output setup (barrier-free): the collector (global I/O rank) writes k-points
      ! as they finish, in strict k-order, without any global collective. Each remote k-point's
      ! subgroup source isends its result to the collector; the collector polls with a
      ! non-blocking mp_testany between its own diagonalisations and flushes the contiguous
      ! completed prefix. n_streams counts the enabled output arrays (quasi and/or DOS); a
      ! k-point is "received" once all its streams have arrived. Tags: quasi = ikp, DOS =
      ! n_kp_dos + ikp.
      n_streams = 0
      n_streams = n_streams + 1
      n_streams = n_streams + 1
      am_collector = para_env%is_source()
      ! Dedicated communicator (duplicate of para_env, same rank order) for all streaming
      ! point-to-point, so its tags can never collide with cp_fm/cp_cfm or any other traffic on
      ! para_env.
      CALL comm_stream%from_dup(para_env)
      max_printed_ikp = 0
      n_recv = 0
      n_send = 0
      ALLOCATE (received(n_kp_dos), source=.FALSE.)
      ALLOCATE (remaining(n_kp_dos), source=n_streams)
      IF (am_collector) THEN
         ! Count the k-points computed elsewhere (one message per stream); the collector's own
         ! subgroup produces the rest locally (remaining = 0, no message expected).
         DO ikp = 1, n_kp_dos
            IF (group_distribution(MOD(ikp - 1, para_env%num_pe)) == my_group) THEN
               remaining(ikp) = 0
            ELSE
               n_recv = n_recv + n_streams
            END IF
         END DO
         ALLOCATE (recv_reqs(MAX(n_recv, 1)), req_ikp(MAX(n_recv, 1)))
         ! Pre-post one irecv per remote k-point per stream, directly into its result slot.
         n_recv = 0
         DO ikp = 1, n_kp_dos
            IF (group_distribution(MOD(ikp - 1, para_env%num_pe)) == my_group) CYCLE
            n_recv = n_recv + 1
            CALL comm_stream%irecv(all_quasi_energies(:, :, ikp), mp_any_source, recv_reqs(n_recv), tag=ikp)
            req_ikp(n_recv) = ikp
            n_recv = n_recv + 1
            CALL comm_stream%irecv(all_a_k(:, :, ikp), mp_any_source, recv_reqs(n_recv), tag=n_kp_dos + ikp)
            req_ikp(n_recv) = ikp
         END DO
      ELSE
         ! Non-collector subgroup sources send their finished k-points; size to the worst case.
         ALLOCATE (send_reqs(MAX(n_kp_dos*n_streams, 1)))
      END IF

      ! Phase 2a (global): compute the band data for ALL DOS k-points once, distributed across
      ! ranks (k-point ikp owned by global rank MOD(ikp-1,num_pe), stored at CEILING(ikp/num_pe)).
      ! Called by all ranks; the large real-space arrays are built and freed inside, so they do
      ! not coexist with the distributed Floquet matrices allocated below.
      CALL floquet_precompute_kp_data(qs_env, bs_env, e_k, de_dk, dipole)

      ! Phase 2b (subgroup loop): each DOS k-point is processed by the subgroup that contains
      ! the rank owning its data; that owner broadcasts the k-point's band data within the
      ! subgroup, then the whole subgroup builds and diagonalises H_F(k). The collector drains
      ! its mailbox and writes finished k-points at the top of every iteration (local, no barrier).
      DO ikp = 1, n_kp_dos
         ! Every rank touches MPI briefly each iteration so rendezvous transfers advance during
         ! the fast non-owned (CYCLE) iterations, not only at the terminal waitall: the collector
         ! reaps + writes, the senders progress their outstanding isends.
         IF (am_collector) THEN
            CALL stream_poll_and_write()
         ELSE IF (n_send > 0) THEN
            all_sent = mp_testall(send_reqs(1:n_send))
         END IF

         owner = MOD(ikp - 1, para_env%num_pe)
         ! Only the subgroup that owns this k-point's data works on it. The test is identical
         ! for all ranks of a subgroup, so they enter the subgroup-collective calls together.
         IF (group_distribution(owner) /= my_group) CYCLE

         ! Identify the owner's rank within the subgroup; the per-spin band data is broadcast
         ! from it to all subgroup ranks inside the spin loop below.
         local_src = -1
         IF (para_env%mepos == owner) local_src = para_env_sub%mepos
         CALL para_env_sub%max(local_src)

         CALL cp_cfm_create(floquet_matrix, floquet_struct, set_zero=.TRUE.)
         CALL cp_cfm_create(cfm_eigenvectors, floquet_struct, set_zero=.TRUE.)

         DO ispin = 1, n_spin
            ! The owner copies its single-spin slice of this k-point's distributed band data
            ! into the broadcast buffers, then broadcasts them to the whole subgroup.
            IF (para_env%mepos == owner) THEN
               loc_idx = CEILING(REAL(ikp)/para_env%num_pe)
               e_k_kp_spin(:) = e_k(ispin, loc_idx, :)
               de_dk_kp_spin(:, :) = de_dk(ispin, loc_idx, :, :)
               dipole_kp_spin(:, :, :) = dipole(ispin, loc_idx, :, :, :)
            END IF
            CALL para_env_sub%bcast(e_k_kp_spin, local_src)
            CALL para_env_sub%bcast(de_dk_kp_spin, local_src)
            CALL para_env_sub%bcast(dipole_kp_spin, local_src)

            ! Build the Floquet-Bloch Hamiltonian H_F(k) on this subgroup's matrix from the
            ! broadcast band data of this k-point - no global communication involved.
            CALL build_floquet_matrix_local(bs_env, e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, &
                                            floquet_matrix)

            ! Diagonalize the Floquet-Bloch matrix H_F(k) F_{α,k} = ε_α F_{α,k} to obtain the
            ! eigenvalues ε_α and the Floquet-Bloch coefficients F^{nu}_{α,k}. cp_cfm_heevd
            ! overwrites its input, but floquet_matrix is not needed afterwards (both the
            ! quasi-energies and the DOS are built from the eigenpairs), so diagonalize in place.
            CALL cp_cfm_set_all(cfm_eigenvectors, z_zero)
            CALL cp_cfm_heevd(floquet_matrix, cfm_eigenvectors, eigenvalues)

            ! Extract the quasi-energies (the nao eigenvalues folded into the first Floquet
            ! Brillouin zone) from the full eigenvalue spectrum.
            CALL get_quasi_energies(eigenvalues, quasi_energies)

            ! Run a conversion check to ensure that MAX_FLOQUET_INDEX,
            ! the size of the truncated matrix H_F(k), is sufficiently large.
            CALL check_floquet_conversion(qs_env, bs_env, cfm_eigenvectors)

            ! Compute the physical spectral function A_phys(k,ω) = -(1/π) Im[ Tr_KS G^R_{00}(k,ω) ]
            ! directly from the eigenpairs (spectral representation of the central sector).
            CALL calculate_floquet_spectral_function(qs_env, bs_env, eigenvalues, cfm_eigenvectors, a_k)

            ! Store this (k-point, spin) result on the subgroup source; it is sent to the
            ! collector once both spins of this k-point are done (see after the spin loop).
            IF (para_env_sub%is_source()) THEN
               all_quasi_energies(:, ispin, ikp) = quasi_energies(:)
               all_a_k(:, ispin, ikp) = a_k(:)
            END IF
         END DO

         ! This k-point is complete (all spins). The subgroup source ships it to the collector;
         ! the collector's own k-points need no message and are just marked received.
         IF (para_env_sub%is_source()) THEN
            IF (am_collector) THEN
               received(ikp) = .TRUE.
            ELSE
               n_send = n_send + 1
               CALL comm_stream%isend(all_quasi_energies(:, :, ikp), 0, send_reqs(n_send), tag=ikp)
               n_send = n_send + 1
               CALL comm_stream%isend(all_a_k(:, :, ikp), 0, send_reqs(n_send), tag=n_kp_dos + ikp)
            END IF
         END IF

         CALL cp_cfm_release(floquet_matrix)
         CALL cp_cfm_release(cfm_eigenvectors)
      END DO

      DEALLOCATE (e_k, de_dk, dipole)

      ! Phase 3 (finish streaming): drain the last transfers and flush any remaining k-points.
      ! No global collective is used - senders wait only on their own isends, the collector on
      ! its own irecvs; once the irecvs complete every result is present, so all k-points can be
      ! marked received and the remaining prefix written in order.
      IF (am_collector) THEN
         IF (n_recv > 0) CALL mp_waitall(recv_reqs(1:n_recv))
         received(:) = .TRUE.
         CALL stream_poll_and_write()
         IF (unit_nr > 0) WRITE (unit_nr, '(A)') ""
         DEALLOCATE (recv_reqs, req_ikp)
      ELSE
         IF (n_send > 0) CALL mp_waitall(send_reqs(1:n_send))
         DEALLOCATE (send_reqs)
      END IF
      CALL comm_stream%free()

      ! Note: the distributed cfm matrices (floquet_matrix, cfm_eigenvectors)
      ! are released inside the k-point loop above; the distributed band arrays (e_k, de_dk,
      ! dipole) are released right after it. Only the persistent objects remain here.
      DEALLOCATE (eigenvalues, quasi_energies, a_k, all_quasi_energies, all_a_k, group_distribution)
      DEALLOCATE (e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, received, remaining)
      CALL cp_fm_struct_release(floquet_struct)
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL mp_para_env_release(para_env_sub)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief Collector-only, purely local (no MPI collective): reap any completed result messages
!>        with a non-blocking mp_testany, then write every newly-complete k-point in strict
!>        k-order, advancing the high-water mark max_printed_ikp. Host-associates the streaming
!>        state (received, remaining, req_ikp, recv_reqs, n_recv, max_printed_ikp) and the
!>        scratch/output buffers of the enclosing routine.
! **************************************************************************************************
      SUBROUTINE stream_poll_and_write()
      INTEGER                                            :: i_done, jkp, jspin
      LOGICAL                                            :: any_done

         ! Reap whatever has arrived so far (non-blocking); a k-point is received once all of
         ! its streams (quasi and/or DOS) have completed. mp_testany reports flag=.TRUE. with an
         ! out-of-range index (MPI_UNDEFINED) when no active request is left, so guard the index.
         IF (n_recv > 0) THEN
            DO
               CALL mp_testany(recv_reqs(1:n_recv), completed=i_done, flag=any_done)
               IF (.NOT. any_done) EXIT
               IF (i_done < 1 .OR. i_done > n_recv) EXIT
               jkp = req_ikp(i_done)
               remaining(jkp) = remaining(jkp) - 1
               IF (remaining(jkp) == 0) received(jkp) = .TRUE.
            END DO
         END IF

         ! Write the contiguous run of completed k-points that starts just after the last one
         ! written (the output files' REPLACE/APPEND logic requires strict k-order).
         DO WHILE (max_printed_ikp < n_kp_dos)
            IF (.NOT. received(max_printed_ikp + 1)) EXIT
            jkp = max_printed_ikp + 1
            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, nkp_start + jkp)
            DO jspin = 1, n_spin
               quasi_energies(:) = all_quasi_energies(:, jspin, jkp)
               CALL write_quasi_energy(bs_env, jspin, jkp, xkp, quasi_energies)
               a_k(:) = all_a_k(:, jspin, jkp)
               CALL write_floquet_spectral_function(qs_env, bs_env, jspin, jkp, xkp, a_k)
               IF (unit_nr > 0) WRITE (unit_nr, '(2X,I6,T12,I4,T18,3F12.6)') jkp, jspin, xkp(1:3)
            END DO
            max_printed_ikp = jkp
         END DO

      END SUBROUTINE stream_poll_and_write

   END SUBROUTINE floquet

END MODULE floquet_main

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Floquet-Bloch band-structure calculations: builds and diagonalises the
!>        truncated Floquet-Bloch Hamiltonian H_F(k) for each DOS k-point and writes the m=0
!>        band structure, the quasi-energies and the Floquet density of states.
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
   USE floquet_utils,                   ONLY: build_floquet_matrix,&
                                              calculate_floquet_observables,&
                                              check_floquet_convergence,&
                                              floquet_determine_subgroup_size,&
                                              floquet_precompute_kp_data,&
                                              floquet_sector_weights,&
                                              write_floquet_results
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_any_source,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_request_type
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

   PUBLIC :: floquet

CONTAINS

! **************************************************************************************************
!> \brief Computes the Floquet-Bloch band structure, density of states and quasi-energies.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE floquet(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'floquet'
      INTEGER, PARAMETER                                 :: floquet_done_tag = 1

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: dipole_kp_spin
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: dipole
      INTEGER :: group_size, handle, ikp, ispin, loc_idx, local_src, max_f_index, my_group, n_E, &
         n_f_size, n_fbands, n_kp_dos, n_own, n_received, n_remote, n_spin, nao, next_to_print, &
         ngroups, nkp, nkp_start, owner, stride_kp, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: group_distribution
      INTEGER, DIMENSION(1)                              :: recv_buf
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: ikp_done
      REAL(KIND=dp)                                      :: energy_step, energy_window
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a_k, e_k_kp_spin, eigenvalues, &
                                                            m0_energies, m0_weights, &
                                                            quasi_energies, w0, wE
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: de_dk_kp_spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: all_a_k, all_m0_energies, &
                                                            all_m0_weights, all_quasi_energies, e_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: de_dk
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors, floquet_matrix
      TYPE(cp_fm_struct_type), POINTER                   :: floquet_struct
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(mp_request_type)                              :: recv_request

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

      n_E = NINT(energy_window/energy_step)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(T2,A)') ' '
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A,A51,A27)') '-', 'FLOQUET BANDSTRUCTURE CALCULATION', '-'
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A)') ' '

         WRITE (unit_nr, '(T2,A,T67,ES12.2E2)') &
            "FLOQUET PARAMETERS                 Amplitude [V/m]:", bs_env%floquet_amplitude*evolt/a_bohr
         WRITE (unit_nr, '(T2,A,T67,F12.4)') &
            "                                   Frequency [eV]:", bs_env%floquet_omega*evolt
         WRITE (unit_nr, '(T2,A)') "                                   "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T55,3F8.4)') &
            "                                   Polarisation:", bs_env%floquet_polarisation(1:3)
         WRITE (unit_nr, '(T2,A,T55,3F8.4)') &
            "                                   Phase offsets:", pi*bs_env%floquet_phi(1:3)
         WRITE (unit_nr, '(T2,A)') "                                   "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T67,I12)') &
            "                                   Max Floquet index:", bs_env%max_floquet_index
         WRITE (unit_nr, '(T2,A,T67,I12)') &
            "                                   Floquet Hamiltonian Size:", n_f_size
         WRITE (unit_nr, '(T2,A)') "                                   "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T67,F12.4)') &
            "                                   Energy window [eV]:", bs_env%energy_window_floquet*evolt
         WRITE (unit_nr, '(T2,A,T67,F12.4)') &
            "                                   Energy step [eV]:", bs_env%energy_step_floquet*evolt
         WRITE (unit_nr, '(T2,A,T67,F12.4)') &
            "                                   Broadening [eV]:", bs_env%broadening_floquet*evolt
         WRITE (unit_nr, '(T2,A)') "                                   "//REPEAT("-", 44)
         WRITE (unit_nr, '(A)') ""

         WRITE (unit_nr, '(T2,A)') &
            "We construct the Floquet-Bloch Hamiltonian and diagonalise it. Projecting the"
         WRITE (unit_nr, '(T2,A)') &
            "eigenvectors onto the Floquet sectors gives the weights w = Σ_n |<n,m|α>|²,"
         WRITE (unit_nr, '(T2,A)') &
            "from which all of the following are obtained."
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(T2,A)') &
            "The m=0 bands are the eigenvectors with the largest central-sector weight; they"
         WRITE (unit_nr, '(T2,A)') &
            "reduce to the equilibrium bands at zero field and are stored, with their"
         WRITE (unit_nr, '(T2,A)') &
            " weights, in FLOQUET_BANDSTRUCTURE.bs"
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(T2,A)') &
            "Folding those bands into the first Floquet Brillouin zone, relative to the VBM,"
         WRITE (unit_nr, '(T2,A)') &
            "gives the quasi-energies stored in QUASI_ENERGIES.bs"
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

      ! Compute the dipole and band data for ALL DOS k-points once, distributed across ranks.
      CALL floquet_precompute_kp_data(qs_env, bs_env, e_k, de_dk, dipole)

      ALLOCATE (eigenvalues(n_f_size), quasi_energies(nao), a_k(n_E), w0(n_f_size), wE(n_f_size))
      ALLOCATE (m0_energies(nao), m0_weights(nao), e_k_kp_spin(nao), de_dk_kp_spin(3, nao))
      ALLOCATE (dipole_kp_spin(3, nao, nao))

      ALLOCATE (all_quasi_energies(nao, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_a_k(n_E, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_m0_energies(nao, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_m0_weights(nao, n_spin, n_kp_dos), source=0.0_dp)

      ! We split the global communicator into subgroups with the Floquet Hamiltonian distributed
      ! across the processes in each subgroup and the k-points distributed across subgroups.
      group_size = floquet_determine_subgroup_size(qs_env, bs_env)

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
         WRITE (unit_nr, '(/,T2,A,I5,A,I1,A)') "FLOQUET CALCULATIONS PROGRESS OUT OF", n_kp_dos, &
            " K-POINTS AND ", n_spin, " SPINS"
         WRITE (unit_nr, '(T2,A)') "KPOINT  SPIN"
      END IF

      NULLIFY (blacs_env_sub)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      ! The struct (and its BLACS context) is k-point independent but the cfm matrices are
      ! (de)allocated per k-point.
      CALL cp_fm_struct_create(floquet_struct, &
                               nrow_global=n_f_size, &
                               ncol_global=n_f_size, &
                               context=blacs_env_sub, &
                               para_env=para_env_sub)

      ! Ordered progress printing: the collector (global rank 0) tracks which DOS k-points have
      ! finished and prints them in ascending order, live from inside the loop.
      IF (para_env%is_source()) THEN
         ALLOCATE (ikp_done(n_kp_dos))
         ikp_done(:) = .FALSE.

         n_own = 0
         IF (para_env_sub%is_source()) THEN
            DO ikp = 1, n_kp_dos
               IF (group_distribution(MOD(ikp - 1, para_env%num_pe)) == my_group) n_own = n_own + 1
            END DO
         END IF
         n_remote = n_kp_dos - n_own

         ! Post the first non-blocking receive for a completion from any subgroup source; it is
         ! reposted after each one completes until all n_remote have been received.
         n_received = 0
         next_to_print = 1
         IF (n_remote > 0) THEN
            recv_buf(1) = 0
            CALL para_env%irecv(recv_buf, mp_any_source, recv_request, tag=floquet_done_tag)
         END IF
      END IF

      DO ikp = 1, n_kp_dos
         ! Collector: drain whatever completion messages have arrived without blocking,
         ! then print every contiguous finished k-point from the front.
         IF (para_env%is_source()) THEN
            DO WHILE (n_received < n_remote)
               IF (.NOT. recv_request%test()) EXIT
               ikp_done(recv_buf(1)) = .TRUE.
               n_received = n_received + 1
               IF (n_received < n_remote) THEN
                  recv_buf(1) = 0
                  CALL para_env%irecv(recv_buf, mp_any_source, recv_request, tag=floquet_done_tag)
               END IF
            END DO
            DO WHILE (next_to_print <= n_kp_dos)
               IF (.NOT. ikp_done(next_to_print)) EXIT
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,I6,I6)') next_to_print
               END IF
               next_to_print = next_to_print + 1
            END DO
         END IF

         ! Only the subgroup that owns this k-point's data works on it
         owner = MOD(ikp - 1, para_env%num_pe)
         IF (group_distribution(owner) /= my_group) CYCLE

         ! Identify the owner's rank within the subgroup.
         local_src = -1
         IF (para_env%mepos == owner) local_src = para_env_sub%mepos
         CALL para_env_sub%max(local_src)

         CALL cp_cfm_create(floquet_matrix, floquet_struct, set_zero=.TRUE.)
         CALL cp_cfm_create(cfm_eigenvectors, floquet_struct, set_zero=.TRUE.)

         DO ispin = 1, n_spin
            ! The subgroup owner copies this k-point's distributed band data to broadcast
            IF (para_env%mepos == owner) THEN
               loc_idx = CEILING(REAL(ikp)/para_env%num_pe)
               e_k_kp_spin(:) = e_k(ispin, loc_idx, :)
               de_dk_kp_spin(:, :) = de_dk(ispin, loc_idx, :, :)
               dipole_kp_spin(:, :, :) = dipole(ispin, loc_idx, :, :, :)
            END IF
            ! The data is then broadcast from it to all subgroup ranks
            CALL para_env_sub%bcast(e_k_kp_spin, local_src)
            CALL para_env_sub%bcast(de_dk_kp_spin, local_src)
            CALL para_env_sub%bcast(dipole_kp_spin, local_src)

            ! Build the Floquet-Bloch Hamiltonian H_F(k) on this subgroup's matrix from the
            ! broadcast band data of this k-point - no global communication involved.
            CALL build_floquet_matrix(bs_env, e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, &
                                      floquet_matrix)

            ! Diagonalize the Floquet-Bloch matrix H_F(k) F_{α,k} = ε_α F_{α,k} to obtain the
            ! eigenvalues ε_α and the Floquet-Bloch coefficients F^{nu}_{α,k}.
            CALL cp_cfm_set_all(cfm_eigenvectors, z_zero)
            CALL cp_cfm_heevd(floquet_matrix, cfm_eigenvectors, eigenvalues)

            ! Project every eigenvector onto the central (m=0) and outermost (m=±M) Floquet
            ! sectors, w_α^(m) = Σ_n |<n,m|α>|^2.
            CALL floquet_sector_weights(nao, max_f_index, cfm_eigenvectors, w0, wE)

            ! Check that MAX_FLOQUET_INDEX, the size of the truncated matrix H_F(k), is large
            ! enough: the m=0 states must carry no appreciable weight on the outermost rungs.
            CALL check_floquet_convergence(bs_env, w0, wE)

            ! Compute the DOS A(k,ω) = -(1/π) Im[Tr_KS G^R_{00}(k,ω)], the m=0 bands (the
            ! eigenvectors with the largest central-sector weight), and the quasi-energies.
            CALL calculate_floquet_observables(bs_env, eigenvalues, w0, &
                                               a_k, m0_energies, m0_weights, quasi_energies)

            ! Store this result on the subgroup source only, so the global sum below picks up
            ! each work item only once
            IF (para_env_sub%is_source()) THEN
               all_quasi_energies(:, ispin, ikp) = quasi_energies(:)
               all_a_k(:, ispin, ikp) = a_k(:)
               all_m0_energies(:, ispin, ikp) = m0_energies(:)
               all_m0_weights(:, ispin, ikp) = m0_weights(:)
            END IF
         END DO

         ! Report this k-point (all spins done) to the collector, global rank 0, which prints the
         ! k-points in ascending order
         IF (para_env_sub%is_source()) THEN
            IF (para_env%is_source()) THEN
               ikp_done(ikp) = .TRUE.
            ELSE
               CALL para_env%send(ikp, 0, floquet_done_tag)
            END IF
         END IF

         CALL cp_cfm_release(floquet_matrix)
         CALL cp_cfm_release(cfm_eigenvectors)
      END DO

      ! Collector: the loop is over, so any k-point still unfinished is one whose completion
      ! message has not yet been received. Block-receive the remaining ones.
      IF (para_env%is_source()) THEN
         DO WHILE (n_received < n_remote)
            CALL recv_request%wait()
            ikp_done(recv_buf(1)) = .TRUE.
            n_received = n_received + 1
            IF (n_received < n_remote) THEN
               recv_buf(1) = 0
               CALL para_env%irecv(recv_buf, mp_any_source, recv_request, tag=floquet_done_tag)
            END IF
         END DO
         DO WHILE (next_to_print <= n_kp_dos)
            IF (.NOT. ikp_done(next_to_print)) EXIT
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,I6,I6)') next_to_print
            END IF
            next_to_print = next_to_print + 1
         END DO
         DEALLOCATE (ikp_done)
      END IF

      DEALLOCATE (e_k, de_dk, dipole)

      ! Collect all results on the global communicator and write them
      CALL para_env%sum(all_quasi_energies)
      CALL para_env%sum(all_a_k)
      CALL para_env%sum(all_m0_energies)
      CALL para_env%sum(all_m0_weights)

      CALL write_floquet_results(bs_env, all_quasi_energies, all_a_k, &
                                 all_m0_energies, all_m0_weights)

      DEALLOCATE (eigenvalues, quasi_energies, a_k, all_quasi_energies, all_a_k, group_distribution)
      DEALLOCATE (w0, wE, m0_energies, m0_weights, all_m0_energies, all_m0_weights)
      DEALLOCATE (e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin)
      CALL cp_fm_struct_release(floquet_struct)
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL mp_para_env_release(para_env_sub)

      CALL timestop(handle)

   END SUBROUTINE floquet

END MODULE floquet_main

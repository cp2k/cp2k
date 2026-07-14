!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Driver for Floquet-Bloch band-structure calculations: builds and diagonalises the
!>        truncated Floquet-Bloch Hamiltonian H_F(k) for each DOS k-point and writes the m=0
!>        band structure, the quasi-energies and the Floquet density of states. The heavy
!>        physics lives in floquet_utils; this module only orchestrates the MPI-subgroup work
!>        distribution across k-points.
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
                                              calculate_floquet_observables,&
                                              check_floquet_convergence,&
                                              floquet_determine_subgroup_size,&
                                              floquet_precompute_kp_data,&
                                              floquet_sector_weights,&
                                              write_floquet_results
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
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
!> \brief Compute the Floquet-Bloch band structure in three phases. (1) Split the global
!>        communicator into subgroups sized to hold the distributed H_F(k) working set and to
!>        keep at least one subgroup per DOS k-point. (2) Precompute the band data for all DOS
!>        k-points once (distributed round-robin across ranks), then loop over k-points: each is
!>        handled by the subgroup owning its data, which broadcasts the band data internally,
!>        assembles and diagonalises H_F(k), and derives the sector weights and observables.
!>        (3) Gather the per-(k-point,spin) results across the global communicator and write the
!>        three output files on the I/O rank.
!> \param qs_env the Quickstep environment (SCF results, k-points, parallel environment)
!> \param bs_env the post-SCF band-structure environment holding all Floquet input parameters
! **************************************************************************************************
   SUBROUTINE floquet(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'floquet'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: dipole_kp_spin
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: dipole
      INTEGER :: group_size, handle, ikp, ispin, loc_idx, local_src, max_f_index, my_group, n_E, &
                 n_f_size, n_fbands, n_kp_dos, n_spin, nao, ngroups, nkp, nkp_start, owner, &
                 stride_kp, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: group_distribution
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
            "We construct the Floquet-Bloch Hamiltonian and diagonalise it. Projecting the"
         WRITE (unit_nr, '(T2,A)') &
            "eigenvectors onto the Floquet sectors gives the weights w = Σ_n |<n,m|α>|², from"
         WRITE (unit_nr, '(T2,A)') &
            "which all of the following are obtained."
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(T2,A)') &
            "The m=0 bands are the eigenvectors with the largest central-sector weight; they"
         WRITE (unit_nr, '(T2,A)') &
            "reduce to the equilibrium bands at zero field and are stored, with their weights,"
         WRITE (unit_nr, '(T2,A)') &
            "in FLOQUET_BANDSTRUCTURE.bs"
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

      ! Central-sector (m=0) and outermost-sector (m=±M) weights of the eigenvectors of H_F(k)
      ! for one (k-point, spin). Every Floquet observable below is a function of w0 and the
      ! eigenvalues alone; wE is consumed only by the truncation check.
      ALLOCATE (w0(n_f_size), wE(n_f_size))
      ALLOCATE (m0_energies(nao), m0_weights(nao))

      ! Per-single-(k-point,spin) broadcast buffers: the owner rank fills its slice and
      ! broadcasts them within the owning subgroup so every subgroup rank can build its piece
      ! of the Floquet matrix. These same buffers are passed to the Floquet matrix builder.
      ALLOCATE (e_k_kp_spin(nao))
      ALLOCATE (de_dk_kp_spin(3, nao))
      ALLOCATE (dipole_kp_spin(3, nao, nao))

      ! Per-(spin,k-point) result slots. Each is filled once by the owning subgroup's source
      ! rank; the full arrays are gathered across the global communicator once, after the loop.
      ALLOCATE (all_quasi_energies(nao, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_a_k(n_E, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_m0_energies(nao, n_spin, n_kp_dos), source=0.0_dp)
      ALLOCATE (all_m0_weights(nao, n_spin, n_kp_dos), source=0.0_dp)

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

            ! Project every eigenvector onto the central (m=0) and outermost (m=±M) Floquet
            ! sectors, w_α^(m) = Σ_n |<n,m|α>|^2. Read straight off the distributed local data,
            ! touching only the rows of those three blocks; one small real reduction follows.
            CALL floquet_sector_weights(nao, max_f_index, cfm_eigenvectors, w0, wE)

            ! Check that MAX_FLOQUET_INDEX, the size of the truncated matrix H_F(k), is large
            ! enough: the m=0 states must carry no appreciable weight on the outermost rungs.
            ! Aborts the whole job (mpi_abort) if that weight exceeds EPS_FLOQUET, so no
            ! reduction of the leakage across ranks is needed.
            CALL check_floquet_convergence(bs_env, w0, wE)

            ! Spectral function A_phys(k,ω) = -(1/π) Im[Tr_KS G^R_{00}(k,ω)], the m=0 bands (the
            ! nao eigenvectors with the largest central-sector weight, reducing to ε_nk at zero
            ! field), and the quasi-energies (those bands folded into the first Floquet BZ).
            CALL calculate_floquet_observables(bs_env, eigenvalues, w0, &
                                               a_k, m0_energies, m0_weights, quasi_energies)

            ! Store this (k-point, spin) result on the subgroup source only, so the global sum
            ! below picks up each work item exactly once.
            IF (para_env_sub%is_source()) THEN
               all_quasi_energies(:, ispin, ikp) = quasi_energies(:)
               all_a_k(:, ispin, ikp) = a_k(:)
               all_m0_energies(:, ispin, ikp) = m0_energies(:)
               all_m0_weights(:, ispin, ikp) = m0_weights(:)
            END IF
         END DO

         CALL cp_cfm_release(floquet_matrix)
         CALL cp_cfm_release(cfm_eigenvectors)
      END DO

      DEALLOCATE (e_k, de_dk, dipole)

      ! Phase 3 (gather + write): collect all results on the global communicator and write them
      ! out in k-point/spin order on the global I/O rank. Each work item was stored by exactly
      ! one subgroup source, so the sum reconstructs the full arrays without double counting.
      CALL para_env%sum(all_quasi_energies)
      CALL para_env%sum(all_a_k)
      CALL para_env%sum(all_m0_energies)
      CALL para_env%sum(all_m0_weights)

      CALL write_floquet_results(bs_env, all_quasi_energies, all_a_k, &
                                 all_m0_energies, all_m0_weights)

      ! Note: the distributed cfm matrices (floquet_matrix, cfm_eigenvectors) are released inside
      ! the k-point loop above; the distributed band arrays (e_k, de_dk, dipole) are released
      ! right after it. Only the persistent objects remain here.
      DEALLOCATE (eigenvalues, quasi_energies, a_k, all_quasi_energies, all_a_k, group_distribution)
      DEALLOCATE (w0, wE, m0_energies, m0_weights, all_m0_energies, all_m0_weights)
      DEALLOCATE (e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin)
      CALL cp_fm_struct_release(floquet_struct)
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL mp_para_env_release(para_env_sub)

      CALL timestop(handle)

   END SUBROUTINE floquet

END MODULE floquet_main

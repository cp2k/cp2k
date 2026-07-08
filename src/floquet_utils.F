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
MODULE floquet_utils
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_k_r_trafo_simple,         ONLY: replicate_rs_matrices,&
                                              rs_to_kp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE machine,                         ONLY: m_hostnm,&
                                              m_memory_details
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: geeig_right,&
                                              gemm_square
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: qs_moment_kpoints_deep
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'floquet_utils'

   ! Number of full Floquet-Hamiltonian copies (each n_f_size^2 complex(dp), 16 bytes per
   ! element) that the subgroup must be able to hold, distributed across its ranks. A few are
   ! live at once (floquet_matrix + cfm_eigenvectors + eigensolver scratch); we size for 10 to
   ! leave headroom.
   INTEGER, PARAMETER, PRIVATE :: n_fm_work_copies = 10

   ! Fraction of the per-rank free memory the Floquet working set is allowed to occupy. The rest
   ! is a reserve for everything the working-set estimate does NOT count: the SCF/post-SCF
   ! baseline, the ScaLAPACK eigensolver's internal workspace, the distributed dipole/band data,
   ! and the replicated real-space arrays built transiently in the precompute.
   REAL(KIND=dp), PARAMETER, PRIVATE :: mem_fill_fraction = 0.4_dp

   ! Public subroutines (the block builders, calculate_epsilon_derivative and
   ! floquet_ranks_per_node are private helpers used only within this module)
   PUBLIC :: build_floquet_matrix_local, &
             floquet_precompute_kp_data, &
             floquet_determine_subgroup_size, &
             check_floquet_conversion, &
             calculate_floquet_spectral_function, &
             get_floquet_energies, &
             write_floquet_results

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xkp ...
!> \param e_k ...
!> \param de_dk ...
!> \param do_parallel option to distribute the results e_k/de_dk across MPI ranks. When
!>        .TRUE. k-point ikp is computed and stored only on rank MOD(ikp-1,num_pe), at the
!>        compressed index CEILING(ikp/num_pe). This MUST match the layout used by
!>        qs_moment_kpoints_deep so that e_k/de_dk co-locate with the dipole of the same
!>        k-point. Default .FALSE. -> every rank computes/stores all k-points (replicated).
! **************************************************************************************************
   SUBROUTINE calculate_epsilon_derivative(qs_env, xkp, e_k, de_dk, do_parallel)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: e_k
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), OPTIONAL                 :: de_dk
      LOGICAL, OPTIONAL                                  :: do_parallel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_epsilon_derivative'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: C_dH_C, C_dS_C, C_k, dH_dk_i, dS_dk_i, &
                                                            H_k, S_k
      INTEGER                                            :: handle, i_dir, ikp, ispin, mepos, n, &
                                                            n_img_all, n_spin, nao, nkp, num_copy, &
                                                            num_pe
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_all
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_all
      LOGICAL                                            :: my_do_parallel, present_dedk, present_ek
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: H_rs, S_rs
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp
      TYPE(kpoint_type), POINTER                         :: kpoints_all, kpoints_scf
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all

      CALL timeset(routineN, handle)

      present_ek = PRESENT(e_k)  ! calculate band energies, ε_k for all kpoints xkp
      present_dedk = PRESENT(de_dk)  ! calculate derivative, ∇_k ε_k of band energies for all kpoints xkp

      IF (.NOT. (present_ek .OR. present_dedk)) CPABORT("Subroutine needs either e_k or de_dk")

      my_do_parallel = .FALSE.
      IF (PRESENT(do_parallel)) my_do_parallel = do_parallel

      CALL get_qs_env(qs_env, &
                      matrix_ks_kp=matrix_ks_kp, &
                      matrix_s_kp=matrix_s_kp, &
                      sab_all=sab_all, &
                      cell=cell, &
                      kpoints=kpoints_scf, &
                      para_env=para_env, &
                      mos=mos)

      CALL get_mo_set(mo_set=mos(1), nao=nao)
      CALL get_cell(cell=cell, h=hmat)

      n_spin = SIZE(matrix_ks_kp, 1)
      nkp = SIZE(xkp, 2)

      ! Distribution of k-points across ranks (mirrors qs_moment_kpoints_deep). With
      ! do_parallel=.FALSE. every rank computes all k-points (replicated). With .TRUE.
      ! rank mepos computes only k-points with MOD(ikp-1,num_pe)==mepos, stored at the
      ! compressed index CEILING(ikp/num_pe).
      mepos = 0
      num_pe = 1
      num_copy = nkp
      IF (my_do_parallel) THEN
         mepos = para_env%mepos
         num_pe = para_env%num_pe
         num_copy = CEILING(REAL(nkp)/num_pe)
      END IF

      ! create kpoint environment kpoints_all which contains all neighbor cells R
      ! without considering any lattice symmetry
      NULLIFY (kpoints_all)
      CALL kpoint_create(kpoints_all)
      CALL kpoint_init_cell_index(kpoints_all, sab_all, para_env, n_img_all)
      CALL get_kpoint_info(kpoints_all, cell_to_index=cell_to_index_all, index_to_cell=index_to_cell_all)

      ALLOCATE (S_rs(1, nao, nao, n_img_all), H_rs(n_spin, nao, nao, n_img_all), source=0.0_dp)

      ! Convert real-space dbcsr matrices into arrays
      CALL replicate_rs_matrices(matrix_s_kp, kpoints_scf, S_rs, cell_to_index_all)
      CALL replicate_rs_matrices(matrix_ks_kp, kpoints_scf, H_rs, cell_to_index_all)

      IF (present_dedk) ALLOCATE (de_dk(n_spin, num_copy, 3, nao), source=0.0_dp)
      IF (present_ek) ALLOCATE (e_k(n_spin, num_copy, nao), source=0.0_dp)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(ikp, ispin, S_k, H_k, eigenvals, C_k, dS_dk_i, dH_dk_i, C_dS_C, C_dH_C) &
!$OMP SHARED(nao, n_spin, de_dk, e_k, present_ek, present_dedk, mepos, num_pe, &
!$OMP nkp, xkp, S_rs, H_rs, index_to_cell_all, hmat)
      IF (present_dedk) ALLOCATE (dS_dk_i(nao, nao), C_dS_C(nao, nao), &
                                  dH_dk_i(nao, nao), C_dH_C(nao, nao), source=z_zero)
      ALLOCATE (C_k(nao, nao), S_k(nao, nao), H_k(nao, nao), source=z_zero)
      ALLOCATE (eigenvals(nao), source=0.0_dp)
!$OMP DO COLLAPSE(2)
      DO ispin = 1, n_spin
         DO ikp = 1, nkp
            IF (MOD(ikp - 1, num_pe) /= mepos) CYCLE

            ! S^R -> S(k), H^R -> H(k)
            S_k = 0
            H_k = 0
            CALL rs_to_kp(S_rs(1, :, :, :), S_k, index_to_cell_all, xkp(:, ikp))
            CALL rs_to_kp(H_rs(ispin, :, :, :), H_k, index_to_cell_all, xkp(:, ikp))

            ! Diagonalize H(k)C(k) = S(k)C(k)ε(k)
            CALL geeig_right(H_k, S_k, eigenvals, C_k)
            IF (present_ek) e_k(ispin, CEILING(REAL(ikp)/num_pe), :) = eigenvals(:)

            IF (present_dedk) THEN
               ! Evaluate the derivatives using
               ! ∇ ε_k = C^H(k) ∇ H_k C(k) - ε_k C^H(k) ∇ S_k C(k)
               DO i_dir = 1, 3
                  CALL rs_to_kp(S_rs(1, :, :, :), dS_dk_i, index_to_cell_all, xkp(:, ikp), i_dir, hmat)
                  CALL rs_to_kp(H_rs(ispin, :, :, :), dH_dk_i, index_to_cell_all, xkp(:, ikp), i_dir, hmat)

                  CALL gemm_square(C_k, 'C', dS_dk_i, 'N', C_k, 'N', C_dS_C)
                  CALL gemm_square(C_k, 'C', dH_dk_i, 'N', C_k, 'N', C_dH_C)

                  DO n = 1, nao
                     de_dk(ispin, CEILING(REAL(ikp)/num_pe), i_dir, n) = &
                        DBLE(C_dH_C(n, n)) - DBLE(eigenvals(n)*C_dS_C(n, n))
                  END DO
               END DO
            END IF
         END DO
      END DO
!$OMP END DO
      IF (present_dedk) DEALLOCATE (dS_dk_i, C_dS_C, dH_dk_i, C_dH_C)
      DEALLOCATE (S_k, H_k, C_k, eigenvals)
!$OMP END PARALLEL
      DEALLOCATE (S_rs, H_rs)
      CALL kpoint_release(kpoints_all)

      CALL timestop(handle)

   END SUBROUTINE calculate_epsilon_derivative

! **************************************************************************************************
!> \brief Momentum matrix elements p_nm(k) for one k-point and spin, from precomputed band data:
!>        p_nm = i d_nm (ε_n - ε_m) + ∇_k ε_n δ_nm (= <ψ_n|-iħ∇_r|ψ_m>). No communication.
!> \param e_k_kp_spin band energies for this k-point/spin (nao)
!> \param de_dk_kp_spin band-energy k-derivatives for this k-point/spin (3, nao)
!> \param dipole_kp_spin dipole matrix elements for this k-point/spin (3, nao, nao)
!> \param momentum momentum matrix elements p (3, nao, nao) for this k-point/spin
! **************************************************************************************************
   SUBROUTINE build_momentum_matrix(e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, momentum)
      REAL(KIND=dp), DIMENSION(:)                        :: e_k_kp_spin
      REAL(KIND=dp), DIMENSION(:, :)                     :: de_dk_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :, :)               :: dipole_kp_spin, momentum

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_momentum_matrix'

      INTEGER                                            :: handle, i, i_dir, j, nao

      CALL timeset(routineN, handle)

      nao = SIZE(e_k_kp_spin)

      ! We calculate momentum matrix elements p_nm = <ψ_n|-iħ∇_r|ψ_m>
      ! p_nm = <u_n|(ħk - iħ∇_r)|u_m>
      ! p_nm = i d_nm (ε_n - ε_m) + ∇_k ε_n δ_nm

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i_dir, i, j) &
!$OMP SHARED(nao, momentum, dipole_kp_spin, e_k_kp_spin, de_dk_kp_spin)
!$OMP DO COLLAPSE(3)
      DO i_dir = 1, 3
         DO i = 1, nao
            DO j = 1, nao
               IF (j == i) THEN
                  momentum(i_dir, i, j) = de_dk_kp_spin(i_dir, i)
               ELSE
                  momentum(i_dir, i, j) = gaussi*(e_k_kp_spin(i) - e_k_kp_spin(j))* &
                                          dipole_kp_spin(i_dir, i, j)
               END IF
            END DO
         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE build_momentum_matrix

! **************************************************************************************************
!> \brief Off-diagonal coupling block of H_F(k) for one k-point and spin, from precomputed data:
!>        off_diag_uv = Σ_α p_uv^α(k) · E_factor(α), E_factor(α) = (i E(α) e^{iφ(α)})/(2ω).
!>        The block is accumulated into off_diag_m (which must be pre-zeroed by the caller).
!> \param bs_env ...
!> \param e_k_kp_spin band energies for this k-point/spin (nao)
!> \param de_dk_kp_spin band-energy k-derivatives for this k-point/spin (3, nao)
!> \param dipole_kp_spin dipole matrix elements for this k-point/spin (3, nao, nao)
!> \param off_diag_m off-diagonal coupling block (nao, nao)
! **************************************************************************************************
   SUBROUTINE build_off_diagonal_matrix(bs_env, e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, off_diag_m)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:)                        :: e_k_kp_spin
      REAL(KIND=dp), DIMENSION(:, :)                     :: de_dk_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :, :)               :: dipole_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: off_diag_m

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_off_diagonal_matrix'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: momentum
      COMPLEX(KIND=dp), DIMENSION(3)                     :: efactor
      INTEGER                                            :: handle, i, i_dir, j, nao
      REAL(KIND=dp)                                      :: amplitude, omega
      REAL(KIND=dp), DIMENSION(3)                        :: e_vec, phi, polarisation

      CALL timeset(routineN, handle)

      ! Builds the matrix that occupies the off diagonal blocks in the Floquet Matrix H_F

      nao = SIZE(e_k_kp_spin)

      polarisation(:) = bs_env%floquet_polarisation(:)
      IF (SQRT(SUM(polarisation**2)) < EPSILON(0.0_dp)) THEN
         CPABORT("Invalid (too small) polarisation vector specified for POLARISATION_PUMP")
      END IF

      amplitude = bs_env%floquet_amplitude
      e_vec(:) = amplitude*polarisation
      phi(:) = pi*bs_env%floquet_phi(:)
      omega = bs_env%floquet_omega

      ALLOCATE (momentum(3, nao, nao), source=z_zero)
      CALL build_momentum_matrix(e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, momentum)

      ! E_factor(α) = (i E(α) exp(i·φ(α)))/(2ω) where α = x,y,z
      DO i_dir = 1, 3
         efactor(i_dir) = gaussi*e_vec(i_dir)*CMPLX(COS(phi(i_dir)), SIN(phi(i_dir)), KIND=dp)/(2*omega)
      END DO

      DO i_dir = 1, 3
         DO i = 1, nao
            DO j = 1, nao
               off_diag_m(i, j) = off_diag_m(i, j) + &
                                  momentum(i_dir, i, j)*efactor(i_dir)
            END DO
         END DO
      END DO
      DEALLOCATE (momentum)

      CALL timestop(handle)

   END SUBROUTINE build_off_diagonal_matrix

! **************************************************************************************************
!> \brief Diagonal block of H_F(k) for Floquet sector f_index: equilibrium band energies shifted
!>        by f_index·ħΩ, diag_e(n,n) = ε_{nk} + f_index·ħΩ. Only the diagonal is written.
!> \param bs_env ...
!> \param e_k_kp_spin band energies for this k-point/spin (nao)
!> \param f_index Floquet sector index
!> \param diag_e diagonal block (nao, nao)
! **************************************************************************************************
   SUBROUTINE build_diagonal_matrix(bs_env, e_k_kp_spin, f_index, diag_e)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:)                        :: e_k_kp_spin
      INTEGER                                            :: f_index
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: diag_e

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_diagonal_matrix'

      INTEGER                                            :: handle, i, nao
      REAL(KIND=dp)                                      :: omega

      CALL timeset(routineN, handle)

      ! Builds the equilibrium-energy diagonal block of Floquet sector f_index

      nao = SIZE(e_k_kp_spin)
      omega = bs_env%floquet_omega
      DO i = 1, nao
         diag_e(i, i) = e_k_kp_spin(i) + f_index*omega
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_diagonal_matrix

! **************************************************************************************************
!> \brief Communication-free assembly of the Floquet-Bloch Hamiltonian H_F(k) for a single
!>        k-point and spin channel, from band data already broadcast to every rank of the
!>        owning subgroup (e_k_kp, de_dk_kp, dipole_kp). Performs no global MPI collectives, so
!>        it runs independently inside the subgroup whose BLACS context owns floquet_matrix.
!>        Assembled from the modular block builders build_off_diagonal_matrix /
!>        build_diagonal_matrix (which in turn use build_momentum_matrix).
!> \param bs_env ...
!> \param e_k_kp_spin band energies for this k-point/spin (nao)
!> \param de_dk_kp_spin band-energy k-derivatives for this k-point/spin (3, nao)
!> \param dipole_kp_spin dipole matrix elements for this k-point/spin (3, nao, nao)
!> \param floquet_matrix the (subgroup-distributed) Floquet matrix to fill
! **************************************************************************************************
   SUBROUTINE build_floquet_matrix_local(bs_env, e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, floquet_matrix)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:)                        :: e_k_kp_spin
      REAL(KIND=dp), DIMENSION(:, :)                     :: de_dk_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :, :)               :: dipole_kp_spin
      TYPE(cp_cfm_type)                                  :: floquet_matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_floquet_matrix_local'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: conj_off_diag_m, diag_e, off_diag_m
      INTEGER                                            :: f_index, handle, i, i_f, max_f_index, &
                                                            n_fbands, nao

      CALL timeset(routineN, handle)

      nao = SIZE(e_k_kp_spin)
      max_f_index = bs_env%max_floquet_index
      n_fbands = 1 + 2*max_f_index

      ! Creates the floquet matrix H_F by placing the diagonal and off diagonal blocks
      ALLOCATE (diag_e(nao, nao), source=z_zero)
      ALLOCATE (off_diag_m(nao, nao), source=z_zero)
      ALLOCATE (conj_off_diag_m(nao, nao), source=z_zero)

      CALL build_off_diagonal_matrix(bs_env, e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, off_diag_m)
      conj_off_diag_m(:, :) = CONJG(TRANSPOSE(off_diag_m(:, :)))

      floquet_matrix%local_data(:, :) = z_zero
      DO i = 1, n_fbands
         i_f = 1 + (i - 1)*nao
         f_index = i - max_f_index - 1
         CALL build_diagonal_matrix(bs_env, e_k_kp_spin, f_index, diag_e)
         CALL cp_cfm_set_submatrix(floquet_matrix, diag_e, i_f, i_f)
         IF (i > 1) THEN
            CALL cp_cfm_set_submatrix(floquet_matrix, off_diag_m, i_f, i_f - nao)
            CALL cp_cfm_set_submatrix(floquet_matrix, conj_off_diag_m, i_f - nao, i_f)
         END IF
      END DO
      DEALLOCATE (diag_e, off_diag_m, conj_off_diag_m)

      CALL timestop(handle)
   END SUBROUTINE build_floquet_matrix_local

! **************************************************************************************************
!> \brief Precompute, on the global communicator, the band quantities needed to assemble the
!>        Floquet-Bloch Hamiltonian for all DOS k-points (all spins), distributed across ranks.
!>
!>        The expensive part of building H_F(k) is the conversion of the distributed
!>        real-space DBCSR matrices into dense real-space arrays (done inside
!>        calculate_epsilon_derivative and qs_moment_kpoints_deep via the global collective
!>        group%sum). Those routines are therefore called by ALL ranks of the global
!>        communicator (do NOT call this from inside a subgroup).
!>
!>        Both are called with do_parallel=.TRUE., so the RESULTS are DISTRIBUTED, not
!>        replicated: k-point ikp is computed and stored only on rank MOD(ikp-1,num_pe), at
!>        the compressed local index CEILING(ikp/num_pe). e_k, de_dk and dipole all use the
!>        same layout, so a rank that owns a k-point owns all three for that k-point. This is
!>        what lets the caller hold only ~nkp/num_pe k-points per rank (the full set cannot
!>        be stored for large systems) and route each k-point to the subgroup that owns it.
!>
!>        The large real-space arrays (H_rs, S_rs, D_rs) are allocated and freed inside these
!>        routines, so they do not coexist with the distributed Floquet matrices the caller
!>        allocates only after this call returns.
!> \param qs_env ...
!> \param bs_env ...
!> \param e_k band energies ε_{nk} (nspin, num_copy, nao), num_copy = CEILING(nkp_DOS/num_pe)
!> \param de_dk k-derivatives ∇_k ε_{nk} (nspin, num_copy, 3, nao)
!> \param dipole dipole matrix elements d_nm(k) (nspin, num_copy, 3, nao, nao)
! **************************************************************************************************
   SUBROUTINE floquet_precompute_kp_data(qs_env, bs_env, e_k, de_dk, dipole)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: e_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: de_dk
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: dipole

      CHARACTER(LEN=*), PARAMETER :: routineN = 'floquet_precompute_kp_data'

      INTEGER                                            :: handle, ikp, n_kp_dos, nkp_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp_all

      CALL timeset(routineN, handle)

      nkp_start = bs_env%nkp_only_DOS
      n_kp_dos = bs_env%nkp_bs_and_DOS - nkp_start

      ! Collect the coordinates of all DOS k-points (global storage index nkp_start + local).
      ALLOCATE (xkp_all(3, n_kp_dos))
      DO ikp = 1, n_kp_dos
         xkp_all(:, ikp) = bs_env%kpoints_DOS%xkp(1:3, nkp_start + ikp)
      END DO

      ! Distribute the results across ranks (one k-point per rank, round-robin by global
      ! rank). Both routines must use the SAME distribution so e_k/de_dk co-locate with the
      ! dipole of each k-point.
      CALL calculate_epsilon_derivative(qs_env, xkp_all, e_k=e_k, de_dk=de_dk, do_parallel=.TRUE.)
      CALL qs_moment_kpoints_deep(qs_env, xkp_all, dipole, do_parallel=.TRUE.)

      DEALLOCATE (xkp_all)

      CALL timestop(handle)
   END SUBROUTINE floquet_precompute_kp_data

! **************************************************************************************************
!> \brief Number of MPI ranks that share a physical node, determined by splitting the given
!>        communicator into node-local communicators keyed on the hostname.
!>        Each rank's color is the global rank index of the first rank with the same hostname,
!>        so same-node ranks always get the same color with no hash collisions possible.
!> \param para_env the (global) parallel environment
!> \return the number of ranks on the calling rank's node
! **************************************************************************************************
   FUNCTION floquet_ranks_per_node(para_env) RESULT(ranks_per_node)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: ranks_per_node

      CHARACTER(LEN=default_string_length)               :: hostname
      INTEGER                                            :: color, i, nranks, str_len
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: hostnames_ascii, my_hostname_ascii
      TYPE(mp_comm_type)                                 :: node_comm

      CALL m_hostnm(hostname)

      ! All-gather every rank's hostname as ASCII integers, then assign color = the index
      ! of the first rank that shares this rank's hostname. Ranks on the same node get
      ! identical color; ranks on different nodes get different colors. No hash needed.
      str_len = default_string_length
      nranks = para_env%num_pe
      ALLOCATE (my_hostname_ascii(str_len), hostnames_ascii(nranks*str_len))
      DO i = 1, str_len
         my_hostname_ascii(i) = ICHAR(hostname(i:i))
      END DO
      hostnames_ascii(:) = 0
      hostnames_ascii(para_env%mepos*str_len + 1:(para_env%mepos + 1)*str_len) = my_hostname_ascii(:)
      CALL para_env%sum(hostnames_ascii)

      color = para_env%mepos
      DO i = 0, para_env%mepos - 1
         IF (ALL(hostnames_ascii(i*str_len + 1:(i + 1)*str_len) == my_hostname_ascii(:))) THEN
            color = i
            EXIT
         END IF
      END DO
      DEALLOCATE (my_hostname_ascii, hostnames_ascii)

      CALL node_comm%from_split(comm=para_env, color=color, key=para_env%mepos)
      ranks_per_node = MAX(1, node_comm%num_pe)
      CALL node_comm%free()

   END FUNCTION floquet_ranks_per_node

! **************************************************************************************************
!> \brief Choose the number of MPI ranks per subgroup as the larger of two floors:
!>          1. Memory floor G_mem: the subgroup must jointly hold the Floquet working set,
!>             taken to be n_fm_work_copies copies of the n_f_size^2 complex(dp) Floquet
!>             Hamiltonian (needed_bytes = 16 * n_fm_work_copies * n_f_size^2), so
!>             G_mem = ceil(needed_bytes / usable per-rank memory), where the usable memory is
!>             mem_fill_fraction of the per-rank FREE memory (MemLikelyFree / ranks-per-node) -
!>             the rest is kept as a reserve for the baseline, eigensolver scratch, band data,
!>             and the precompute's replicated arrays, none of which are in needed_bytes.
!>          2. K-point floor G_load = floor(num_pe / n_kp_dos): the LARGEST G that still
!>             yields at least n_kp_dos subgroups (floor(num_pe/G) >= n_kp_dos). When ranks
!>             outnumber k-points this hands every k-point its own (larger, hence faster)
!>             subgroup instead of leaving spare ranks idle; the caller stride-interleaves the
!>             split so the data owners actually spread across those subgroups. When
!>             n_kp_dos >= num_pe this floor is 0, so G_mem alone applies (unchanged behaviour).
!>             (floor, not ceil: ceil would give one subgroup too few, forcing an extra round.)
!> \param qs_env ...
!> \param bs_env ...
!> \return the chosen number of ranks per subgroup (1 .. num_pe)
! **************************************************************************************************
   FUNCTION floquet_determine_subgroup_size(qs_env, bs_env) RESULT(group_size)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: group_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'floquet_determine_subgroup_size'

      INTEGER                                            :: g_load, g_mem, handle, n_f_size, &
                                                            n_kp_dos, nao, ranks_per_node, unit_nr
      INTEGER(KIND=int_8)                                :: Buffers, Cached, MemFree, MemLikelyFree, &
                                                            MemTotal, needed_bytes, Slab, &
                                                            SReclaimable, usable_per_rank
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      unit_nr = bs_env%unit_nr
      n_f_size = nao*(1 + 2*bs_env%max_floquet_index)
      n_kp_dos = bs_env%nkp_bs_and_DOS - bs_env%nkp_only_DOS

      ! Memory the Floquet working set may use on one rank: a fraction of the per-rank FREE
      ! memory (node free memory / ranks per node). Free (not total) so the already-resident
      ! SCF baseline is excluded; the fraction leaves a reserve for the many consumers not in
      ! needed_bytes (eigensolver scratch, band data, replicated precompute arrays).
      CALL m_memory_details(MemTotal, MemFree, Buffers, Cached, Slab, SReclaimable, MemLikelyFree)
      ranks_per_node = floquet_ranks_per_node(para_env)
      usable_per_rank = INT(mem_fill_fraction*REAL(MemLikelyFree, dp), int_8)/INT(ranks_per_node, int_8)

      ! Total memory the subgroup must hold: n_fm_work_copies copies of the Floquet
      ! Hamiltonian, each n_f_size^2 complex(dp) numbers at 16 bytes each.
      needed_bytes = 16_int_8*INT(n_fm_work_copies, int_8)*INT(n_f_size, int_8)**2

      IF (usable_per_rank <= 0_int_8) THEN
         ! Memory info unavailable (e.g. no /proc on this platform): fall back to one
         ! subgroup spanning all ranks, which minimises the per-rank footprint.
         g_mem = para_env%num_pe
      ELSE
         ! Memory floor: smallest G with needed_bytes/G <= usable_per_rank. Clamp to num_pe in
         ! int_8 arithmetic before converting, to avoid overflowing a default integer.
         g_mem = INT(MIN((needed_bytes + usable_per_rank - 1_int_8)/usable_per_rank, &
                         INT(para_env%num_pe, int_8)))
      END IF

      ! K-point floor: floor(num_pe / n_kp_dos) is the largest G that still gives at least
      ! n_kp_dos subgroups (floor(num_pe/G) >= n_kp_dos), so each k-point can run in its own
      ! subgroup (one round) when memory allows, with the spare ranks growing each subgroup
      ! rather than sitting idle. Zero when n_kp_dos >= num_pe (then G_mem alone applies).
      g_load = para_env%num_pe/MAX(n_kp_dos, 1)

      group_size = MAX(g_mem, g_load)
      group_size = MIN(group_size, para_env%num_pe)
      group_size = MAX(group_size, 1)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') ""
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Detected MPI ranks per node:", ranks_per_node
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Free memory per rank:", &
            (MemLikelyFree/INT(ranks_per_node, int_8))/1048576_int_8, " MiB"
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Usable per rank (reserve applied):", &
            usable_per_rank/1048576_int_8, " MiB"
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Floquet copies needed:", n_fm_work_copies
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Total Floquet working set:", &
            needed_bytes/1048576_int_8, " MiB"
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Floquet working set per rank:", &
            (needed_bytes/INT(MAX(group_size, 1), int_8))/1048576_int_8, " MiB"
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Memory floor (ranks/subgroup):", g_mem
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | K-point floor (ranks/subgroup):", g_load
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Chosen MPI ranks per subgroup:", group_size
      END IF

      CALL timestop(handle)
   END FUNCTION floquet_determine_subgroup_size

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param cfm_eigenvectors ...
! **************************************************************************************************
   SUBROUTINE check_floquet_conversion(qs_env, bs_env, cfm_eigenvectors)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_floquet_conversion'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vector_block_0, vector_block_1
      INTEGER                                            :: handle, max_f_index, nao, start_row_0, &
                                                            start_row_1
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      max_f_index = bs_env%max_floquet_index

      start_row_0 = 1 + nao*max_f_index
      start_row_1 = 1 + nao*(max_f_index + 1)

      ALLOCATE (vector_block_0(3*nao, nao), vector_block_1(3*nao, nao), source=z_zero)
      CALL cp_cfm_get_submatrix(cfm_eigenvectors, vector_block_0, start_row_0, start_row_0)
      CALL cp_cfm_get_submatrix(cfm_eigenvectors, vector_block_1, start_row_1, start_row_1)
      IF (bs_env%eps_floquet > 0) THEN
         IF (ABS(MAXVAL(ABS(vector_block_0)) - MAXVAL(ABS(vector_block_1))) > bs_env%eps_floquet) THEN
            CPABORT("FLOQUET DIAGONALIZATION DID NOT CONVERGE. MAX_FLOQUET_INDEX TOO SMALL")
         END IF
      END IF
      DEALLOCATE (vector_block_0, vector_block_1)
      CALL timestop(handle)

   END SUBROUTINE check_floquet_conversion

! **************************************************************************************************
!> \brief k-resolved Floquet density of states A(k,E) = -(1/π) Im Tr_KS G^R_00(k,E) from the
!>        spectral representation of the retarded Green's function. With H_F = Σ_α λ_α |α><α|,
!>          A(E) = -(1/π) Σ_α w_α Im[ 1/(E + iη - λ_α) ],   η = broadening/2,
!>        where w_α = Σ_{n in central (f=0) Floquet sector} |<n|α>|^2 is the weight of Floquet
!>        eigenvector α on that sector. This reuses the eigenpairs from cp_cfm_heevd, so it
!>        avoids the per-energy dense matrix inversion of the full Floquet Hamiltonian.
!> \param qs_env ...
!> \param bs_env ...
!> \param eigenvalues Floquet eigenvalues λ_α (n_f_size)
!> \param cfm_eigenvectors Floquet eigenvectors |α> (n_f_size x n_f_size)
!> \param a_k the spectral function A(k,E) on the energy grid (n_E)
! **************************************************************************************************
   SUBROUTINE calculate_floquet_spectral_function(qs_env, bs_env, eigenvalues, cfm_eigenvectors, a_k)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_floquet_spectral_function'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: v0_block
      INTEGER                                            :: handle, i_E, max_f_index, n_E, n_f_size, &
                                                            nao, start_row_g00
      REAL(KIND=dp)                                      :: broad, E_min, energy, energy_step, mu
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weight
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos)

      CALL get_mo_set(mo_set=mos(1), mu=mu, nao=nao)

      broad = bs_env%broadening_floquet
      energy_step = bs_env%energy_step_floquet
      E_min = mu - bs_env%energy_window_floquet
      max_f_index = bs_env%max_floquet_index
      n_E = SIZE(a_k)
      n_f_size = SIZE(eigenvalues)
      start_row_g00 = 1 + nao*max_f_index

      ! Central-sector rows of the eigenvectors, <n|α> for n in the f=0 block and all α, then
      ! their per-eigenvector weight w_α = Σ_n |<n|α>|^2 (Σ_α w_α = nao by unitarity).
      ALLOCATE (v0_block(nao, n_f_size))
      CALL cp_cfm_get_submatrix(cfm_eigenvectors, v0_block, start_row_g00, 1)
      ALLOCATE (weight(n_f_size))
      weight(:) = SUM(ABS(v0_block)**2, DIM=1)

      DO i_E = 1, n_E
         energy = E_min + i_E*energy_step
         a_k(i_E) = -SUM(weight(:)*AIMAG(z_one/(energy + gaussi*broad/2.0_dp - eigenvalues(:))))/pi
      END DO

      DEALLOCATE (v0_block, weight)

      CALL timestop(handle)

   END SUBROUTINE calculate_floquet_spectral_function

! **************************************************************************************************
!> \brief Extract the nao quasi-energies from the full set of Floquet eigenvalues: the
!>        quasi-energies are the nao eigenvalues of smallest magnitude (folded into the first
!>        Floquet Brillouin zone -ħΩ/2 < ε ≤ ħΩ/2), returned ordered by their position in the
!>        spectrum.
!> \param bs_env ...
!> \param eigenvalues the full set of Floquet eigenvalues
!> \param floquet_bands the energies of the m=0 sector
!> \param quasi_energies the extracted quasi-energies
! **************************************************************************************************
   SUBROUTINE get_floquet_energies(bs_env, eigenvalues, floquet_bands, quasi_energies)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues, floquet_bands, &
                                                            quasi_energies

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_floquet_energies'

      INTEGER                                            :: handle, m, n, nao
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: work
      REAL(KIND=dp)                                      :: omega

      CALL timeset(routineN, handle)

      omega = bs_env%floquet_omega
      n = SIZE(eigenvalues)
      nao = SIZE(quasi_energies)
      m = (n/nao - 1)/2

      ALLOCATE (work(nao))
      floquet_bands(:) = eigenvalues(nao*m + 1:nao*(m + 1))
      quasi_energies(:) = floquet_bands(:) - omega*(0.5_dp + FLOOR(floquet_bands/omega))
      CALL sort(quasi_energies, nao, work)
      DEALLOCATE (work)

      CALL timestop(handle)

   END SUBROUTINE get_floquet_energies

! **************************************************************************************************
!> \brief Write all Floquet results in one pass: the quasi-energies to <floquet_qe_file>.bs and
!>        the k-resolved density of states to <floquet_dos_file>.out. Each file is opened once
!>        (REPLACE) and written for every DOS k-point and spin, on the global I/O rank.
!> \param qs_env ...
!> \param bs_env ...
!> \param all_floquet_bands ...
!> \param all_quasi_energies quasi-energies for every (band, spin, DOS k-point) (nao, n_spin, n_kp_dos)
!> \param all_a_k spectral function A(k,E) for every (energy, spin, DOS k-point) (n_E, n_spin, n_kp_dos)
! **************************************************************************************************
   SUBROUTINE write_floquet_results(qs_env, bs_env, all_floquet_bands, all_quasi_energies, all_a_k)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: all_floquet_bands, all_quasi_energies, &
                                                            all_a_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_floquet_results'

      CHARACTER(LEN=default_string_length)               :: fname
      INTEGER                                            :: bunit, handle, i, i_E, ikp_for_file, &
                                                            ispin, n_E, n_kp_dos, n_spin, nao, &
                                                            nkp_start, qunit, wunit
      REAL(KIND=dp)                                      :: E_min, energy, energy_step, mu
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      IF (bs_env%para_env%is_source()) THEN

         nao = SIZE(all_floquet_bands, 1)
         n_spin = SIZE(all_floquet_bands, 2)
         n_kp_dos = SIZE(all_floquet_bands, 3)
         n_E = SIZE(all_a_k, 1)
         nkp_start = bs_env%nkp_only_DOS

         CALL get_qs_env(qs_env, mos=mos)
         CALL get_mo_set(mo_set=mos(1), mu=mu)
         energy_step = bs_env%energy_step_floquet
         E_min = mu - bs_env%energy_window_floquet

         WRITE (fname, "(2A)") TRIM(bs_env%floquet_bs_file), ".out"
         CALL open_file(TRIM(fname), unit_number=bunit, file_status="REPLACE", file_action="WRITE")
         WRITE (bunit, "(A)") "# Floquet engineered bands of the m=0 sector"

         WRITE (fname, "(2A)") TRIM(bs_env%floquet_qe_file), ".bs"
         CALL open_file(TRIM(fname), unit_number=qunit, file_status="REPLACE", file_action="WRITE")
         WRITE (qunit, "(A)") "# Quasi-energies obtained by diagonalising the Floquet Hamiltonian"
         WRITE (qunit, "(A)") "# (in units of eV, folded to Floquet Brillouin zone)"

         WRITE (fname, "(2A)") TRIM(bs_env%floquet_dos_file), ".out"
         CALL open_file(TRIM(fname), unit_number=wunit, file_status="REPLACE", file_action="WRITE")
         WRITE (wunit, "(A)") "# Floquet Density of States: D(ω,k) = -1/π*Im[Tr_KS(G^R(ω,k))]"

         DO ikp_for_file = 1, n_kp_dos
            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, nkp_start + ikp_for_file)
            DO ispin = 1, n_spin

               ! Floquet bands for this k-point and spin
               WRITE (bunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
                  "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
               WRITE (bunit, "(A)") "#  Floquet band  Band-energy [eV]"
               DO i = 1, nao
                  WRITE (bunit, "(I8,F21.8)") i, all_floquet_bands(i, ispin, ikp_for_file)*evolt
               END DO

               ! Quasi-energies (folded to -ħΩ/2 < ε ≤ ħΩ/2) for this k-point and spin
               WRITE (qunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
                  "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
               WRITE (qunit, "(A)") "#  Floquet band  Quasi-energy [eV]"
               DO i = 1, nao
                  WRITE (qunit, "(I8,F21.8)") i, all_quasi_energies(i, ispin, ikp_for_file)*evolt
               END DO

               ! DOS A(k,ω) for this k-point and spin
               WRITE (wunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
                  "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
               WRITE (wunit, "(A)") "#Energy-E_F (eV)  A(ω,k) = DOS (1/eV)"
               DO i_E = 1, n_E
                  energy = E_min + i_E*energy_step
                  WRITE (wunit, "(2X,2G13.4)") energy*evolt, all_a_k(i_E, ispin, ikp_for_file)/evolt
               END DO

            END DO
         END DO

         CALL close_file(bunit)
         CALL close_file(qunit)
         CALL close_file(wunit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_floquet_results

END MODULE floquet_utils

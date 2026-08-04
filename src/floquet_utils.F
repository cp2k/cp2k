!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Helper routines for the Floquet-Bloch band-structure calculation (floquet_main).
!> \par History
!> \author Shridhar Shanbhag (27.01.2026)
! **************************************************************************************************
MODULE floquet_utils
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_submatrix,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE floquet_types,                   ONLY: floquet_env_type
   USE input_constants,                 ONLY: small_cell_full_kp
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
   USE machine,                         ONLY: m_memory_details
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: geeig_right,&
                                              gemm_square
   USE message_passing,                 ONLY: mp_comm_split_type_shared,&
                                              mp_comm_type,&
                                              mp_para_env_type
   USE physcon,                         ONLY: a_bohr,&
                                              evolt,&
                                              kelvin
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
   ! element) that the subgroup must be able to hold, distributed across its ranks.
   ! We size for 8 to leave headroom.
   INTEGER, PARAMETER, PRIVATE :: n_fm_work_copies = 8

   ! Public subroutines
   PUBLIC :: build_floquet_matrix, &
             compute_e_k_de_dk_dipole, &
             make_floquet_subgroups, &
             distribute_floquet_kp_data, &
             floquet_sector_weights, &
             check_floquet_convergence, &
             calculate_floquet_observables, &
             write_floquet_header, &
             write_floquet_results

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xkp ...
!> \param e_k ...
!> \param de_dk ...
!> \param do_parallel the option to distribute the results (e_k/de_dk) across MPI ranks.
!>        When .TRUE. k-point ikp is computed and stored only on rank MOD(ikp-1,num_pe)
!>        Default .FALSE. -> every rank computes/stores all k-points (replicated).
! **************************************************************************************************
   SUBROUTINE calculate_epsilon_derivative(qs_env, xkp, e_k, de_dk, do_parallel)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT), OPTIONAL                           :: e_k
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT), OPTIONAL    :: de_dk
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_parallel

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

      ! Distribution of k-points across ranks (mirrors qs_moment_kpoints_deep)
      ! do_parallel = .TRUE. => ikp stored in mepos==MOD(ikp-1,num_pe)
      ! do_parallel = .FALSE. => every rank computes all k-points (replicated)
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
!> \brief Momentum matrix elements p_nm(k) for one k-point and one spin.
!> \param e_k_kp_spin ...
!> \param de_dk_kp_spin ...
!> \param dipole_kp_spin ...
!> \param momentum ...
! **************************************************************************************************
   SUBROUTINE build_momentum_matrix(e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, momentum)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: e_k_kp_spin
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: de_dk_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: dipole_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: momentum

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
!> \brief Off-diagonal coupling block of H_F(k) for one k-point and spin, from precomputed data.
!> \param bs_env ...
!> \param e_k_kp_spin ...
!> \param de_dk_kp_spin ...
!> \param dipole_kp_spin ...
!> \param off_diag_m ...
! **************************************************************************************************
   SUBROUTINE build_off_diagonal_matrix(bs_env, e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, off_diag_m)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: e_k_kp_spin
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: de_dk_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: dipole_kp_spin
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: off_diag_m

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
         CPABORT("Invalid (too small) polarisation vector specified for POLARISATION")
      END IF

      amplitude = bs_env%floquet_amplitude
      e_vec(:) = amplitude*polarisation
      phi(:) = pi*bs_env%floquet_phi(:)
      omega = bs_env%floquet_omega

      ALLOCATE (momentum(3, nao, nao), source=z_zero)
      CALL build_momentum_matrix(e_k_kp_spin, de_dk_kp_spin, dipole_kp_spin, momentum)

      ! E_factor(α) = (i E(α) exp(i·φ(α)))/(2ω), where α = x,y,z
      DO i_dir = 1, 3
         efactor(i_dir) = gaussi*e_vec(i_dir)*CMPLX(COS(phi(i_dir)), SIN(phi(i_dir)), KIND=dp)/(2*omega)
      END DO

      ! off_diag_uv = Σ_α [p_uv^α(k) · E_factor(α)], where α = x,y,z
      off_diag_m(:, :) = z_zero
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
!> \brief Diagonal block of H_F(k) for Floquet sector f_index.
!> \param bs_env ...
!> \param e_k_kp_spin ...
!> \param f_index Floquet sector index
!> \param diag_e ...
! **************************************************************************************************
   SUBROUTINE build_diagonal_matrix(bs_env, e_k_kp_spin, f_index, diag_e)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: e_k_kp_spin
      INTEGER, INTENT(IN)                                :: f_index
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: diag_e

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_diagonal_matrix'

      INTEGER                                            :: handle, i, nao
      REAL(KIND=dp)                                      :: omega

      CALL timeset(routineN, handle)
      ! Builds the diagonal block of the H_F for Floquet sector index f_index
      ! Equilibrium band energies shifted by f_index·ħΩ, i.e., diag_e(n,n) = ε_{nk} + f_index·ħΩ.
      nao = SIZE(e_k_kp_spin)
      omega = bs_env%floquet_omega
      diag_e(:, :) = z_zero
      DO i = 1, nao
         diag_e(i, i) = e_k_kp_spin(i) + f_index*omega
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_diagonal_matrix

! **************************************************************************************************
!> \brief Builds the Floquet-Bloch Hamiltonian H_F(k) for a single k-point and spin channel from
!>        the band data currently held in floquet_env.
!> \param bs_env ...
!> \param floquet_env ...
!> \param floquet_matrix ...
! **************************************************************************************************
   SUBROUTINE build_floquet_matrix(bs_env, floquet_env, floquet_matrix)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(floquet_env_type), INTENT(IN)                 :: floquet_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: floquet_matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_floquet_matrix'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: conj_off_diag_m, diag_e, off_diag_m
      INTEGER                                            :: f_index, handle, i, i_f, max_f_index, &
                                                            n_fbands, nao

      CALL timeset(routineN, handle)

      nao = floquet_env%nao
      max_f_index = floquet_env%max_f_index
      n_fbands = 1 + 2*max_f_index

      ! Creates the floquet matrix H_F by placing the diagonal and off diagonal blocks
      ALLOCATE (diag_e(nao, nao), source=z_zero)
      ALLOCATE (off_diag_m(nao, nao), source=z_zero)
      ALLOCATE (conj_off_diag_m(nao, nao), source=z_zero)

      CALL build_off_diagonal_matrix(bs_env, floquet_env%e_k_kp_spin, floquet_env%de_dk_kp_spin, &
                                     floquet_env%dipole_kp_spin, off_diag_m)
      conj_off_diag_m(:, :) = CONJG(TRANSPOSE(off_diag_m(:, :)))

      floquet_matrix%local_data(:, :) = z_zero
      DO i = 1, n_fbands
         i_f = 1 + (i - 1)*nao
         f_index = i - max_f_index - 1
         CALL build_diagonal_matrix(bs_env, floquet_env%e_k_kp_spin, f_index, diag_e)
         CALL cp_cfm_set_submatrix(floquet_matrix, diag_e, i_f, i_f)
         IF (i > 1) THEN
            CALL cp_cfm_set_submatrix(floquet_matrix, off_diag_m, i_f, i_f - nao)
            CALL cp_cfm_set_submatrix(floquet_matrix, conj_off_diag_m, i_f - nao, i_f)
         END IF
      END DO
      DEALLOCATE (diag_e, off_diag_m, conj_off_diag_m)

      CALL timestop(handle)
   END SUBROUTINE build_floquet_matrix

! **************************************************************************************************
!> \brief Make this k-point/spin's band data available to every rank of the owning subgroup.
!> \param para_env the global parallel environment
!> \param para_env_sub the subgroup parallel environment
!> \param ispin the spin channel
!> \param ikp the DOS k-point index
!> \param e_k band energies for all DOS k-points, distributed
!> \param de_dk band-energy k-derivatives, distributed
!> \param dipole dipole matrix elements, distributed
!> \param floquet_env ...
! **************************************************************************************************
   SUBROUTINE distribute_floquet_kp_data(para_env, para_env_sub, ispin, ikp, e_k, de_dk, dipole, &
                                         floquet_env)
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: ispin, ikp
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: e_k
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: de_dk
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :, :), &
         INTENT(IN)                                      :: dipole
      TYPE(floquet_env_type), INTENT(INOUT)              :: floquet_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'distribute_floquet_kp_data'

      INTEGER                                            :: handle, loc_idx, local_src, owner

      CALL timeset(routineN, handle)

      ! Identify the owner's rank within the subgroup (all subgroup ranks agree on local_src).
      owner = MOD(ikp - 1, para_env%num_pe)

      local_src = -1
      IF (para_env%mepos == owner) local_src = para_env_sub%mepos
      CALL para_env_sub%max(local_src)

      ! The owner copies its single-spin slice of this k-point's distributed band data into the
      ! broadcast buffers, then broadcasts them to the whole subgroup.
      IF (para_env%mepos == owner) THEN
         loc_idx = CEILING(REAL(ikp)/para_env%num_pe)
         floquet_env%e_k_kp_spin(:) = e_k(ispin, loc_idx, :)
         floquet_env%de_dk_kp_spin(:, :) = de_dk(ispin, loc_idx, :, :)
         floquet_env%dipole_kp_spin(:, :, :) = dipole(ispin, loc_idx, :, :, :)
      END IF
      CALL para_env_sub%bcast(floquet_env%e_k_kp_spin, local_src)
      CALL para_env_sub%bcast(floquet_env%de_dk_kp_spin, local_src)
      CALL para_env_sub%bcast(floquet_env%dipole_kp_spin, local_src)

      CALL timestop(handle)

   END SUBROUTINE distribute_floquet_kp_data

! **************************************************************************************************
!> \brief Precompute, on the global communicator, the band quantities needed to assemble the
!>        Floquet-Bloch Hamiltonian for all DOS k-points (all spins), distributed across ranks.
!> \param qs_env ...
!> \param bs_env ...
!> \param e_k ...
!> \param de_dk ...
!> \param dipole ...
! **************************************************************************************************
   SUBROUTINE compute_e_k_de_dk_dipole(qs_env, bs_env, e_k, de_dk, dipole)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: e_k
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: de_dk
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :), INTENT(OUT)           :: dipole

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_e_k_de_dk_dipole'

      INTEGER                                            :: handle, ikp, nkp_only_bs, nkp_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp_all

      CALL timeset(routineN, handle)

      ! The band-structure k-points are the last nkp_only_bs entries of kpoints_DOS%xkp, stored
      ! after the nkp_only_DOS DOS-only k-points.
      nkp_start = bs_env%nkp_only_DOS
      nkp_only_bs = bs_env%nkp_only_bs

      ! Collect the coordinates of all band-structure k-points (global index nkp_start + local).
      ALLOCATE (xkp_all(3, nkp_only_bs))
      DO ikp = 1, nkp_only_bs
         xkp_all(:, ikp) = bs_env%kpoints_DOS%xkp(1:3, nkp_start + ikp)
      END DO

      ! Calculate and distribute the results across ranks
      ! One k-point per rank, round-robin by global rank
      ! So the results for ikp are stored in mepos==MOD(ikp-1,num_pe)
      CALL calculate_epsilon_derivative(qs_env, xkp_all, e_k=e_k, de_dk=de_dk, do_parallel=.TRUE.)
      CALL qs_moment_kpoints_deep(qs_env, xkp_all, dipole, do_parallel=.TRUE.)

      DEALLOCATE (xkp_all)

      CALL timestop(handle)
   END SUBROUTINE compute_e_k_de_dk_dipole

! **************************************************************************************************
!> \brief Finds the number of MPI ranks that share a physical node, determined by splitting the
!>        given communicator into node-local communicators via MPI's shared-memory split type.
!> \param para_env ...
!> \return the number of ranks on the calling rank's node
! **************************************************************************************************
   FUNCTION find_ranks_per_node(para_env) RESULT(ranks_per_node)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: ranks_per_node

      TYPE(mp_comm_type)                                 :: node_comm

      CALL node_comm%from_split_type(para_env, mp_comm_split_type_shared, key=para_env%mepos)
      ranks_per_node = MAX(1, node_comm%num_pe)
      CALL node_comm%free()

   END FUNCTION find_ranks_per_node

! **************************************************************************************************
!> \brief Choose the optimal number of MPI ranks per subgroup for Floquet calculations.
!> \param qs_env ...
!> \param bs_env ...
!> \return the chosen number of ranks per subgroup (1 .. num_pe)
! **************************************************************************************************
   FUNCTION floquet_determine_subgroup_size(qs_env, bs_env) RESULT(group_size)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: group_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'floquet_determine_subgroup_size'

      INTEGER                                            :: g_load, g_mem, handle, n_f_size, nao, &
                                                            nkp_only_bs, ranks_per_node, unit_nr
      INTEGER(KIND=int_8)                                :: Buffers, Cached, MemFree, MemLikelyFree, &
                                                            MemTotal, needed_bytes, Slab, &
                                                            SReclaimable, usable_per_rank
      REAL(KIND=dp)                                      :: mem_fill_fraction
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      unit_nr = bs_env%unit_nr
      n_f_size = nao*(1 + 2*bs_env%max_floquet_index)
      nkp_only_bs = bs_env%nkp_only_bs
      mem_fill_fraction = bs_env%floquet_mem_fill_fraction

      ! We choose the subgroup size to be the larger of two floors:
      ! 1. Memory floor G_mem: the subgroup must jointly hold the Floquet working set,
      !    G_mem = ceil(needed_bytes/usable per-rank memory), where the usable memory is
      !    the mem_fill_fraction of the per-rank FREE memory.
      ! 2. K-point floor G_load = floor(num_pe / nkp_only_bs): the LARGEST G that still
      !    yields at least nkp_only_bs subgroups (floor(num_pe/G) >= nkp_only_bs). When ranks
      !    outnumber k-points this hands every k-point its own (larger, hence faster)
      !    subgroup instead of leaving spare ranks idle.

      ! usable_per_rank = a fraction of the per-rank FREE memory, not including SCF and other data
      CALL m_memory_details(MemTotal, MemFree, Buffers, Cached, Slab, SReclaimable, MemLikelyFree)
      ranks_per_node = find_ranks_per_node(para_env)
      usable_per_rank = INT(mem_fill_fraction*REAL(MemFree, dp), int_8)/INT(ranks_per_node, int_8)

      ! Total memory the subgroup must hold: n_fm_work_copies copies of the Floquet
      ! Hamiltonian, each n_f_size^2 complex(dp) numbers at 16 bytes each.
      needed_bytes = 16_int_8*INT(n_fm_work_copies, int_8)*INT(n_f_size, int_8)**2

      IF (usable_per_rank <= 0_int_8) THEN
         ! Memory info unavailable (e.g. no /proc found) use one subgroup spanning all ranks
         g_mem = para_env%num_pe
      ELSE
         ! Memory floor: smallest G with needed_bytes/G <= usable_per_rank
         g_mem = INT(MIN((needed_bytes + usable_per_rank - 1_int_8)/usable_per_rank, &
                         INT(para_env%num_pe, int_8)))
      END IF

      IF (g_mem > para_env%num_pe) CPWARN("Total Memory Likely Insufficient, process may be killed")

      ! K-point floor: floor(num_pe / nkp_only_bs) = largest G that gives at least nkp_only_bs subgroups
      g_load = para_env%num_pe/MAX(nkp_only_bs, 1)

      group_size = MAX(g_mem, g_load)
      group_size = MIN(group_size, para_env%num_pe)
      group_size = MAX(group_size, 1)
      CALL para_env%max(group_size)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') ""
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Detected MPI ranks per node:", ranks_per_node
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Free memory per rank:", &
            (MemLikelyFree/INT(ranks_per_node, int_8))/1048576_int_8, " MB"
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Usable per rank (reserve applied):", &
            usable_per_rank/1048576_int_8, " MB"
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Floquet copies needed:", n_fm_work_copies
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Total Floquet working set:", &
            needed_bytes/1048576_int_8, " MB"
         WRITE (unit_nr, '(T2,A,T66,I12,A)') "FLOQUET MEMORY | Floquet working set per rank:", &
            (needed_bytes/INT(MAX(group_size, 1), int_8))/1048576_int_8, " MB"
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Memory floor (ranks/subgroup):", g_mem
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | K-point floor (ranks/subgroup):", g_load
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Chosen MPI ranks per subgroup:", group_size
      END IF

      CALL timestop(handle)
   END FUNCTION floquet_determine_subgroup_size

! **************************************************************************************************
!> \brief Split the para_env into subgroups so that the Floquet Hamiltonian is distributed across
!>        the ranks of each subgroup. Also creates the subgroup BLACS context.
!> \param qs_env ...
!> \param bs_env ...
!> \param para_env_sub the created subgroup parallel environment
!> \param blacs_env_sub the created subgroup BLACS context
!> \param group_distribution subgroup index of every global rank, 0:num_pe-1
!> \param ngroups the number of subgroups created
! **************************************************************************************************
   SUBROUTINE make_floquet_subgroups(qs_env, bs_env, para_env_sub, blacs_env_sub, &
                                     group_distribution, ngroups)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: group_distribution
      INTEGER, INTENT(OUT)                               :: ngroups

      CHARACTER(LEN=*), PARAMETER :: routineN = 'make_floquet_subgroups'

      INTEGER                                            :: group_size, handle, n_spin, nkp_only_bs, &
                                                            stride_kp, unit_nr
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env)
      unit_nr = bs_env%unit_nr
      nkp_only_bs = bs_env%nkp_only_bs
      n_spin = bs_env%n_spin

      ! We split the global communicator into subgroups with the Floquet Hamiltonian distributed
      ! across the processes in each subgroup and the k-points distributed across subgroups.
      group_size = floquet_determine_subgroup_size(qs_env, bs_env)

      IF (nkp_only_bs < para_env%num_pe) THEN
         stride_kp = para_env%num_pe/group_size
      ELSE
         stride_kp = 1
      END IF
      ALLOCATE (group_distribution(0:para_env%num_pe - 1))
      ALLOCATE (para_env_sub)
      CALL para_env_sub%from_split(comm=para_env, ngroups=ngroups, &
                                   group_distribution=group_distribution, &
                                   subgroup_min_size=group_size, stride=stride_kp)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Subgroup rank stride:", stride_kp
         WRITE (unit_nr, '(T2,A,T66,I15)') "FLOQUET MEMORY | Number of subgroups:", ngroups
         WRITE (unit_nr, '(/,T2,A,I5,A,I1,A)') "FLOQUET CALCULATIONS PROGRESS OUT OF", nkp_only_bs, &
            " K-POINTS AND ", n_spin, " SPINS"
      END IF

      NULLIFY (blacs_env_sub)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      CALL timestop(handle)

   END SUBROUTINE make_floquet_subgroups

! **************************************************************************************************
!> \brief Central and boundary-sector weights of every Floquet eigenvector, stored into
!>        floquet_env%w0 and floquet_env%wE.
!> \param floquet_env ...
!> \param cfm_eigenvectors Floquet eigenvectors
! **************************************************************************************************
   SUBROUTINE floquet_sector_weights(floquet_env, cfm_eigenvectors)
      TYPE(floquet_env_type), INTENT(INOUT)              :: floquet_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_eigenvectors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'floquet_sector_weights'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: v_block
      INTEGER                                            :: handle, max_f_index, n_f_size, nao

      CALL timeset(routineN, handle)

      nao = floquet_env%nao
      max_f_index = floquet_env%max_f_index
      n_f_size = floquet_env%n_f_size

      ALLOCATE (v_block(nao, n_f_size))

      ! Evaluate w0_α = Σ_{n=1..nao} |<n,m=0|α>|^2 (central sector)
      CALL cp_cfm_get_submatrix(cfm_eigenvectors, v_block, 1 + nao*max_f_index, 1)
      floquet_env%w0(:) = SUM(ABS(v_block)**2, DIM=1)

      ! Evaluate wE_α = Σ_{n=1..nao} (|<n,m=-M|α>|^2 + |<n,m=+M|α>|^2) (outermost rungs)
      floquet_env%wE(:) = 0.0_dp
      IF (max_f_index >= 1) THEN
         CALL cp_cfm_get_submatrix(cfm_eigenvectors, v_block, 1, 1)
         floquet_env%wE(:) = SUM(ABS(v_block)**2, DIM=1)

         CALL cp_cfm_get_submatrix(cfm_eigenvectors, v_block, 1 + nao*2*max_f_index, 1)
         floquet_env%wE(:) = floquet_env%wE(:) + SUM(ABS(v_block)**2, DIM=1)
      END IF

      DEALLOCATE (v_block)

      ! The orthonormal rows of a unitary matrix are such that Σ_α w0_α = nao.
      CPASSERT(ABS(SUM(floquet_env%w0) - REAL(nao, dp)) < 1.0E-6_dp*REAL(nao, dp))

      CALL timestop(handle)

   END SUBROUTINE floquet_sector_weights

! **************************************************************************************************
!> \brief Checks that MAX_FLOQUET_INDEX is large enough and the Floquet Hamiltonian was
!>        truncated far enough away from the central sector to leave it unperturbed.
!> \param bs_env ...
!> \param floquet_env holds the central- (w0) and outermost-sector (wE) weights
! **************************************************************************************************
   SUBROUTINE check_floquet_convergence(bs_env, floquet_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(floquet_env_type), INTENT(IN)                 :: floquet_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_floquet_convergence'

      CHARACTER(LEN=default_string_length)               :: msg
      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: leak

      ! The truncation is converged if m=0 states have small weight on outermost rungs m = ±M
      CALL timeset(routineN, handle)

      ! MAX_FLOQUET_INDEX = 0 means no sidebands at all
      IF (bs_env%max_floquet_index < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! No cleanly m=0-dominated state exists: the drive has hybridised every band with
      ! its sidebands, strong field but not a truncation failure.
      IF (.NOT. ANY(floquet_env%w0 > 0.5_dp)) THEN
         CPWARN("Floquet: no m=0-dominated state; cannot assess MAX_FLOQUET_INDEX convergence.")
         CALL timestop(handle)
         RETURN
      END IF

      leak = MAXVAL(floquet_env%wE, MASK=(floquet_env%w0 > 0.5_dp))

      IF (bs_env%eps_floquet > 0.0_dp .AND. leak > bs_env%eps_floquet) THEN
         WRITE (msg, '(A,ES10.2E2,A,ES10.2E2)') &
            "MAX_FLOQUET_INDEX is too small. Leak: ", leak, "exceeds EPS_FLOQUET: ", bs_env%eps_floquet
         CPABORT(TRIM(msg))
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_floquet_convergence

! **************************************************************************************************
!> \brief Energy origin for the Floquet calculations, chosen to be the valence band maximum.
!> \param bs_env ...
!> \return the valence band maximum
! **************************************************************************************************
   FUNCTION floquet_reference_energy(bs_env) RESULT(mu)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: mu

      IF (bs_env%do_gw .OR. &
          bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp) THEN
         mu = bs_env%band_edges_scf%VBM
      ELSE IF (bs_env%n_spin == 1) THEN
         mu = bs_env%eigenval_scf_Gamma(bs_env%n_occ(1), 1)
      ELSE
         mu = MAX(bs_env%eigenval_scf_Gamma(bs_env%n_occ(1), 1), &
                  bs_env%eigenval_scf_Gamma(bs_env%n_occ(2), 2))
      END IF

   END FUNCTION floquet_reference_energy

! **************************************************************************************************
!> \brief Calculate all Floquet observables for one k-point and spin: the DOS (a_k), the m=0 bands,
!>        the quasi-energies, and store in floquet_env.
!> \param bs_env ...
!> \param para_env_sub ...
!> \param ispin ...
!> \param ikp ...
!> \param floquet_env holds eigenvalues, w0 (input) and a_k, m0_energies, m0_weights,
!>        quasi_energies (output)
! **************************************************************************************************
   SUBROUTINE calculate_floquet_observables(bs_env, para_env_sub, ispin, ikp, floquet_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, INTENT(IN)                                :: ispin, ikp
      TYPE(floquet_env_type), INTENT(INOUT)              :: floquet_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_floquet_observables'

      INTEGER                                            :: handle, i, i_E, j, n_E, n_f_size, nao
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: work
      REAL(KIND=dp)                                      :: broad, cum, E_min, energy, energy_step, &
                                                            mu, omega, target_level

      CALL timeset(routineN, handle)

      n_f_size = floquet_env%n_f_size
      nao = floquet_env%nao
      n_E = floquet_env%n_E

      mu = floquet_reference_energy(bs_env)
      omega = bs_env%floquet_omega
      broad = bs_env%broadening_floquet
      energy_step = bs_env%energy_step_floquet
      E_min = mu - bs_env%energy_window_floquet

      ! (1) Spectral function / DOS, A(k,E) = -(1/π) Im Tr_KS G^R_00(k,E),
      ! where G^R_00(k,E) = (E + iη - H_F)^(-1) is the retarded Green's function
      ! With H_F = Σ_α λ_α |α><α|, we simplify
      ! A(k,E) = -(1/π) Σ_α w0_α Im[1/(E + iη - λ_α)], η = broadening/2.
      ! This avoids the per-energy dense inversion of the full H_F.

      DO i_E = 1, n_E
         energy = E_min + i_E*energy_step
         floquet_env%a_k(i_E) = -SUM(floquet_env%w0(:)* &
                                     AIMAG(z_one/(energy + gaussi*broad/2.0_dp - &
                                                  floquet_env%eigenvalues(:))))/pi
      END DO

      ! (2) m=0 band structure by weighted count
      cum = 0.0_dp
      j = 0
      DO i = 1, nao
         target_level = REAL(i, dp) - 0.5_dp
         DO WHILE (cum < target_level .AND. j < n_f_size)
            j = j + 1
            cum = cum + floquet_env%w0(j)
         END DO
         floquet_env%m0_energies(i) = floquet_env%eigenvalues(j) - mu
         floquet_env%m0_weights(i) = floquet_env%w0(j)
      END DO

      ! (3) Fold the m=0 bands (already relative to the VBM) into the first Floquet Brillouin zone.
      floquet_env%quasi_energies(:) = floquet_env%m0_energies &
                                      - omega*REAL(CEILING(floquet_env%m0_energies/omega - 0.5_dp), dp)

      ! Sort the quasi-energies
      ALLOCATE (work(nao))
      CALL sort(floquet_env%quasi_energies, nao, work)
      DEALLOCATE (work)

      ! Store this result on the subgroup source only, so the global sum below picks up
      ! each work item only once
      IF (para_env_sub%is_source()) THEN
         floquet_env%all_quasi_energies(:, ispin, ikp) = floquet_env%quasi_energies(:)
         floquet_env%all_a_k(:, ispin, ikp) = floquet_env%a_k(:)
         floquet_env%all_m0_energies(:, ispin, ikp) = floquet_env%m0_energies(:)
         floquet_env%all_m0_weights(:, ispin, ikp) = floquet_env%m0_weights(:)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_floquet_observables

! **************************************************************************************************
!> \brief Print the Floquet header with the input parameters and a short description of the output.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE write_floquet_header(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_floquet_header'

      INTEGER                                            :: handle, n_f_size, unit_nr

      CALL timeset(routineN, handle)

      unit_nr = bs_env%unit_nr
      n_f_size = bs_env%n_ao*(1 + 2*bs_env%max_floquet_index)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(T2,A)') ' '
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A,A51,A27)') '-', 'FLOQUET BANDSTRUCTURE CALCULATION', '-'
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A)') ' '

         WRITE (unit_nr, '(T2,A,T37,A,T67,ES12.2E2)') "FLOQUET PARAMETERS", "Amplitude [V/m]:", &
            bs_env%floquet_amplitude*evolt/a_bohr
         WRITE (unit_nr, '(T37,A,T67,F12.4)') "Frequency [eV]:", bs_env%floquet_omega*evolt
         WRITE (unit_nr, '(T37,A)') REPEAT("-", 44)
         WRITE (unit_nr, '(T37,A,T55,3F8.4)') "Polarisation:", bs_env%floquet_polarisation(1:3)
         WRITE (unit_nr, '(T37,A,T55,3F8.4)') "Phase offsets:", pi*bs_env%floquet_phi(1:3)
         WRITE (unit_nr, '(T37,A)') REPEAT("-", 44)
         WRITE (unit_nr, '(T37,A,T67,I12)') "Max Floquet index:", bs_env%max_floquet_index
         WRITE (unit_nr, '(T37,A,T67,I12)') "Floquet Hamiltonian Size:", n_f_size
         WRITE (unit_nr, '(T37,A)') REPEAT("-", 44)
         WRITE (unit_nr, '(T37,A,T67,F12.4)') &
            "Energy window [eV]:", bs_env%energy_window_floquet*evolt
         WRITE (unit_nr, '(T37,A,T67,F12.4)') "Energy step [eV]:", bs_env%energy_step_floquet*evolt
         WRITE (unit_nr, '(T37,A,T67,F12.4)') "Broadening [eV]:", bs_env%broadening_floquet*evolt
         WRITE (unit_nr, '(T37,A)') REPEAT("-", 44)
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

      CALL timestop(handle)

   END SUBROUTINE write_floquet_header

! **************************************************************************************************
!> \brief Sum the accumulated results across MPI ranks and write the m=0 band structure,
!>        the quasi-energies and the Floquet DOS to their files.
!> \param bs_env ...
!> \param floquet_env ...
! **************************************************************************************************
   SUBROUTINE write_floquet_results(bs_env, floquet_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(floquet_env_type), INTENT(INOUT)              :: floquet_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_floquet_results'

      CHARACTER(LEN=default_string_length)               :: fname
      INTEGER                                            :: bunit, handle, i, i_E, ikp_for_file, &
                                                            ispin, n_E, n_spin, nao, nkp_only_bs, &
                                                            nkp_start, qunit, wunit
      REAL(KIND=dp)                                      :: E_min, energy, energy_step, f_occ, kT, &
                                                            mu, x
      REAL(KIND=dp), DIMENSION(3)                        :: xkp

      CALL timeset(routineN, handle)

      ! Collect all results on the global communicator and write them
      CALL bs_env%para_env%sum(floquet_env%all_quasi_energies)
      CALL bs_env%para_env%sum(floquet_env%all_a_k)
      CALL bs_env%para_env%sum(floquet_env%all_m0_energies)
      CALL bs_env%para_env%sum(floquet_env%all_m0_weights)

      IF (bs_env%para_env%is_source()) THEN

         nao = floquet_env%nao
         n_spin = floquet_env%n_spin
         nkp_only_bs = floquet_env%nkp_only_bs
         n_E = floquet_env%n_E
         nkp_start = bs_env%nkp_only_DOS

         ! All three share one energy origin, the VBM.
         mu = floquet_reference_energy(bs_env)
         energy_step = bs_env%energy_step_floquet
         E_min = mu - bs_env%energy_window_floquet

         ! k_B T in Hartree (only used when TEMPERATURE > 0); E[K] = E[Hartree]*kelvin.
         kT = bs_env%floquet_temperature/kelvin

         ! Each file is opened once (REPLACE) and written for every band-structure k-point and spin.

         ! The m=0 band structure
         WRITE (fname, "(2A)") TRIM(bs_env%floquet_bs_file), ".bs"
         CALL open_file(TRIM(fname), unit_number=bunit, file_status="REPLACE", file_action="WRITE")
         WRITE (bunit, "(A)") "# Floquet m=0 (central-sector) band structure"
         WRITE (bunit, "(A)") "# (in units of eV, relative to the VBM, not folded)"
         WRITE (bunit, "(A)") "# w = sum_n |<n,m=0|alpha>|^2 in [0,1]: w~1 clean m=0 replica,"
         WRITE (bunit, "(A)") "# w~0.5 hybridised with a sideband (drive near resonance)"

         ! Quasi-energies
         WRITE (fname, "(2A)") TRIM(bs_env%floquet_qe_file), ".bs"
         CALL open_file(TRIM(fname), unit_number=qunit, file_status="REPLACE", file_action="WRITE")
         WRITE (qunit, "(A)") "# Quasi-energies obtained by diagonalising the Floquet Hamiltonian"
         WRITE (qunit, "(A)") "# (in units of eV, the m=0 bands relative to the VBM, folded to"
         WRITE (qunit, "(A)") "# the first Floquet Brillouin zone -hbar*Omega/2 < e <= hbar*Omega/2)"

         ! DOS
         WRITE (fname, "(2A)") TRIM(bs_env%floquet_dos_file), ".out"
         CALL open_file(TRIM(fname), unit_number=wunit, file_status="REPLACE", file_action="WRITE")
         WRITE (wunit, "(A)") "# Floquet Density of States: D(ω,k) = -1/π*Im[Tr_KS(G^R(ω,k))]"

         DO ikp_for_file = 1, nkp_only_bs
            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, nkp_start + ikp_for_file)
            DO ispin = 1, n_spin

               ! <floquet_bs_file>.bs the m=0 band structure (absolute energies) and its weights
               WRITE (bunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
                  "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
               WRITE (bunit, "(A)") "#  Floquet band  Energy [eV]         m=0 weight"
               DO i = 1, nao
                  WRITE (bunit, "(I8,F21.8,F17.5)") i, &
                     floquet_env%all_m0_energies(i, ispin, ikp_for_file)*evolt, &
                     floquet_env%all_m0_weights(i, ispin, ikp_for_file)
               END DO

               ! <floquet_qe_file>.bs Quasi-energies (folded to -ħΩ/2 < ε ≤ ħΩ/2)
               WRITE (qunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
                  "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
               WRITE (qunit, "(A)") "#  Floquet band  Quasi-energy [eV]"
               DO i = 1, nao
                  WRITE (qunit, "(I8,F21.8)") i, &
                     floquet_env%all_quasi_energies(i, ispin, ikp_for_file)*evolt
               END DO

               ! <floquet_dos_file>.out the k-resolved DOS A(k,ω), and (if TEMPERATURE > 0) the
               ! occupied spectral weight f(E)*A(k,ω), f the Fermi-Dirac occupation at the VBM.
               WRITE (wunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
                  "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
               IF (bs_env%floquet_temperature > 0.0_dp) THEN
                  WRITE (wunit, "(A)") "#Energy-VBM (eV)  A(ω,k) = DOS (1/eV)   f*A = occupied DOS (1/eV)"
                  DO i_E = 1, n_E
                     energy = E_min + i_E*energy_step
                     x = (energy - mu)/kT                    ! (E - VBM)/k_B T, dimensionless
                     IF (x > 40.0_dp) THEN                   ! guard EXP overflow at low T
                        f_occ = 0.0_dp
                     ELSE IF (x < -40.0_dp) THEN
                        f_occ = 1.0_dp
                     ELSE
                        f_occ = 1.0_dp/(EXP(x) + 1.0_dp)
                     END IF
                     WRITE (wunit, "(2X,3G13.4)") (energy - mu)*evolt, &
                        floquet_env%all_a_k(i_E, ispin, ikp_for_file)/evolt, &
                        f_occ*floquet_env%all_a_k(i_E, ispin, ikp_for_file)/evolt
                  END DO
               ELSE
                  WRITE (wunit, "(A)") "#Energy-VBM (eV)  A(ω,k) = DOS (1/eV)"
                  DO i_E = 1, n_E
                     energy = E_min + i_E*energy_step
                     WRITE (wunit, "(2X,2G13.4)") (energy - mu)*evolt, &
                        floquet_env%all_a_k(i_E, ispin, ikp_for_file)/evolt
                  END DO
               END IF

            END DO
         END DO

         CALL close_file(bunit)
         CALL close_file(qunit)
         CALL close_file(wunit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_floquet_results

END MODULE floquet_utils

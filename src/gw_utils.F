!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 07.2023
! **************************************************************************************************
MODULE gw_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Graml2024,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig_canon
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_p_type, dbcsr_set, &
        dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_create, dbt_destroy, dbt_filter, dbt_iterator_blocks_left, &
        dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, dbt_iterator_type, &
        dbt_mp_environ_pgrid, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE gw_integrals,                    ONLY: build_3c_integral_block
   USE input_constants,                 ONLY: do_potential_truncated,&
                                              large_cell_Gamma,&
                                              ri_rpa_g0w0_crossing_newton,&
                                              small_cell_full_kp,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_type
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup,&
                                              cp_libint_static_init
   USE machine,                         ONLY: m_memory,&
                                              m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: gcd
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_para_env_type
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
   USE minimax_exp_gw,                  ONLY: get_exp_minimax_coeff_gw
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff,&
                                              get_rpa_minimax_coeff_larger_grid
   USE mp2_gpw,                         ONLY: create_mat_munu
   USE mp2_grids,                       ONLY: get_l_sq_wghts_cos_tf_t_to_w,&
                                              get_l_sq_wghts_cos_tf_w_to_t,&
                                              get_l_sq_wghts_sin_tf_t_to_w
   USE mp2_ri_2c,                       ONLY: trunc_coulomb_for_exchange
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE post_scf_bandstructure_types,    ONLY: band_edges_type,&
                                              post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: get_fname
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_part_release,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tensors,                      ONLY: build_3c_integrals,&
                                              build_3c_neighbor_lists,&
                                              get_tensor_occupancy,&
                                              neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_type
   USE rpa_gw,                          ONLY: continuation_pade
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: create_and_init_bs_env_for_gw, de_init_bs_env, get_i_j_atoms, get_VBM_CBM_bandgaps, &
             kpoint_init_cell_index_simple, compute_xkp, time_to_freq, analyt_conti_and_print, &
             add_R, is_cell_in_index_to_cell

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param bs_sec ...
! **************************************************************************************************
   SUBROUTINE create_and_init_bs_env_for_gw(qs_env, bs_env, bs_sec)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: bs_sec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_and_init_bs_env_for_gw'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cite_reference(Graml2024)

      CALL read_gw_input_parameters(bs_env, bs_sec)

      CALL print_header_and_input_parameters(bs_env)

      CALL setup_AO_and_RI_basis_set(qs_env, bs_env)

      CALL get_RI_basis_and_basis_function_indices(qs_env, bs_env)

      CALL set_heuristic_parameters(bs_env, qs_env)

      CALL cp_libint_static_init()

      CALL setup_kpoints_chi_eps_W(bs_env, bs_env%kpoints_chi_eps_W)

      IF (bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp) THEN
         CALL setup_cells_3c(qs_env, bs_env)
      END IF

      CALL set_parallelization_parameters(qs_env, bs_env)

      CALL allocate_matrices(qs_env, bs_env)

      CALL compute_V_xc(qs_env, bs_env)

      CALL create_tensors(qs_env, bs_env)

      SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
      CASE (large_cell_Gamma)

         CALL allocate_GW_eigenvalues(bs_env)

         CALL check_sparsity_3c(qs_env, bs_env)

         CALL set_sparsity_parallelization_parameters(bs_env)

         CALL check_for_restart_files(qs_env, bs_env)

      CASE (small_cell_full_kp)

         CALL compute_3c_integrals(qs_env, bs_env)

         CALL setup_cells_Delta_R(bs_env)

         CALL setup_parallelization_Delta_R(bs_env)

         CALL compute_cfm_mo_coeff_kp_and_eigenval_scf_kp(qs_env, bs_env)

         CALL allocate_matrices_small_cell_full_kp(qs_env, bs_env)

         CALL trafo_V_xc_R_to_kp(qs_env, bs_env)

      END SELECT

      CALL setup_time_and_frequency_minimax_grid(bs_env)

      ! free memory in qs_env; only if one is not calculating the LDOS because
      ! we need real-space grid operations in pw_env, task_list for the LDOS
      ! Recommendation in case of memory issues: first perform GW calculation without calculating
      !                                          LDOS (to safe memor). Then, use GW restart files
      !                                          in a subsequent calculation to calculate the LDOS
      IF (.NOT. bs_env%do_ldos) THEN
         CALL qs_env_part_release(qs_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE create_and_init_bs_env_for_gw

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE de_init_bs_env(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'de_init_bs_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      ! deallocate quantities here which:
      ! 1. cannot be deallocated in bs_env_release due to circular dependencies
      ! 2. consume a lot of memory and should not be kept until the quantity is
      !    deallocated in bs_env_release

      IF (ASSOCIATED(bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(bs_env%nl_3c)

      CALL cp_libint_static_cleanup()

      CALL timestop(handle)

   END SUBROUTINE de_init_bs_env

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param bs_sec ...
! **************************************************************************************************
   SUBROUTINE read_gw_input_parameters(bs_env, bs_sec)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: bs_sec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_gw_input_parameters'

      INTEGER                                            :: handle
      TYPE(section_vals_type), POINTER                   :: gw_sec

      CALL timeset(routineN, handle)

      NULLIFY (gw_sec)
      gw_sec => section_vals_get_subs_vals(bs_sec, "GW")

      CALL section_vals_val_get(gw_sec, "NUM_TIME_FREQ_POINTS", i_val=bs_env%num_time_freq_points)
      CALL section_vals_val_get(gw_sec, "EPS_FILTER", r_val=bs_env%eps_filter)
      CALL section_vals_val_get(gw_sec, "MEMORY_PER_PROC", r_val=bs_env%input_memory_per_proc_GB)
      CALL section_vals_val_get(gw_sec, "APPROX_KP_EXTRAPOL", l_val=bs_env%approx_kp_extrapol)

      CALL timestop(handle)

   END SUBROUTINE read_gw_input_parameters

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_AO_and_RI_basis_set(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_AO_and_RI_basis_set'

      INTEGER                                            :: handle, natom, nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      natom=natom, nkind=nkind)

      ! set up basis
      ALLOCATE (bs_env%sizes_RI(natom), bs_env%sizes_AO(natom))
      ALLOCATE (bs_env%basis_set_RI(nkind), bs_env%basis_set_AO(nkind))

      CALL basis_set_list_setup(bs_env%basis_set_RI, "RI_AUX", qs_kind_set)
      CALL basis_set_list_setup(bs_env%basis_set_AO, "ORB", qs_kind_set)

      CALL get_particle_set(particle_set, qs_kind_set, nsgf=bs_env%sizes_RI, &
                            basis=bs_env%basis_set_RI)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=bs_env%sizes_AO, &
                            basis=bs_env%basis_set_AO)

      CALL timestop(handle)

   END SUBROUTINE setup_AO_and_RI_basis_set

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_RI_basis_and_basis_function_indices(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_RI_basis_and_basis_function_indices'

      INTEGER                                            :: handle, i_RI, iatom, ikind, iset, &
                                                            max_AO_bf_per_atom, n_ao_test, n_atom, &
                                                            n_kind, n_RI, nset, nsgf, u
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(:), POINTER                     :: l_max, l_min, nsgf_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! determine RI basis set size
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)

      n_kind = SIZE(qs_kind_set)
      n_atom = bs_env%n_atom

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO ikind = 1, n_kind
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, &
                          basis_type="RI_AUX")
         CPASSERT(ASSOCIATED(basis_set_a))
      END DO

      ALLOCATE (bs_env%i_RI_start_from_atom(n_atom))
      ALLOCATE (bs_env%i_RI_end_from_atom(n_atom))
      ALLOCATE (bs_env%i_ao_start_from_atom(n_atom))
      ALLOCATE (bs_env%i_ao_end_from_atom(n_atom))

      n_RI = 0
      DO iatom = 1, n_atom
         bs_env%i_RI_start_from_atom(iatom) = n_RI + 1
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="RI_AUX")
         n_RI = n_RI + nsgf
         bs_env%i_RI_end_from_atom(iatom) = n_RI
      END DO
      bs_env%n_RI = n_RI

      max_AO_bf_per_atom = 0
      n_ao_test = 0
      DO iatom = 1, n_atom
         bs_env%i_ao_start_from_atom(iatom) = n_ao_test + 1
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="ORB")
         n_ao_test = n_ao_test + nsgf
         bs_env%i_ao_end_from_atom(iatom) = n_ao_test
         max_AO_bf_per_atom = MAX(max_AO_bf_per_atom, nsgf)
      END DO
      CPASSERT(n_ao_test == bs_env%n_ao)
      bs_env%max_AO_bf_per_atom = max_AO_bf_per_atom

      ALLOCATE (bs_env%l_RI(n_RI))
      i_RI = 0
      DO iatom = 1, n_atom
         ikind = kind_of(iatom)

         nset = bs_env%basis_set_RI(ikind)%gto_basis_set%nset
         l_max => bs_env%basis_set_RI(ikind)%gto_basis_set%lmax
         l_min => bs_env%basis_set_RI(ikind)%gto_basis_set%lmin
         nsgf_set => bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf_set

         DO iset = 1, nset
            CPASSERT(l_max(iset) == l_min(iset))
            bs_env%l_RI(i_RI + 1:i_RI + nsgf_set(iset)) = l_max(iset)
            i_RI = i_RI + nsgf_set(iset)
         END DO

      END DO
      CPASSERT(i_RI == n_RI)

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, FMT="(T2,A)") " "
         WRITE (u, FMT="(T2,2A,T75,I8)") "Number of auxiliary Gaussian basis functions ", &
            "for χ, ε, W", n_RI
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_RI_basis_and_basis_function_indices

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_chi_eps_W(bs_env, kpoints)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_chi_eps_W'

      INTEGER                                            :: exp_two, handle, i_dim, n_dim, nkp, &
                                                            nkp_extra, nkp_orig, u
      INTEGER, DIMENSION(3)                              :: nkp_grid, nkp_grid_extra, periodic
      REAL(KIND=dp)                                      :: exp_s_p, n_dim_inv

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      kpoints%kp_scheme = "GENERAL"

      periodic(1:3) = bs_env%periodic(1:3)

      DO i_dim = 1, 3

         CPASSERT(periodic(i_dim) == 0 .OR. periodic(i_dim) == 1)

         SELECT CASE (periodic(i_dim))
         CASE (0)

            nkp_grid(i_dim) = 1
            nkp_grid_extra(i_dim) = 1

         CASE (1)

            SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
            CASE (large_cell_Gamma)
               nkp_grid(i_dim) = 4
               nkp_grid_extra(i_dim) = 6
            CASE (small_cell_full_kp)
               ! for small cell, use the smallest 2^n x 2^m x 2^k k-mesh which is larger than
               ! the SCF k-point grid; extrapolate with 2^(n+1) x 2^(m+1) x 2^(k+1) k-mesh
               exp_two = 4
               DO WHILE (.TRUE.)
                  IF (exp_two > bs_env%kpoints_scf_desymm%nkp_grid(i_dim)) THEN
                     nkp_grid(i_dim) = exp_two
                     nkp_grid_extra(i_dim) = exp_two*2
                     EXIT
                  ELSE
                     exp_two = exp_two*2
                  END IF
               END DO
            END SELECT

         CASE DEFAULT

            CPABORT("Error in periodicity.")

         END SELECT

      END DO

      nkp_orig = MAX(nkp_grid(1)*nkp_grid(2)*nkp_grid(3)/2, 1)

      nkp_extra = nkp_grid_extra(1)*nkp_grid_extra(2)*nkp_grid_extra(3)/2

      nkp = nkp_orig + nkp_extra

      kpoints%nkp_grid(1:3) = nkp_grid(1:3)
      kpoints%nkp = nkp

      bs_env%nkp_grid_chi_eps_W_orig(1:3) = nkp_grid(1:3)
      bs_env%nkp_grid_chi_eps_W_extra(1:3) = nkp_grid_extra(1:3)
      bs_env%nkp_chi_eps_W_orig = nkp_orig
      bs_env%nkp_chi_eps_W_extra = nkp_extra
      bs_env%nkp_chi_eps_W_orig_plus_extra = nkp

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))
      ALLOCATE (bs_env%wkp_no_extra(nkp), bs_env%wkp_s_p(nkp))

      CALL compute_xkp(kpoints%xkp, 1, nkp_orig, nkp_grid)
      CALL compute_xkp(kpoints%xkp, nkp_orig + 1, nkp, nkp_grid_extra)

      n_dim = SUM(periodic)
      IF (n_dim == 0) THEN
         ! molecules
         kpoints%wkp(1) = 1.0_dp
         bs_env%wkp_s_p(1) = 1.0_dp
         bs_env%wkp_no_extra(1) = 1.0_dp
      ELSE

         n_dim_inv = 1.0_dp/REAL(n_dim, KIND=dp)

         ! k-point weights are chosen to automatically extrapolate the k-point mesh
         CALL compute_wkp(kpoints%wkp(1:nkp_orig), nkp_orig, nkp_extra, n_dim_inv)
         CALL compute_wkp(kpoints%wkp(nkp_orig + 1:nkp), nkp_extra, nkp_orig, n_dim_inv)

         bs_env%wkp_no_extra(1:nkp_orig) = 0.0_dp
         bs_env%wkp_no_extra(nkp_orig + 1:nkp) = 1.0_dp/REAL(nkp_extra, KIND=dp)

         IF (n_dim == 3) THEN
            ! W_PQ(k) for an s-function P and a p-function Q diverges as 1/k at k=0
            ! (instead of 1/k^2 for P and Q both being s-functions).
            exp_s_p = 2.0_dp*n_dim_inv
            CALL compute_wkp(bs_env%wkp_s_p(1:nkp_orig), nkp_orig, nkp_extra, exp_s_p)
            CALL compute_wkp(bs_env%wkp_s_p(nkp_orig + 1:nkp), nkp_extra, nkp_orig, exp_s_p)
         ELSE
            bs_env%wkp_s_p(1:nkp) = bs_env%wkp_no_extra(1:nkp)
         END IF

      END IF

      IF (bs_env%approx_kp_extrapol) THEN
         bs_env%wkp_orig = 1.0_dp/REAL(nkp_orig, KIND=dp)
      END IF

      ! heuristic parameter: how many k-points for χ, ε, and W are used simultaneously
      ! (less simultaneous k-points: less memory, but more computational effort because of
      !  recomputation of V(k))
      bs_env%nkp_chi_eps_W_batch = 4

      bs_env%num_chi_eps_W_batches = (bs_env%nkp_chi_eps_W_orig_plus_extra - 1)/ &
                                     bs_env%nkp_chi_eps_W_batch + 1

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, FMT="(T2,A)") " "
         WRITE (u, FMT="(T2,1A,T71,3I4)") "K-point mesh 1 for χ, ε, W", nkp_grid(1:3)
         WRITE (u, FMT="(T2,2A,T71,3I4)") "K-point mesh 2 for χ, ε, W ", &
            "(for k-point extrapolation of W)", nkp_grid_extra(1:3)
         WRITE (u, FMT="(T2,A,T80,L)") "Approximate the k-point extrapolation", &
            bs_env%approx_kp_extrapol
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_chi_eps_W

! **************************************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE kpoint_init_cell_index_simple(kpoints, qs_env)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_init_cell_index_simple'

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, para_env, sab_orb)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, dft_control=dft_control, sab_orb=sab_orb)
      CALL kpoint_init_cell_index(kpoints, sab_orb, para_env, dft_control)

      CALL timestop(handle)

   END SUBROUTINE kpoint_init_cell_index_simple

! **************************************************************************************************
!> \brief ...
!> \param xkp ...
!> \param ikp_start ...
!> \param ikp_end ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE compute_xkp(xkp, ikp_start, ikp_end, grid)

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      INTEGER                                            :: ikp_start, ikp_end
      INTEGER, DIMENSION(3)                              :: grid

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_xkp'

      INTEGER                                            :: handle, i, ix, iy, iz

      CALL timeset(routineN, handle)

      i = ikp_start
      DO ix = 1, grid(1)
         DO iy = 1, grid(2)
            DO iz = 1, grid(3)

               IF (i > ikp_end) CYCLE

               xkp(1, i) = REAL(2*ix - grid(1) - 1, KIND=dp)/(2._dp*REAL(grid(1), KIND=dp))
               xkp(2, i) = REAL(2*iy - grid(2) - 1, KIND=dp)/(2._dp*REAL(grid(2), KIND=dp))
               xkp(3, i) = REAL(2*iz - grid(3) - 1, KIND=dp)/(2._dp*REAL(grid(3), KIND=dp))
               i = i + 1

            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_xkp

! **************************************************************************************************
!> \brief ...
!> \param wkp ...
!> \param nkp_1 ...
!> \param nkp_2 ...
!> \param exponent ...
! **************************************************************************************************
   SUBROUTINE compute_wkp(wkp, nkp_1, nkp_2, exponent)
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER                                            :: nkp_1, nkp_2
      REAL(KIND=dp)                                      :: exponent

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_wkp'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: nkp_ratio

      CALL timeset(routineN, handle)

      nkp_ratio = REAL(nkp_2, KIND=dp)/REAL(nkp_1, KIND=dp)

      wkp(:) = 1.0_dp/REAL(nkp_1, KIND=dp)/(1.0_dp - nkp_ratio**exponent)

      CALL timestop(handle)

   END SUBROUTINE compute_wkp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_matrices(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'allocate_matrices'

      INTEGER                                            :: handle, i_t
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_tensor
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_RI_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      fm_struct => bs_env%fm_ks_Gamma(1)%matrix_struct

      CALL cp_fm_create(bs_env%fm_Gocc, fm_struct)
      CALL cp_fm_create(bs_env%fm_Gvir, fm_struct)

      NULLIFY (fm_struct_RI_global)
      CALL cp_fm_struct_create(fm_struct_RI_global, context=blacs_env, nrow_global=bs_env%n_RI, &
                               ncol_global=bs_env%n_RI, para_env=para_env)
      CALL cp_fm_create(bs_env%fm_RI_RI, fm_struct_RI_global)
      CALL cp_fm_create(bs_env%fm_chi_Gamma_freq, fm_struct_RI_global)
      CALL cp_fm_create(bs_env%fm_W_MIC_freq, fm_struct_RI_global)
      IF (bs_env%approx_kp_extrapol) THEN
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_1_extra, fm_struct_RI_global)
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_1_no_extra, fm_struct_RI_global)
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_1_extra, 0.0_dp)
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_1_no_extra, 0.0_dp)
      END IF
      CALL cp_fm_struct_release(fm_struct_RI_global)

      ! create blacs_env for subgroups of tensor operations
      NULLIFY (blacs_env_tensor)
      CALL cp_blacs_env_create(blacs_env=blacs_env_tensor, para_env=bs_env%para_env_tensor)

      ! allocate dbcsr matrices in the tensor subgroup; actually, one only needs a small
      ! subset of blocks in the tensor subgroup, however, all atomic blocks are allocated.
      ! One might think of creating a dbcsr matrix with only the blocks that are needed
      ! in the tensor subgroup
      CALL create_mat_munu(bs_env%mat_ao_ao_tensor, qs_env, bs_env%eps_atom_grid_2d_mat, &
                           blacs_env_tensor, do_ri_aux_basis=.FALSE.)

      CALL create_mat_munu(bs_env%mat_RI_RI_tensor, qs_env, bs_env%eps_atom_grid_2d_mat, &
                           blacs_env_tensor, do_ri_aux_basis=.TRUE.)

      CALL create_mat_munu(bs_env%mat_RI_RI, qs_env, bs_env%eps_atom_grid_2d_mat, &
                           blacs_env, do_ri_aux_basis=.TRUE.)

      CALL cp_blacs_env_release(blacs_env_tensor)

      NULLIFY (bs_env%mat_chi_Gamma_tau)
      CALL dbcsr_allocate_matrix_set(bs_env%mat_chi_Gamma_tau, bs_env%num_time_freq_points)

      DO i_t = 1, bs_env%num_time_freq_points
         ALLOCATE (bs_env%mat_chi_Gamma_tau(i_t)%matrix)
         CALL dbcsr_create(bs_env%mat_chi_Gamma_tau(i_t)%matrix, template=bs_env%mat_RI_RI%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_matrices

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_GW_eigenvalues(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_GW_eigenvalues'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (bs_env%eigenval_G0W0(bs_env%n_ao, bs_env%nkp_bs_and_DOS, bs_env%n_spin))
      ALLOCATE (bs_env%eigenval_HF(bs_env%n_ao, bs_env%nkp_bs_and_DOS, bs_env%n_spin))

      CALL timestop(handle)

   END SUBROUTINE allocate_GW_eigenvalues

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE create_tensors(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_tensors'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL init_interaction_radii(bs_env)

      ! split blocks does not improve load balancing/efficienfy for tensor contraction, so we go
      ! with the standard atomic blocks
      CALL create_3c_t(bs_env%t_RI_AO__AO, bs_env%para_env_tensor, "(RI AO | AO)", [1, 2], [3], &
                       bs_env%sizes_RI, bs_env%sizes_AO, &
                       create_nl_3c=.TRUE., nl_3c=bs_env%nl_3c, qs_env=qs_env)
      CALL create_3c_t(bs_env%t_RI__AO_AO, bs_env%para_env_tensor, "(RI | AO AO)", [1], [2, 3], &
                       bs_env%sizes_RI, bs_env%sizes_AO)

      CALL create_2c_t(bs_env, bs_env%sizes_RI, bs_env%sizes_AO)

      CALL timestop(handle)

   END SUBROUTINE create_tensors

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE check_sparsity_3c(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_sparsity_3c'

      INTEGER                                            :: handle, n_atom_step, RI_atom
      INTEGER(int_8)                                     :: mem, non_zero_elements_sum, nze
      REAL(dp)                                           :: max_dist_AO_atoms, occ, occupation_sum
      REAL(KIND=dp)                                      :: t1, t2
      TYPE(dbt_type)                                     :: t_3c_global
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_global_array
      TYPE(neighbor_list_3c_type)                        :: nl_3c_global

      CALL timeset(routineN, handle)

      ! check the sparsity of 3c integral tensor (µν|P); calculate maximum distance between
      ! AO atoms µ, ν where at least a single integral (µν|P) is larger than the filter threshold
      CALL create_3c_t(t_3c_global, bs_env%para_env, "(RI AO | AO)", [1, 2], [3], &
                       bs_env%sizes_RI, bs_env%sizes_AO, &
                       create_nl_3c=.TRUE., nl_3c=nl_3c_global, qs_env=qs_env)

      CALL m_memory(mem)
      CALL bs_env%para_env%max(mem)

      ALLOCATE (t_3c_global_array(1, 1))
      CALL dbt_create(t_3c_global, t_3c_global_array(1, 1))

      CALL bs_env%para_env%sync()
      t1 = m_walltime()

      occupation_sum = 0.0_dp
      non_zero_elements_sum = 0
      max_dist_AO_atoms = 0.0_dp
      n_atom_step = INT(SQRT(REAL(bs_env%n_atom, KIND=dp)))
      ! do not compute full 3c integrals at once because it may cause out of memory
      DO RI_atom = 1, bs_env%n_atom, n_atom_step

         CALL build_3c_integrals(t_3c_global_array, &
                                 bs_env%eps_filter, &
                                 qs_env, &
                                 nl_3c_global, &
                                 int_eps=bs_env%eps_filter, &
                                 basis_i=bs_env%basis_set_RI, &
                                 basis_j=bs_env%basis_set_AO, &
                                 basis_k=bs_env%basis_set_AO, &
                                 bounds_i=[RI_atom, MIN(RI_atom + n_atom_step - 1, bs_env%n_atom)], &
                                 potential_parameter=bs_env%ri_metric, &
                                 desymmetrize=.FALSE.)

         CALL dbt_filter(t_3c_global_array(1, 1), bs_env%eps_filter)

         CALL bs_env%para_env%sync()

         CALL get_tensor_occupancy(t_3c_global_array(1, 1), nze, occ)
         non_zero_elements_sum = non_zero_elements_sum + nze
         occupation_sum = occupation_sum + occ

         CALL get_max_dist_AO_atoms(t_3c_global_array(1, 1), max_dist_AO_atoms, qs_env)

         CALL dbt_clear(t_3c_global_array(1, 1))

      END DO

      t2 = m_walltime()

      bs_env%occupation_3c_int = occupation_sum
      bs_env%max_dist_AO_atoms = max_dist_AO_atoms

      CALL dbt_destroy(t_3c_global)
      CALL dbt_destroy(t_3c_global_array(1, 1))
      DEALLOCATE (t_3c_global_array)

      CALL neighbor_list_3c_destroy(nl_3c_global)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') ''
         WRITE (bs_env%unit_nr, '(T2,A,F27.1,A)') &
            'Computed 3-center integrals (µν|P), execution time', t2 - t1, ' s'
         WRITE (bs_env%unit_nr, '(T2,A,F48.3,A)') 'Percentage of non-zero (µν|P)', &
            occupation_sum*100, ' %'
         WRITE (bs_env%unit_nr, '(T2,A,F33.1,A)') 'Max. distance between µ,ν in non-zero (µν|P)', &
            max_dist_AO_atoms*angstrom, ' A'
         WRITE (bs_env%unit_nr, '(T2,2A,I20,A)') 'Required memory if storing all 3-center ', &
            'integrals (µν|P)', INT(REAL(non_zero_elements_sum, KIND=dp)*8.0E-9_dp), ' GB'
      END IF

      CALL timestop(handle)

   END SUBROUTINE check_sparsity_3c

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param sizes_RI ...
!> \param sizes_AO ...
! **************************************************************************************************
   SUBROUTINE create_2c_t(bs_env, sizes_RI, sizes_AO)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI, sizes_AO

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_2c_t'

      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      TYPE(dbt_pgrid_type)                               :: pgrid_2d

      CALL timeset(routineN, handle)

      ! inspired from rpa_im_time.F / hfx_types.F

      pdims_2d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_2d, pgrid_2d)

      CALL create_2c_tensor(bs_env%t_G, dist_1, dist_2, pgrid_2d, sizes_AO, sizes_AO, &
                            name="(AO | AO)")
      DEALLOCATE (dist_1, dist_2)
      CALL create_2c_tensor(bs_env%t_chi, dist_1, dist_2, pgrid_2d, sizes_RI, sizes_RI, &
                            name="(RI | RI)")
      DEALLOCATE (dist_1, dist_2)
      CALL create_2c_tensor(bs_env%t_W, dist_1, dist_2, pgrid_2d, sizes_RI, sizes_RI, &
                            name="(RI | RI)")
      DEALLOCATE (dist_1, dist_2)
      CALL dbt_pgrid_destroy(pgrid_2d)

      CALL timestop(handle)

   END SUBROUTINE create_2c_t

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param para_env ...
!> \param tensor_name ...
!> \param map1 ...
!> \param map2 ...
!> \param sizes_RI ...
!> \param sizes_AO ...
!> \param create_nl_3c ...
!> \param nl_3c ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE create_3c_t(tensor, para_env, tensor_name, map1, map2, sizes_RI, sizes_AO, &
                          create_nl_3c, nl_3c, qs_env)
      TYPE(dbt_type)                                     :: tensor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=12)                                  :: tensor_name
      INTEGER, DIMENSION(:)                              :: map1, map2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI, sizes_AO
      LOGICAL, OPTIONAL                                  :: create_nl_3c
      TYPE(neighbor_list_3c_type), OPTIONAL              :: nl_3c
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_3c_t'

      INTEGER                                            :: handle, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_AO_1, dist_AO_2, dist_RI
      INTEGER, DIMENSION(3)                              :: pcoord, pdims, pdims_3d
      LOGICAL                                            :: my_create_nl_3c
      TYPE(dbt_pgrid_type)                               :: pgrid_3d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(mp_cart_type)                                 :: mp_comm_t3c_2
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      pdims_3d = 0
      CALL dbt_pgrid_create(para_env, pdims_3d, pgrid_3d)
      CALL create_3c_tensor(tensor, dist_RI, dist_AO_1, dist_AO_2, &
                            pgrid_3d, sizes_RI, sizes_AO, sizes_AO, &
                            map1=map1, map2=map2, name=tensor_name)

      IF (PRESENT(create_nl_3c)) THEN
         my_create_nl_3c = create_nl_3c
      ELSE
         my_create_nl_3c = .FALSE.
      END IF

      IF (my_create_nl_3c) THEN
         CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set)
         CALL dbt_mp_environ_pgrid(pgrid_3d, pdims, pcoord)
         CALL mp_comm_t3c_2%create(pgrid_3d%mp_comm_2d, 3, pdims)
         CALL distribution_3d_create(dist_3d, dist_RI, dist_AO_1, dist_AO_2, &
                                     nkind, particle_set, mp_comm_t3c_2, own_comm=.TRUE.)

         CALL build_3c_neighbor_lists(nl_3c, &
                                      qs_env%bs_env%basis_set_RI, &
                                      qs_env%bs_env%basis_set_AO, &
                                      qs_env%bs_env%basis_set_AO, &
                                      dist_3d, qs_env%bs_env%ri_metric, &
                                      "GW_3c_nl", qs_env, own_dist=.TRUE.)
      END IF

      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
      CALL dbt_pgrid_destroy(pgrid_3d)

      CALL timestop(handle)

   END SUBROUTINE create_3c_t

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE init_interaction_radii(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_interaction_radii'

      INTEGER                                            :: handle, ibasis
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis

      CALL timeset(routineN, handle)

      DO ibasis = 1, SIZE(bs_env%basis_set_AO)

         orb_basis => bs_env%basis_set_AO(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(orb_basis, bs_env%eps_filter)

         ri_basis => bs_env%basis_set_RI(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(ri_basis, bs_env%eps_filter)

      END DO

      CALL timestop(handle)

   END SUBROUTINE init_interaction_radii

! **************************************************************************************************
!> \brief ...
!> \param t_3c_int ...
!> \param max_dist_AO_atoms ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_max_dist_AO_atoms(t_3c_int, max_dist_AO_atoms, qs_env)
      TYPE(dbt_type)                                     :: t_3c_int
      REAL(KIND=dp)                                      :: max_dist_AO_atoms
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_max_dist_AO_atoms'

      INTEGER                                            :: atom_1, atom_2, handle, num_cells
      INTEGER, DIMENSION(3)                              :: atom_ind
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp)                                      :: abs_rab
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set, para_env)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set, para_env=para_env)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_3c_int, max_dist_AO_atoms, num_cells, index_to_cell, particle_set, cell) &
!$OMP PRIVATE(iter,atom_ind,rab, abs_rab, atom_1, atom_2)
      CALL dbt_iterator_start(iter, t_3c_int)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, atom_ind)

         atom_1 = atom_ind(2)
         atom_2 = atom_ind(3)

         rab = pbc(particle_set(atom_1)%r(1:3), particle_set(atom_2)%r(1:3), cell)

         abs_rab = SQRT(rab(1)**2 + rab(2)**2 + rab(3)**2)

         max_dist_AO_atoms = MAX(max_dist_AO_atoms, abs_rab)

      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL para_env%max(max_dist_AO_atoms)

      CALL timestop(handle)

   END SUBROUTINE get_max_dist_AO_atoms

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE set_sparsity_parallelization_parameters(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_sparsity_parallelization_parameters'

      INTEGER :: handle, i_ivl, IL_ivl, j_ivl, n_atom_per_IL_ivl, n_atom_per_ivl, n_intervals_i, &
         n_intervals_inner_loop_atoms, n_intervals_j, u
      INTEGER(KIND=int_8)                                :: input_memory_per_proc

      CALL timeset(routineN, handle)

      ! heuristic parameter to prevent out of memory
      bs_env%safety_factor_memory = 0.10_dp

      input_memory_per_proc = INT(bs_env%input_memory_per_proc_GB*1.0E9_dp, KIND=int_8)

      ! choose atomic range for λ ("i_atom"), ν ("j_atom") in
      ! M_λνP(iτ) = sum_µ (µν|P) G^occ_µλ(i|τ|,k=0)
      ! N_νλQ(iτ) = sum_σ (σλ|Q) G^vir_σν(i|τ|,k=0)
      ! such that M and N fit into the memory
      n_atom_per_ivl = INT(SQRT(bs_env%safety_factor_memory*input_memory_per_proc &
                                *bs_env%group_size_tensor/24/bs_env%n_RI &
                                /SQRT(bs_env%occupation_3c_int)))/bs_env%max_AO_bf_per_atom

      n_intervals_i = (bs_env%n_atom_i - 1)/n_atom_per_ivl + 1
      n_intervals_j = (bs_env%n_atom_j - 1)/n_atom_per_ivl + 1

      bs_env%n_atom_per_interval_ij = n_atom_per_ivl
      bs_env%n_intervals_i = n_intervals_i
      bs_env%n_intervals_j = n_intervals_j

      ALLOCATE (bs_env%i_atom_intervals(2, n_intervals_i))
      ALLOCATE (bs_env%j_atom_intervals(2, n_intervals_j))

      DO i_ivl = 1, n_intervals_i
         bs_env%i_atom_intervals(1, i_ivl) = (i_ivl - 1)*n_atom_per_ivl + bs_env%atoms_i(1)
         bs_env%i_atom_intervals(2, i_ivl) = MIN(i_ivl*n_atom_per_ivl + bs_env%atoms_i(1) - 1, &
                                                 bs_env%atoms_i(2))
      END DO

      DO j_ivl = 1, n_intervals_j
         bs_env%j_atom_intervals(1, j_ivl) = (j_ivl - 1)*n_atom_per_ivl + bs_env%atoms_j(1)
         bs_env%j_atom_intervals(2, j_ivl) = MIN(j_ivl*n_atom_per_ivl + bs_env%atoms_j(1) - 1, &
                                                 bs_env%atoms_j(2))
      END DO

      ALLOCATE (bs_env%skip_Sigma_occ(n_intervals_i, n_intervals_j))
      ALLOCATE (bs_env%skip_Sigma_vir(n_intervals_i, n_intervals_j))
      bs_env%skip_Sigma_occ(:, :) = .FALSE.
      bs_env%skip_Sigma_vir(:, :) = .FALSE.

      ! choose atomic range for µ and σ ("inner loop (IL) atom") in
      ! M_λνP(iτ) = sum_µ (µν|P) G^occ_µλ(i|τ|,k=0)
      ! N_νλQ(iτ) = sum_σ (σλ|Q) G^vir_σν(i|τ|,k=0)
      n_atom_per_IL_ivl = MIN(INT(bs_env%safety_factor_memory*input_memory_per_proc &
                                  *bs_env%group_size_tensor/n_atom_per_ivl &
                                  /bs_env%max_AO_bf_per_atom &
                                  /bs_env%n_RI/8/SQRT(bs_env%occupation_3c_int) &
                                  /bs_env%max_AO_bf_per_atom), bs_env%n_atom)

      n_intervals_inner_loop_atoms = (bs_env%n_atom - 1)/n_atom_per_IL_ivl + 1

      bs_env%n_atom_per_IL_interval = n_atom_per_IL_ivl
      bs_env%n_intervals_inner_loop_atoms = n_intervals_inner_loop_atoms

      ALLOCATE (bs_env%inner_loop_atom_intervals(2, n_intervals_inner_loop_atoms))
      DO IL_ivl = 1, n_intervals_inner_loop_atoms
         bs_env%inner_loop_atom_intervals(1, IL_ivl) = (IL_ivl - 1)*n_atom_per_IL_ivl + 1
         bs_env%inner_loop_atom_intervals(2, IL_ivl) = MIN(IL_ivl*n_atom_per_IL_ivl, bs_env%n_atom)
      END DO

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,I33)') 'Number of i and j atoms in M_λνP(τ), N_νλQ(τ):', n_atom_per_ivl
         WRITE (u, '(T2,A,I18)') 'Number of inner loop atoms for µ in M_λνP = sum_µ (µν|P) G_µλ', &
            n_atom_per_IL_ivl
      END IF

      CALL timestop(handle)

   END SUBROUTINE set_sparsity_parallelization_parameters

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE check_for_restart_files(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_for_restart_files'

      CHARACTER(LEN=9)                                   :: frmt
      CHARACTER(LEN=default_string_length)               :: f_chi, f_S_n, f_S_p, f_S_x, f_W_t, &
                                                            prefix, project_name
      INTEGER                                            :: handle, i_spin, i_t_or_w, ind, n_spin, &
                                                            num_time_freq_points
      LOGICAL                                            :: chi_exists, Sigma_neg_time_exists, &
                                                            Sigma_pos_time_exists, &
                                                            Sigma_x_spin_exists, W_time_exists
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input, print_key

      CALL timeset(routineN, handle)

      num_time_freq_points = bs_env%num_time_freq_points
      n_spin = bs_env%n_spin

      ALLOCATE (bs_env%read_chi(num_time_freq_points))
      ALLOCATE (bs_env%calc_chi(num_time_freq_points))
      ALLOCATE (bs_env%Sigma_c_exists(num_time_freq_points, n_spin))

      CALL get_qs_env(qs_env, input=input)

      logger => cp_get_default_logger()
      print_key => section_vals_get_subs_vals(input, 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART')
      project_name = cp_print_key_generate_filename(logger, print_key, extension="", &
                                                    my_local=.FALSE.)
      WRITE (prefix, '(2A)') TRIM(project_name), "-RESTART_"
      bs_env%prefix = prefix

      bs_env%all_W_exist = .TRUE.

      DO i_t_or_w = 1, num_time_freq_points

         IF (i_t_or_w < 10) THEN
            WRITE (frmt, '(A)') '(3A,I1,A)'
            WRITE (f_chi, frmt) TRIM(prefix), bs_env%chi_name, "_0", i_t_or_w, ".matrix"
            WRITE (f_W_t, frmt) TRIM(prefix), bs_env%W_time_name, "_0", i_t_or_w, ".matrix"
         ELSE IF (i_t_or_w < 100) THEN
            WRITE (frmt, '(A)') '(3A,I2,A)'
            WRITE (f_chi, frmt) TRIM(prefix), bs_env%chi_name, "_", i_t_or_w, ".matrix"
            WRITE (f_W_t, frmt) TRIM(prefix), bs_env%W_time_name, "_", i_t_or_w, ".matrix"
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         INQUIRE (file=TRIM(f_chi), exist=chi_exists)
         INQUIRE (file=TRIM(f_W_t), exist=W_time_exists)

         bs_env%read_chi(i_t_or_w) = chi_exists
         bs_env%calc_chi(i_t_or_w) = .NOT. chi_exists

         bs_env%all_W_exist = bs_env%all_W_exist .AND. W_time_exists

         ! the self-energy is spin-dependent
         DO i_spin = 1, n_spin

            ind = i_t_or_w + (i_spin - 1)*num_time_freq_points

            IF (ind < 10) THEN
               WRITE (frmt, '(A)') '(3A,I1,A)'
               WRITE (f_S_p, frmt) TRIM(prefix), bs_env%Sigma_p_name, "_0", ind, ".matrix"
               WRITE (f_S_n, frmt) TRIM(prefix), bs_env%Sigma_n_name, "_0", ind, ".matrix"
            ELSE IF (i_t_or_w < 100) THEN
               WRITE (frmt, '(A)') '(3A,I2,A)'
               WRITE (f_S_p, frmt) TRIM(prefix), bs_env%Sigma_p_name, "_", ind, ".matrix"
               WRITE (f_S_n, frmt) TRIM(prefix), bs_env%Sigma_n_name, "_", ind, ".matrix"
            END IF

            INQUIRE (file=TRIM(f_S_p), exist=Sigma_pos_time_exists)
            INQUIRE (file=TRIM(f_S_n), exist=Sigma_neg_time_exists)

            bs_env%Sigma_c_exists(i_t_or_w, i_spin) = Sigma_pos_time_exists .AND. &
                                                      Sigma_neg_time_exists

         END DO

      END DO

      IF (bs_env%all_W_exist) THEN
         bs_env%read_chi(:) = .FALSE.
         bs_env%calc_chi(:) = .FALSE.
      END IF

      bs_env%Sigma_x_exists = .TRUE.
      DO i_spin = 1, n_spin
         WRITE (f_S_x, '(3A,I1,A)') TRIM(prefix), bs_env%Sigma_x_name, "_0", i_spin, ".matrix"
         INQUIRE (file=TRIM(f_S_x), exist=Sigma_x_spin_exists)
         bs_env%Sigma_x_exists = bs_env%Sigma_x_exists .AND. Sigma_x_spin_exists
      END DO

      CALL timestop(handle)

   END SUBROUTINE check_for_restart_files

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE set_parallelization_parameters(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_parallelization_parameters'

      INTEGER                                            :: color_sub, dummy_1, dummy_2, handle, &
                                                            num_pe, num_t_groups, u
      INTEGER(KIND=int_8)                                :: mem
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env)

      num_pe = para_env%num_pe
      ! if not already set, use all processors for the group (for large-cell GW, performance
      ! seems to be best for a single group with all MPI processes per group)
      IF (bs_env%group_size_tensor < 0 .OR. bs_env%group_size_tensor > num_pe) &
         bs_env%group_size_tensor = num_pe

      ! group_size_tensor must divide num_pe without rest; otherwise everything will be complicated
      IF (MODULO(num_pe, bs_env%group_size_tensor) .NE. 0) THEN
         CALL find_good_group_size(num_pe, bs_env%group_size_tensor)
      END IF

      ! para_env_tensor for tensor subgroups
      color_sub = para_env%mepos/bs_env%group_size_tensor
      bs_env%tensor_group_color = color_sub

      ALLOCATE (bs_env%para_env_tensor)
      CALL bs_env%para_env_tensor%from_split(para_env, color_sub)

      num_t_groups = para_env%num_pe/bs_env%group_size_tensor
      bs_env%num_tensor_groups = num_t_groups

      CALL get_i_j_atoms(bs_env%atoms_i, bs_env%atoms_j, bs_env%n_atom_i, bs_env%n_atom_j, &
                         color_sub, bs_env)

      ALLOCATE (bs_env%atoms_i_t_group(2, num_t_groups))
      ALLOCATE (bs_env%atoms_j_t_group(2, num_t_groups))
      DO color_sub = 0, num_t_groups - 1
         CALL get_i_j_atoms(bs_env%atoms_i_t_group(1:2, color_sub + 1), &
                            bs_env%atoms_j_t_group(1:2, color_sub + 1), &
                            dummy_1, dummy_2, color_sub, bs_env)
      END DO

      CALL m_memory(mem)
      CALL bs_env%para_env%max(mem)

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A,I47)') 'Group size for tensor operations', bs_env%group_size_tensor
      END IF

      CALL timestop(handle)

   END SUBROUTINE set_parallelization_parameters

! **************************************************************************************************
!> \brief ...
!> \param num_pe ...
!> \param group_size ...
! **************************************************************************************************
   SUBROUTINE find_good_group_size(num_pe, group_size)

      INTEGER                                            :: num_pe, group_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'find_good_group_size'

      INTEGER                                            :: group_size_minus, group_size_orig, &
                                                            group_size_plus, handle, i_diff

      CALL timeset(routineN, handle)

      group_size_orig = group_size

      DO i_diff = 1, num_pe

         group_size_minus = group_size - i_diff

         IF (MODULO(num_pe, group_size_minus) == 0 .AND. group_size_minus > 0) THEN
            group_size = group_size_minus
            EXIT
         END IF

         group_size_plus = group_size + i_diff

         IF (MODULO(num_pe, group_size_plus) == 0 .AND. group_size_plus <= num_pe) THEN
            group_size = group_size_plus
            EXIT
         END IF

      END DO

      IF (group_size_orig == group_size) CPABORT("Group size error")

      CALL timestop(handle)

   END SUBROUTINE find_good_group_size

! **************************************************************************************************
!> \brief ...
!> \param atoms_i ...
!> \param atoms_j ...
!> \param n_atom_i ...
!> \param n_atom_j ...
!> \param color_sub ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_i_j_atoms(atoms_i, atoms_j, n_atom_i, n_atom_j, color_sub, bs_env)

      INTEGER, DIMENSION(2)                              :: atoms_i, atoms_j
      INTEGER                                            :: n_atom_i, n_atom_j, color_sub
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_i_j_atoms'

      INTEGER                                            :: handle, i_atoms_per_group, i_group, &
                                                            ipcol, ipcol_loop, iprow, iprow_loop, &
                                                            j_atoms_per_group, npcol, nprow

      CALL timeset(routineN, handle)

      ! create a square mesh of tensor groups for iatom and jatom; code from blacs_env_create
      CALL square_mesh(nprow, npcol, bs_env%num_tensor_groups)

      i_group = 0
      DO ipcol_loop = 0, npcol - 1
         DO iprow_loop = 0, nprow - 1
            IF (i_group == color_sub) THEN
               iprow = iprow_loop
               ipcol = ipcol_loop
            END IF
            i_group = i_group + 1
         END DO
      END DO

      IF (MODULO(bs_env%n_atom, nprow) == 0) THEN
         i_atoms_per_group = bs_env%n_atom/nprow
      ELSE
         i_atoms_per_group = bs_env%n_atom/nprow + 1
      END IF

      IF (MODULO(bs_env%n_atom, npcol) == 0) THEN
         j_atoms_per_group = bs_env%n_atom/npcol
      ELSE
         j_atoms_per_group = bs_env%n_atom/npcol + 1
      END IF

      atoms_i(1) = iprow*i_atoms_per_group + 1
      atoms_i(2) = MIN((iprow + 1)*i_atoms_per_group, bs_env%n_atom)
      n_atom_i = atoms_i(2) - atoms_i(1) + 1

      atoms_j(1) = ipcol*j_atoms_per_group + 1
      atoms_j(2) = MIN((ipcol + 1)*j_atoms_per_group, bs_env%n_atom)
      n_atom_j = atoms_j(2) - atoms_j(1) + 1

      CALL timestop(handle)

   END SUBROUTINE get_i_j_atoms

! **************************************************************************************************
!> \brief ...
!> \param nprow ...
!> \param npcol ...
!> \param nproc ...
! **************************************************************************************************
   SUBROUTINE square_mesh(nprow, npcol, nproc)
      INTEGER                                            :: nprow, npcol, nproc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'square_mesh'

      INTEGER                                            :: gcd_max, handle, ipe, jpe

      CALL timeset(routineN, handle)

      gcd_max = -1
      DO ipe = 1, CEILING(SQRT(REAL(nproc, dp)))
         jpe = nproc/ipe
         IF (ipe*jpe .NE. nproc) CYCLE
         IF (gcd(ipe, jpe) >= gcd_max) THEN
            nprow = ipe
            npcol = jpe
            gcd_max = gcd(ipe, jpe)
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE square_mesh

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE set_heuristic_parameters(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_heuristic_parameters'

      INTEGER                                            :: handle, u
      LOGICAL                                            :: do_BvK_cell

      CALL timeset(routineN, handle)

      ! Determines number of cells used for summing the cells R in the Coulomb matrix,
      ! V_PQ(k) = \sum_R <P,cell=0 | 1/r | Q,cell=R>. SIZE_LATTICE_SUM_V 3 gives
      ! good convergence
      bs_env%size_lattice_sum_V = 3

      ! for generating numerically stable minimax Fourier integration weights
      bs_env%num_points_per_magnitude = 200

      ! for periodic systems and for 20 minimax points, we use a regularized minimax mesh
      ! (from experience: regularized minimax meshes converges faster for periodic systems
      !  and for 20 pts)
      IF (SUM(bs_env%periodic) .NE. 0 .OR. bs_env%num_time_freq_points == 20) THEN
         bs_env%regularization_minimax = 1.0E-6_dp
      ELSE
         bs_env%regularization_minimax = 0.0_dp
      END IF

      bs_env%stabilize_exp = 70.0_dp
      bs_env%eps_atom_grid_2d_mat = 1.0E-50_dp

      ! only use interval ω in [0, 1 Ha] (1 Hartree = 27.211 eV) for virt, and ω in [-1 Ha, 0]
      ! for occ for use in analytic continuation of self-energy Σ^c_n(iω,k) -> Σ^c_n(ϵ,k)
      bs_env%freq_max_fit = 1.0_dp

      ! use a 16-parameter Padé fit
      bs_env%nparam_pade = 16

      ! resolution of the identity with the truncated Coulomb metric, cutoff radius 3 Angström
      bs_env%ri_metric%potential_type = do_potential_truncated
      bs_env%ri_metric%omega = 0.0_dp
      ! cutoff radius = 3 Angström
      bs_env%ri_metric%cutoff_radius = 3.0_dp/angstrom
      bs_env%ri_metric%filename = "t_c_g.dat"

      bs_env%eps_eigval_mat_RI = 0.0_dp

      ! for periodic systems, we use the regularized resolution of the identity
      IF (SUM(bs_env%periodic) == 0) THEN
         bs_env%regularization_RI = 0.0_dp
      ELSE
         bs_env%regularization_RI = 1.0E-2_dp
      END IF

      ! truncated Coulomb operator for exchange self-energy
      ! (see details in Guidon, VandeVondele, Hutter, JCTC 5, 3010 (2009) and references therein)
      do_BvK_cell = bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp
      CALL trunc_coulomb_for_exchange(qs_env, bs_env%trunc_coulomb, &
                                      rel_cutoff_trunc_coulomb_ri_x=0.5_dp, &
                                      cell_grid=bs_env%cell_grid_scf_desymm, &
                                      do_BvK_cell=do_BvK_cell)

      ! for small-cell GW, we need more cells than normally used by the filter bs_env%eps_filter
      ! (in particular for computing the self-energy because of higher number of cells needed)
      bs_env%heuristic_filter_factor = 1.0E-4

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, FMT="(T2,2A,F21.1,A)") "Cutoff radius for the truncated Coulomb ", &
            "operator in Σ^x:", bs_env%trunc_coulomb%cutoff_radius*angstrom, " Å"
      END IF

      CALL timestop(handle)

   END SUBROUTINE set_heuristic_parameters

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_header_and_input_parameters(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_header_and_input_parameters'

      INTEGER                                            :: handle, u

      CALL timeset(routineN, handle)

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, *) ' '
         WRITE (u, '(T2,2A)') '------------------------------------------------', &
            '-------------------------------'
         WRITE (u, '(T2,2A)') '-                                               ', &
            '                              -'
         WRITE (u, '(T2,2A)') '-                              GW CALCULATION   ', &
            '                              -'
         WRITE (u, '(T2,2A)') '-                                               ', &
            '                              -'
         WRITE (u, '(T2,2A)') '------------------------------------------------', &
            '-------------------------------'
         WRITE (u, '(T2,A)') ' '
         WRITE (u, '(T2,A,I45)') 'Input: Number of time/freq. points', bs_env%num_time_freq_points
         WRITE (u, '(T2,A,ES27.1)') 'Input: Filter threshold for sparse tensor operations', &
            bs_env%eps_filter
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_header_and_input_parameters

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_V_xc(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_V_xc'

      INTEGER                                            :: handle, img, ispin, myfun, nimages
      REAL(KIND=dp)                                      :: energy_ex, energy_exc, energy_total
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_ks_without_v_xc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input, energy=energy, dft_control=dft_control)

      ! previously, dft_control%nimages set to # neighbor cells, revert for Γ-only KS matrix
      nimages = dft_control%nimages
      dft_control%nimages = bs_env%nimages_scf

      ! we need to reset XC functional, therefore, get XC input
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=xc_none)

      ! save the energy before the energy gets updated
      energy_total = energy%total
      energy_exc = energy%exc
      energy_ex = energy%ex

      SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
      CASE (large_cell_Gamma)

         NULLIFY (mat_ks_without_v_xc)
         CALL dbcsr_allocate_matrix_set(mat_ks_without_v_xc, bs_env%n_spin)

         DO ispin = 1, bs_env%n_spin
            ALLOCATE (mat_ks_without_v_xc(ispin)%matrix)
            CALL dbcsr_create(mat_ks_without_v_xc(ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
         END DO

         ! calculate KS-matrix without XC
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                           ext_ks_matrix=mat_ks_without_v_xc)

         DO ispin = 1, bs_env%n_spin
            ! transfer dbcsr matrix to fm
            CALL cp_fm_create(bs_env%fm_V_xc_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
            CALL copy_dbcsr_to_fm(mat_ks_without_v_xc(ispin)%matrix, bs_env%fm_V_xc_Gamma(ispin))

            ! v_xc = h_ks - h_ks(v_xc = 0)
            CALL cp_fm_scale_and_add(alpha=-1.0_dp, matrix_a=bs_env%fm_V_xc_Gamma(ispin), &
                                     beta=1.0_dp, matrix_b=bs_env%fm_ks_Gamma(ispin))
         END DO

         CALL dbcsr_deallocate_matrix_set(mat_ks_without_v_xc)

      CASE (small_cell_full_kp)

         ! calculate KS-matrix without XC
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)
         CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks_kp)

         ALLOCATE (bs_env%fm_V_xc_R(dft_control%nimages, bs_env%n_spin))
         DO ispin = 1, bs_env%n_spin
            DO img = 1, dft_control%nimages
               ! safe fm_V_xc_R in fm_matrix because saving in dbcsr matrix caused trouble...
               CALL copy_dbcsr_to_fm(matrix_ks_kp(ispin, img)%matrix, bs_env%fm_work_mo(1))
               CALL cp_fm_create(bs_env%fm_V_xc_R(img, ispin), bs_env%fm_work_mo(1)%matrix_struct)
               ! store h_ks(v_xc = 0) in fm_V_xc_R
               CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=bs_env%fm_V_xc_R(img, ispin), &
                                        beta=1.0_dp, matrix_b=bs_env%fm_work_mo(1))
            END DO
         END DO

      END SELECT

      ! set back the energy
      energy%total = energy_total
      energy%exc = energy_exc
      energy%ex = energy_ex

      ! set back nimages
      dft_control%nimages = nimages

      ! set the DFT functional and HF fraction back
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)

      IF (bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp) THEN
         ! calculate KS-matrix again with XC
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)
         DO ispin = 1, bs_env%n_spin
            DO img = 1, dft_control%nimages
               ! store h_ks in fm_work_mo
               CALL copy_dbcsr_to_fm(matrix_ks_kp(ispin, img)%matrix, bs_env%fm_work_mo(1))
               ! v_xc = h_ks - h_ks(v_xc = 0)
               CALL cp_fm_scale_and_add(alpha=-1.0_dp, matrix_a=bs_env%fm_V_xc_R(img, ispin), &
                                        beta=1.0_dp, matrix_b=bs_env%fm_work_mo(1))
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_V_xc

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_time_and_frequency_minimax_grid(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_time_and_frequency_minimax_grid'

      INTEGER                                            :: handle, homo, i_w, ierr, ispin, j_w, &
                                                            n_mo, num_time_freq_points, u
      REAL(KIND=dp)                                      :: E_max, E_max_ispin, E_min, E_min_ispin, &
                                                            E_range, max_error_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: points_and_weights

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao
      num_time_freq_points = bs_env%num_time_freq_points

      ALLOCATE (bs_env%imag_freq_points(num_time_freq_points))
      ALLOCATE (bs_env%imag_time_points(num_time_freq_points))
      ALLOCATE (bs_env%weights_cos_t_to_w(num_time_freq_points, num_time_freq_points))
      ALLOCATE (bs_env%weights_cos_w_to_t(num_time_freq_points, num_time_freq_points))
      ALLOCATE (bs_env%weights_sin_t_to_w(num_time_freq_points, num_time_freq_points))

      ! minimum and maximum difference between eigenvalues of unoccupied and an occupied MOs
      E_min = 1000.0_dp
      E_max = -1000.0_dp
      DO ispin = 1, bs_env%n_spin
         homo = bs_env%n_occ(ispin)
         SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
         CASE (large_cell_Gamma)
            E_min_ispin = bs_env%eigenval_scf_Gamma(homo + 1, ispin) - &
                          bs_env%eigenval_scf_Gamma(homo, ispin)
            E_max_ispin = bs_env%eigenval_scf_Gamma(n_mo, ispin) - &
                          bs_env%eigenval_scf_Gamma(1, ispin)
         CASE (small_cell_full_kp)
            E_min_ispin = MINVAL(bs_env%eigenval_scf(homo + 1, :, ispin)) - &
                          MAXVAL(bs_env%eigenval_scf(homo, :, ispin))
            E_max_ispin = MAXVAL(bs_env%eigenval_scf(n_mo, :, ispin)) - &
                          MINVAL(bs_env%eigenval_scf(1, :, ispin))
         END SELECT
         E_min = MIN(E_min, E_min_ispin)
         E_max = MAX(E_max, E_max_ispin)
      END DO

      E_range = E_max/E_min

      ALLOCATE (points_and_weights(2*num_time_freq_points))

      ! frequency points
      IF (num_time_freq_points .LE. 20) THEN
         CALL get_rpa_minimax_coeff(num_time_freq_points, E_range, points_and_weights, ierr, .FALSE.)
      ELSE
         CALL get_rpa_minimax_coeff_larger_grid(num_time_freq_points, E_range, points_and_weights)
      END IF

      ! one needs to scale the minimax grids, see Azizi, Wilhelm, Golze, Panades-Barrueta,
      ! Giantomassi, Rinke, Draxl, Gonze et al., 2 publications
      bs_env%imag_freq_points(:) = points_and_weights(1:num_time_freq_points)*E_min

      ! determine number of fit points in the interval [0,ω_max] for virt, or [-ω_max,0] for occ
      bs_env%num_freq_points_fit = 0
      DO i_w = 1, num_time_freq_points
         IF (bs_env%imag_freq_points(i_w) < bs_env%freq_max_fit) THEN
            bs_env%num_freq_points_fit = bs_env%num_freq_points_fit + 1
         END IF
      END DO

      ! iω values for the analytic continuation Σ^c_n(iω,k) -> Σ^c_n(ϵ,k)
      ALLOCATE (bs_env%imag_freq_points_fit(bs_env%num_freq_points_fit))
      j_w = 0
      DO i_w = 1, num_time_freq_points
         IF (bs_env%imag_freq_points(i_w) < bs_env%freq_max_fit) THEN
            j_w = j_w + 1
            bs_env%imag_freq_points_fit(j_w) = bs_env%imag_freq_points(i_w)
         END IF
      END DO

      ! reset the number of Padé parameters if smaller than the number of
      ! imaginary-frequency points for the fit
      IF (bs_env%num_freq_points_fit < bs_env%nparam_pade) THEN
         bs_env%nparam_pade = bs_env%num_freq_points_fit
      END IF

      ! time points
      IF (num_time_freq_points .LE. 20) THEN
         CALL get_exp_minimax_coeff(num_time_freq_points, E_range, points_and_weights)
      ELSE
         CALL get_exp_minimax_coeff_gw(num_time_freq_points, E_range, points_and_weights)
      END IF

      bs_env%imag_time_points(:) = points_and_weights(1:num_time_freq_points)/(2.0_dp*E_min)

      DEALLOCATE (points_and_weights)

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,F55.2)') 'SCF direct band gap (eV)', E_min*evolt
         WRITE (u, '(T2,A,F53.2)') 'Max. SCF eigval diff. (eV)', E_max*evolt
         WRITE (u, '(T2,A,F55.2)') 'E-Range for minimax grid', E_range
         WRITE (u, '(T2,A,I27)') 'Number of Padé parameters for analytic continuation:', &
            bs_env%nparam_pade
         WRITE (u, '(T2,A)') ''
      END IF

      ! in minimax grids, Fourier transforms t -> w and w -> t are split using
      ! e^(iwt) = cos(wt) + i sin(wt); we thus calculate weights for trafos with a cos and
      ! sine prefactor; details in Azizi, Wilhelm, Golze, Giantomassi, Panades-Barrueta,
      ! Rinke, Draxl, Gonze et al., 2 publications

      ! cosine transform weights imaginary time to imaginary frequency
      CALL get_l_sq_wghts_cos_tf_t_to_w(num_time_freq_points, &
                                        bs_env%imag_time_points, &
                                        bs_env%weights_cos_t_to_w, &
                                        bs_env%imag_freq_points, &
                                        E_min, E_max, max_error_min, &
                                        bs_env%num_points_per_magnitude, &
                                        bs_env%regularization_minimax)

      ! cosine transform weights imaginary frequency to imaginary time
      CALL get_l_sq_wghts_cos_tf_w_to_t(num_time_freq_points, &
                                        bs_env%imag_time_points, &
                                        bs_env%weights_cos_w_to_t, &
                                        bs_env%imag_freq_points, &
                                        E_min, E_max, max_error_min, &
                                        bs_env%num_points_per_magnitude, &
                                        bs_env%regularization_minimax)

      ! sine transform weights imaginary time to imaginary frequency
      CALL get_l_sq_wghts_sin_tf_t_to_w(num_time_freq_points, &
                                        bs_env%imag_time_points, &
                                        bs_env%weights_sin_t_to_w, &
                                        bs_env%imag_freq_points, &
                                        E_min, E_max, max_error_min, &
                                        bs_env%num_points_per_magnitude, &
                                        bs_env%regularization_minimax)

      CALL timestop(handle)

   END SUBROUTINE setup_time_and_frequency_minimax_grid

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_cells_3c(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'setup_cells_3c'

      INTEGER :: atom_i, atom_j, atom_k, cell_pair_count, handle, i, i_cell_x, i_cell_x_max, &
         i_cell_x_min, i_size, ikind, img, j, j_cell, j_cell_max, j_cell_y, j_cell_y_max, &
         j_cell_y_min, j_size, k_cell, k_cell_max, k_cell_z, k_cell_z_max, k_cell_z_min, k_size, &
         nimage_pairs_3c, nimages_3c, nimages_3c_max, nkind, u
      INTEGER(KIND=int_8)                                :: mem_occ_per_proc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_other_3c_images_max
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_3c_max, nblocks_3c_max
      INTEGER, DIMENSION(3)                              :: cell_index, n_max
      REAL(KIND=dp) :: avail_mem_per_proc_GB, cell_dist, cell_radius_3c, eps, exp_min_ao, &
         exp_min_RI, frobenius_norm, mem_3c_GB, mem_occ_per_proc_GB, radius_ao, radius_ao_product, &
         radius_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: exp_ao, exp_RI

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind)

      exp_min_RI = 10.0_dp
      exp_min_ao = 10.0_dp

      DO ikind = 1, nkind

         CALL get_gto_basis_set(bs_env%basis_set_RI(ikind)%gto_basis_set, zet=exp_RI)
         CALL get_gto_basis_set(bs_env%basis_set_ao(ikind)%gto_basis_set, zet=exp_ao)

         ! we need to remove all exponents lower than a lower bound, e.g. 1E-3, because
         ! for contracted basis sets, there might be exponents = 0 in zet
         DO i = 1, SIZE(exp_RI, 1)
            DO j = 1, SIZE(exp_RI, 2)
               IF (exp_RI(i, j) < exp_min_RI .AND. exp_RI(i, j) > 1E-3_dp) exp_min_RI = exp_RI(i, j)
            END DO
         END DO
         DO i = 1, SIZE(exp_ao, 1)
            DO j = 1, SIZE(exp_ao, 2)
               IF (exp_ao(i, j) < exp_min_ao .AND. exp_ao(i, j) > 1E-3_dp) exp_min_ao = exp_ao(i, j)
            END DO
         END DO

      END DO

      eps = bs_env%eps_filter*bs_env%heuristic_filter_factor

      radius_ao = SQRT(-LOG(eps)/exp_min_ao)
      radius_ao_product = SQRT(-LOG(eps)/(2.0_dp*exp_min_ao))
      radius_RI = SQRT(-LOG(eps)/exp_min_RI)

      ! For a 3c integral (μR υS | P0) we have that cell R and cell S need to be within radius_3c
      cell_radius_3c = radius_ao_product + radius_RI + bs_env%ri_metric%cutoff_radius

      n_max(1:3) = bs_env%periodic(1:3)*30

      nimages_3c_max = 0

      i_cell_x_min = 0
      i_cell_x_max = 0
      j_cell_y_min = 0
      j_cell_y_max = 0
      k_cell_z_min = 0
      k_cell_z_max = 0

      DO i_cell_x = -n_max(1), n_max(1)
         DO j_cell_y = -n_max(2), n_max(2)
            DO k_cell_z = -n_max(3), n_max(3)

               cell_index(1:3) = (/i_cell_x, j_cell_y, k_cell_z/)

               CALL get_cell_dist(cell_index, bs_env%hmat, cell_dist)

               IF (cell_dist < cell_radius_3c) THEN
                  nimages_3c_max = nimages_3c_max + 1
                  i_cell_x_min = MIN(i_cell_x_min, i_cell_x)
                  i_cell_x_max = MAX(i_cell_x_max, i_cell_x)
                  j_cell_y_min = MIN(j_cell_y_min, j_cell_y)
                  j_cell_y_max = MAX(j_cell_y_max, j_cell_y)
                  k_cell_z_min = MIN(k_cell_z_min, k_cell_z)
                  k_cell_z_max = MAX(k_cell_z_max, k_cell_z)
               END IF

            END DO
         END DO
      END DO

      ! get index_to_cell_3c_max for the maximum possible cell range;
      ! compute 3c integrals later in this routine and check really which cell is needed
      ALLOCATE (index_to_cell_3c_max(nimages_3c_max, 3))

      img = 0
      DO i_cell_x = -n_max(1), n_max(1)
         DO j_cell_y = -n_max(2), n_max(2)
            DO k_cell_z = -n_max(3), n_max(3)

               cell_index(1:3) = (/i_cell_x, j_cell_y, k_cell_z/)

               CALL get_cell_dist(cell_index, bs_env%hmat, cell_dist)

               IF (cell_dist < cell_radius_3c) THEN
                  img = img + 1
                  index_to_cell_3c_max(img, 1:3) = cell_index(1:3)
               END IF

            END DO
         END DO
      END DO

      ! get pairs of R and S which have non-zero 3c integral (μR υS | P0)
      ALLOCATE (nblocks_3c_max(nimages_3c_max, nimages_3c_max))
      nblocks_3c_max(:, :) = 0

      cell_pair_count = 0
      DO j_cell = 1, nimages_3c_max
         DO k_cell = 1, nimages_3c_max

            cell_pair_count = cell_pair_count + 1

            ! trivial parallelization over cell pairs
            IF (MODULO(cell_pair_count, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

            DO atom_j = 1, bs_env%n_atom
               DO atom_k = 1, bs_env%n_atom
                  DO atom_i = 1, bs_env%n_atom

                     j_size = bs_env%i_ao_end_from_atom(atom_j) - bs_env%i_ao_start_from_atom(atom_j) + 1
                     k_size = bs_env%i_ao_end_from_atom(atom_k) - bs_env%i_ao_start_from_atom(atom_k) + 1
                     i_size = bs_env%i_RI_end_from_atom(atom_i) - bs_env%i_RI_start_from_atom(atom_i) + 1

                     ALLOCATE (int_3c(j_size, k_size, i_size))

                     ! compute 3-c int. ( μ(atom j) R , ν (atom k) S | P (atom i) 0 )
                     ! ("|": truncated Coulomb operator), inside build_3c_integrals: (j k | i)
                     CALL build_3c_integral_block(int_3c, qs_env, bs_env%ri_metric, &
                                                  basis_j=bs_env%basis_set_AO, &
                                                  basis_k=bs_env%basis_set_AO, &
                                                  basis_i=bs_env%basis_set_RI, &
                                                  cell_j=index_to_cell_3c_max(j_cell, 1:3), &
                                                  cell_k=index_to_cell_3c_max(k_cell, 1:3), &
                                                  atom_k=atom_k, atom_j=atom_j, atom_i=atom_i)

                     frobenius_norm = SQRT(SUM(int_3c(:, :, :)**2))

                     DEALLOCATE (int_3c)

                     ! we use a higher threshold here to safe memory when storing the 3c integrals
                     ! in every tensor group
                     IF (frobenius_norm > eps) THEN
                        nblocks_3c_max(j_cell, k_cell) = nblocks_3c_max(j_cell, k_cell) + 1
                     END IF

                  END DO
               END DO
            END DO

         END DO
      END DO

      CALL bs_env%para_env%sum(nblocks_3c_max)

      ALLOCATE (n_other_3c_images_max(nimages_3c_max))
      n_other_3c_images_max(:) = 0

      nimages_3c = 0
      nimage_pairs_3c = 0

      DO j_cell = 1, nimages_3c_max
         DO k_cell = 1, nimages_3c_max
            IF (nblocks_3c_max(j_cell, k_cell) > 0) THEN
               n_other_3c_images_max(j_cell) = n_other_3c_images_max(j_cell) + 1
               nimage_pairs_3c = nimage_pairs_3c + 1
            END IF
         END DO

         IF (n_other_3c_images_max(j_cell) > 0) nimages_3c = nimages_3c + 1

      END DO

      bs_env%nimages_3c = nimages_3c
      ALLOCATE (bs_env%index_to_cell_3c(nimages_3c, 3))
      ALLOCATE (bs_env%cell_to_index_3c(i_cell_x_min:i_cell_x_max, &
                                        j_cell_y_min:j_cell_y_max, &
                                        k_cell_z_min:k_cell_z_max))
      bs_env%cell_to_index_3c(:, :, :) = -1

      ALLOCATE (bs_env%nblocks_3c(nimages_3c, nimages_3c))
      bs_env%nblocks_3c(nimages_3c, nimages_3c) = 0

      j_cell = 0
      DO j_cell_max = 1, nimages_3c_max
         IF (n_other_3c_images_max(j_cell_max) == 0) CYCLE
         j_cell = j_cell + 1
         cell_index(1:3) = index_to_cell_3c_max(j_cell_max, 1:3)
         bs_env%index_to_cell_3c(j_cell, 1:3) = cell_index(1:3)
         bs_env%cell_to_index_3c(cell_index(1), cell_index(2), cell_index(3)) = j_cell

         k_cell = 0
         DO k_cell_max = 1, nimages_3c_max
            IF (n_other_3c_images_max(k_cell_max) == 0) CYCLE
            k_cell = k_cell + 1

            bs_env%nblocks_3c(j_cell, k_cell) = nblocks_3c_max(j_cell_max, k_cell_max)
         END DO

      END DO

      ! we use: 8*10^-9 GB / double precision number
      mem_3c_GB = REAL(bs_env%n_RI, KIND=dp)*REAL(bs_env%n_ao, KIND=dp)**2 &
                  *REAL(nimage_pairs_3c, KIND=dp)*8E-9_dp

      CALL m_memory(mem_occ_per_proc)
      CALL bs_env%para_env%max(mem_occ_per_proc)

      mem_occ_per_proc_GB = REAL(mem_occ_per_proc, KIND=dp)/1.0E9_dp

      ! number of processors per group that entirely stores the 3c integrals and does tensor ops
      avail_mem_per_proc_GB = bs_env%input_memory_per_proc_GB - mem_occ_per_proc_GB

      ! careful: downconvering real to integer, 1.9 -> 1; thus add 1.0 for upconversion, 1.9 -> 2
      bs_env%group_size_tensor = MAX(INT(mem_3c_GB/avail_mem_per_proc_GB + 1.0_dp), 1)

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, FMT="(T2,A,F52.1,A)") "Radius of atomic orbitals", radius_ao*angstrom, " Å"
         WRITE (u, FMT="(T2,A,F55.1,A)") "Radius of RI functions", radius_RI*angstrom, " Å"
         WRITE (u, FMT="(T2,A,I47)") "Number of cells for 3c integrals", nimages_3c
         WRITE (u, FMT="(T2,A,I42)") "Number of cell pairs for 3c integrals", nimage_pairs_3c
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,F37.1,A)') 'Input: Available memory per MPI process', &
            bs_env%input_memory_per_proc_GB, ' GB'
         WRITE (u, '(T2,A,F35.1,A)') 'Used memory per MPI process before GW run', &
            mem_occ_per_proc_GB, ' GB'
         WRITE (u, '(T2,A,F44.1,A)') 'Memory of three-center integrals', mem_3c_GB, ' GB'
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_cells_3c

! **************************************************************************************************
!> \brief ...
!> \param index_to_cell_1 ...
!> \param index_to_cell_2 ...
!> \param nimages_1 ...
!> \param nimages_2 ...
!> \param index_to_cell ...
!> \param cell_to_index ...
!> \param nimages ...
! **************************************************************************************************
   SUBROUTINE sum_two_R_grids(index_to_cell_1, index_to_cell_2, nimages_1, nimages_2, &
                              index_to_cell, cell_to_index, nimages)

      INTEGER, DIMENSION(:, :)                           :: index_to_cell_1, index_to_cell_2
      INTEGER                                            :: nimages_1, nimages_2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER                                            :: nimages

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sum_two_R_grids'

      INTEGER                                            :: handle, i_dim, img_1, img_2, nimages_max
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_tmp
      INTEGER, DIMENSION(3)                              :: cell_1, cell_2, R, R_max, R_min

      CALL timeset(routineN, handle)

      DO i_dim = 1, 3
         R_min(i_dim) = MINVAL(index_to_cell_1(:, i_dim)) + MINVAL(index_to_cell_2(:, i_dim))
         R_max(i_dim) = MAXVAL(index_to_cell_1(:, i_dim)) + MAXVAL(index_to_cell_2(:, i_dim))
      END DO

      nimages_max = (R_max(1) - R_min(1) + 1)*(R_max(2) - R_min(2) + 1)*(R_max(3) - R_min(3) + 1)

      ALLOCATE (index_to_cell_tmp(nimages_max, 3))
      index_to_cell_tmp(:, :) = -1

      ALLOCATE (cell_to_index(R_min(1):R_max(1), R_min(2):R_max(2), R_min(3):R_max(3)))
      cell_to_index(:, :, :) = -1

      nimages = 0

      DO img_1 = 1, nimages_1

         DO img_2 = 1, nimages_2

            cell_1(1:3) = index_to_cell_1(img_1, 1:3)
            cell_2(1:3) = index_to_cell_2(img_2, 1:3)

            R(1:3) = cell_1(1:3) + cell_2(1:3)

            ! check whether we have found a new cell
            IF (cell_to_index(R(1), R(2), R(3)) == -1) THEN

               nimages = nimages + 1
               cell_to_index(R(1), R(2), R(3)) = nimages
               index_to_cell_tmp(nimages, 1:3) = R(1:3)

            END IF

         END DO

      END DO

      ALLOCATE (index_to_cell(nimages, 3))
      index_to_cell(:, :) = index_to_cell_tmp(1:nimages, 1:3)

      CALL timestop(handle)

   END SUBROUTINE sum_two_R_grids

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_3c_integrals(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_integrals'

      INTEGER                                            :: handle, j_cell, k_cell, nimages_3c

      CALL timeset(routineN, handle)

      nimages_3c = bs_env%nimages_3c
      ALLOCATE (bs_env%t_3c_int(nimages_3c, nimages_3c))
      DO j_cell = 1, nimages_3c
         DO k_cell = 1, nimages_3c
            CALL dbt_create(bs_env%t_RI_AO__AO, bs_env%t_3c_int(j_cell, k_cell))
         END DO
      END DO

      CALL build_3c_integrals(bs_env%t_3c_int, &
                              bs_env%eps_filter, &
                              qs_env, &
                              bs_env%nl_3c, &
                              int_eps=bs_env%eps_filter*0.05_dp, &
                              basis_i=bs_env%basis_set_RI, &
                              basis_j=bs_env%basis_set_AO, &
                              basis_k=bs_env%basis_set_AO, &
                              potential_parameter=bs_env%ri_metric, &
                              desymmetrize=.FALSE., do_kpoints=.TRUE., cell_sym=.TRUE., &
                              cell_to_index_ext=bs_env%cell_to_index_3c)

      CALL bs_env%para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE compute_3c_integrals

! **************************************************************************************************
!> \brief ...
!> \param cell_index ...
!> \param hmat ...
!> \param cell_dist ...
! **************************************************************************************************
   SUBROUTINE get_cell_dist(cell_index, hmat, cell_dist)

      INTEGER, DIMENSION(3)                              :: cell_index
      REAL(KIND=dp)                                      :: hmat(3, 3), cell_dist

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_cell_dist'

      INTEGER                                            :: handle, i_dim
      INTEGER, DIMENSION(3)                              :: cell_index_adj
      REAL(KIND=dp)                                      :: cell_dist_3(3)

      CALL timeset(routineN, handle)

      ! the distance of cells needs to be taken to adjacent neighbors, not
      ! between the center of the cells. We thus need to rescale the cell index
      DO i_dim = 1, 3
         IF (cell_index(i_dim) > 0) cell_index_adj(i_dim) = cell_index(i_dim) - 1
         IF (cell_index(i_dim) < 0) cell_index_adj(i_dim) = cell_index(i_dim) + 1
         IF (cell_index(i_dim) == 0) cell_index_adj(i_dim) = cell_index(i_dim)
      END DO

      cell_dist_3(1:3) = MATMUL(hmat, REAL(cell_index_adj, KIND=dp))

      cell_dist = SQRT(ABS(SUM(cell_dist_3(1:3)**2)))

      CALL timestop(handle)

   END SUBROUTINE get_cell_dist

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param kpoints ...
!> \param do_print ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_scf_desymm(qs_env, bs_env, kpoints, do_print)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_scf_desymm'

      INTEGER                                            :: handle, i_cell_x, i_dim, img, j_cell_y, &
                                                            k_cell_z, nimages, nkp, u
      INTEGER, DIMENSION(3)                              :: cell_grid, cixd, nkp_grid
      TYPE(kpoint_type), POINTER                         :: kpoints_scf

      LOGICAL:: do_print

      CALL timeset(routineN, handle)

      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      CALL get_qs_env(qs_env=qs_env, kpoints=kpoints_scf)

      nkp_grid(1:3) = kpoints_scf%nkp_grid(1:3)
      nkp = nkp_grid(1)*nkp_grid(2)*nkp_grid(3)

      ! we need in periodic directions at least 2 k-points in the SCF
      DO i_dim = 1, 3
         IF (bs_env%periodic(i_dim) == 1) THEN
            CPASSERT(nkp_grid(i_dim) > 1)
         END IF
      END DO

      kpoints%kp_scheme = "GENERAL"
      kpoints%nkp_grid(1:3) = nkp_grid(1:3)
      kpoints%nkp = nkp
      bs_env%nkp_scf_desymm = nkp

      ALLOCATE (kpoints%xkp(1:3, nkp))
      CALL compute_xkp(kpoints%xkp, 1, nkp, nkp_grid)

      ALLOCATE (kpoints%wkp(nkp))
      kpoints%wkp(:) = 1.0_dp/REAL(nkp, KIND=dp)

      ! for example 4x3x6 kpoint grid -> 3x3x5 cell grid because we need the same number of
      ! neighbor cells on both sides of the unit cell
      cell_grid(1:3) = nkp_grid(1:3) - MODULO(nkp_grid(1:3) + 1, 2)
      ! cell index: for example for x: from -n_x/2 to +n_x/2, n_x: number of cells in x direction
      cixd(1:3) = cell_grid(1:3)/2

      nimages = cell_grid(1)*cell_grid(2)*cell_grid(3)

      bs_env%nimages_scf_desymm = nimages

      ALLOCATE (kpoints%cell_to_index(-cixd(1):cixd(1), -cixd(2):cixd(2), -cixd(3):cixd(3)))
      ALLOCATE (kpoints%index_to_cell(nimages, 3))

      img = 0
      DO i_cell_x = -cixd(1), cixd(1)
         DO j_cell_y = -cixd(2), cixd(2)
            DO k_cell_z = -cixd(3), cixd(3)
               img = img + 1
               kpoints%cell_to_index(i_cell_x, j_cell_y, k_cell_z) = img
               kpoints%index_to_cell(img, 1:3) = (/i_cell_x, j_cell_y, k_cell_z/)
            END DO
         END DO
      END DO

      u = bs_env%unit_nr
      IF (u > 0 .AND. do_print) THEN
         WRITE (u, FMT="(T2,A,I49)") "Number of cells for G, χ, W, Σ", nimages
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_scf_desymm

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_cells_Delta_R(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_cells_Delta_R'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! cell sums batch wise for fixed ΔR = S_1 - R_1; for example:
      ! Σ_λσ^R = sum_PR1νS1 M^G_λ0,νS1,PR1 M^W_σR,νS1,PR1

      CALL sum_two_R_grids(bs_env%index_to_cell_3c, &
                           bs_env%index_to_cell_3c, &
                           bs_env%nimages_3c, bs_env%nimages_3c, &
                           bs_env%index_to_cell_Delta_R, &
                           bs_env%cell_to_index_Delta_R, &
                           bs_env%nimages_Delta_R)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, FMT="(T2,A,I61)") "Number of cells ΔR", bs_env%nimages_Delta_R
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_cells_Delta_R

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_parallelization_Delta_R(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_parallelization_Delta_R'

      INTEGER                                            :: handle, i_cell_Delta_R, i_task_local, &
                                                            n_tasks_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: i_cell_Delta_R_group, &
                                                            n_tensor_ops_Delta_R

      CALL timeset(routineN, handle)

      CALL compute_n_tensor_ops_Delta_R(bs_env, n_tensor_ops_Delta_R)

      CALL compute_Delta_R_dist(bs_env, n_tensor_ops_Delta_R, i_cell_Delta_R_group, n_tasks_local)

      bs_env%n_tasks_Delta_R_local = n_tasks_local

      ALLOCATE (bs_env%task_Delta_R(n_tasks_local))

      i_task_local = 0
      DO i_cell_Delta_R = 1, bs_env%nimages_Delta_R

         IF (i_cell_Delta_R_group(i_cell_Delta_R) /= bs_env%tensor_group_color) CYCLE

         i_task_local = i_task_local + 1

         bs_env%task_Delta_R(i_task_local) = i_cell_Delta_R

      END DO

      CALL timestop(handle)

   END SUBROUTINE setup_parallelization_Delta_R

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param n_tensor_ops_Delta_R ...
!> \param i_cell_Delta_R_group ...
!> \param n_tasks_local ...
! **************************************************************************************************
   SUBROUTINE compute_Delta_R_dist(bs_env, n_tensor_ops_Delta_R, i_cell_Delta_R_group, n_tasks_local)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_tensor_ops_Delta_R, &
                                                            i_cell_Delta_R_group
      INTEGER                                            :: n_tasks_local

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Delta_R_dist'

      INTEGER                                            :: handle, i_Delta_R_max_op, i_group_min, &
                                                            nimages_Delta_R, u
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_tensor_ops_Delta_R_in_group

      CALL timeset(routineN, handle)

      nimages_Delta_R = bs_env%nimages_Delta_R

      u = bs_env%unit_nr

      IF (u > 0 .AND. nimages_Delta_R < bs_env%num_tensor_groups) THEN
         WRITE (u, FMT="(T2,A,I5,A,I5,A)") "There are only ", nimages_Delta_R, &
            " tasks to work on but there are ", bs_env%num_tensor_groups, " groups."
         WRITE (u, FMT="(T2,A)") "Please reduce the number of MPI processes."
         WRITE (u, '(T2,A)') ''
      END IF

      ALLOCATE (n_tensor_ops_Delta_R_in_group(bs_env%num_tensor_groups))
      n_tensor_ops_Delta_R_in_group(:) = 0
      ALLOCATE (i_cell_Delta_R_group(nimages_Delta_R))
      i_cell_Delta_R_group(:) = -1

      n_tasks_local = 0

      DO WHILE (ANY(n_tensor_ops_Delta_R(:) .NE. 0))

         ! get largest element of n_tensor_ops_Delta_R
         i_Delta_R_max_op = MAXLOC(n_tensor_ops_Delta_R, 1)

         ! distribute i_Delta_R_max_op to tensor group which has currently the smallest load
         i_group_min = MINLOC(n_tensor_ops_Delta_R_in_group, 1)

         ! the tensor groups are 0-index based; but i_group_min is 1-index based
         i_cell_Delta_R_group(i_Delta_R_max_op) = i_group_min - 1
         n_tensor_ops_Delta_R_in_group(i_group_min) = n_tensor_ops_Delta_R_in_group(i_group_min) + &
                                                      n_tensor_ops_Delta_R(i_Delta_R_max_op)

         ! remove i_Delta_R_max_op from n_tensor_ops_Delta_R
         n_tensor_ops_Delta_R(i_Delta_R_max_op) = 0

         IF (bs_env%tensor_group_color == i_group_min - 1) n_tasks_local = n_tasks_local + 1

      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_Delta_R_dist

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param n_tensor_ops_Delta_R ...
! **************************************************************************************************
   SUBROUTINE compute_n_tensor_ops_Delta_R(bs_env, n_tensor_ops_Delta_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_tensor_ops_Delta_R

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_n_tensor_ops_Delta_R'

      INTEGER :: handle, i_cell_Delta_R, i_cell_R, i_cell_R1, i_cell_R1_minus_R, i_cell_R2, &
         i_cell_R2_m_R1, i_cell_S1, i_cell_S1_m_R1_p_R2, i_cell_S1_minus_R, i_cell_S2, &
         nimages_Delta_R
      INTEGER, DIMENSION(3) :: cell_DR, cell_m_R1, cell_R, cell_R1, cell_R1_minus_R, cell_R2, &
         cell_R2_m_R1, cell_S1, cell_S1_m_R2_p_R1, cell_S1_minus_R, cell_S1_p_S2_m_R1, cell_S2
      LOGICAL                                            :: cell_found

      CALL timeset(routineN, handle)

      nimages_Delta_R = bs_env%nimages_Delta_R

      ALLOCATE (n_tensor_ops_Delta_R(nimages_Delta_R))
      n_tensor_ops_Delta_R(:) = 0

      ! compute number of tensor operations for specific Delta_R
      DO i_cell_Delta_R = 1, nimages_Delta_R

         IF (MODULO(i_cell_Delta_R, bs_env%num_tensor_groups) /= bs_env%tensor_group_color) CYCLE

         DO i_cell_R1 = 1, bs_env%nimages_3c

            cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
            cell_DR(1:3) = bs_env%index_to_cell_Delta_R(i_cell_Delta_R, 1:3)

            ! S_1 = R_1 + ΔR (from ΔR = S_1 - R_1)
            CALL add_R(cell_R1, cell_DR, bs_env%index_to_cell_3c, cell_S1, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_S1)
            IF (.NOT. cell_found) CYCLE

            DO i_cell_R2 = 1, bs_env%nimages_scf_desymm

               cell_R2(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R2, 1:3)

               ! R_2 - R_1
               CALL add_R(cell_R2, -cell_R1, bs_env%index_to_cell_3c, cell_R2_m_R1, &
                          cell_found, bs_env%cell_to_index_3c, i_cell_R2_m_R1)
               IF (.NOT. cell_found) CYCLE

               ! S_1 - R_1 + R_2
               CALL add_R(cell_S1, cell_R2_m_R1, bs_env%index_to_cell_3c, cell_S1_m_R2_p_R1, &
                          cell_found, bs_env%cell_to_index_3c, i_cell_S1_m_R1_p_R2)
               IF (.NOT. cell_found) CYCLE

               n_tensor_ops_Delta_R(i_cell_Delta_R) = n_tensor_ops_Delta_R(i_cell_Delta_R) + 1

            END DO ! i_cell_R2

            DO i_cell_S2 = 1, bs_env%nimages_scf_desymm

               cell_S2(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_S2, 1:3)
               cell_m_R1(1:3) = -cell_R1(1:3)
               cell_S1_p_S2_m_R1(1:3) = cell_S1(1:3) + cell_S2(1:3) - cell_R1(1:3)

               CALL is_cell_in_index_to_cell(cell_m_R1, bs_env%index_to_cell_3c, cell_found)
               IF (.NOT. cell_found) CYCLE

               CALL is_cell_in_index_to_cell(cell_S1_p_S2_m_R1, bs_env%index_to_cell_3c, cell_found)
               IF (.NOT. cell_found) CYCLE

            END DO ! i_cell_S2

            DO i_cell_R = 1, bs_env%nimages_scf_desymm

               cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)

               ! R_1 - R
               CALL add_R(cell_R1, -cell_R, bs_env%index_to_cell_3c, cell_R1_minus_R, &
                          cell_found, bs_env%cell_to_index_3c, i_cell_R1_minus_R)
               IF (.NOT. cell_found) CYCLE

               ! S_1 - R
               CALL add_R(cell_S1, -cell_R, bs_env%index_to_cell_3c, cell_S1_minus_R, &
                          cell_found, bs_env%cell_to_index_3c, i_cell_S1_minus_R)
               IF (.NOT. cell_found) CYCLE

            END DO ! i_cell_R

         END DO ! i_cell_R1

      END DO ! i_cell_Delta_R

      CALL bs_env%para_env%sum(n_tensor_ops_Delta_R)

      CALL timestop(handle)

   END SUBROUTINE compute_n_tensor_ops_Delta_R

! **************************************************************************************************
!> \brief ...
!> \param cell_1 ...
!> \param cell_2 ...
!> \param index_to_cell ...
!> \param cell_1_plus_2 ...
!> \param cell_found ...
!> \param cell_to_index ...
!> \param i_cell_1_plus_2 ...
! **************************************************************************************************
   SUBROUTINE add_R(cell_1, cell_2, index_to_cell, cell_1_plus_2, cell_found, &
                    cell_to_index, i_cell_1_plus_2)

      INTEGER, DIMENSION(3)                              :: cell_1, cell_2
      INTEGER, DIMENSION(:, :)                           :: index_to_cell
      INTEGER, DIMENSION(3)                              :: cell_1_plus_2
      LOGICAL                                            :: cell_found
      INTEGER, DIMENSION(:, :, :), INTENT(IN), &
         OPTIONAL, POINTER                               :: cell_to_index
      INTEGER, INTENT(OUT), OPTIONAL                     :: i_cell_1_plus_2

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_R'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      cell_1_plus_2(1:3) = cell_1(1:3) + cell_2(1:3)

      CALL is_cell_in_index_to_cell(cell_1_plus_2, index_to_cell, cell_found)

      IF (PRESENT(i_cell_1_plus_2)) THEN
         IF (cell_found) THEN
            CPASSERT(PRESENT(cell_to_index))
            i_cell_1_plus_2 = cell_to_index(cell_1_plus_2(1), cell_1_plus_2(2), cell_1_plus_2(3))
         ELSE
            i_cell_1_plus_2 = -1000
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE add_R

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param index_to_cell ...
!> \param cell_found ...
! **************************************************************************************************
   SUBROUTINE is_cell_in_index_to_cell(cell, index_to_cell, cell_found)
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :)                           :: index_to_cell
      LOGICAL                                            :: cell_found

      CHARACTER(LEN=*), PARAMETER :: routineN = 'is_cell_in_index_to_cell'

      INTEGER                                            :: handle, i_cell, nimg
      INTEGER, DIMENSION(3)                              :: cell_i

      CALL timeset(routineN, handle)

      nimg = SIZE(index_to_cell, 1)

      cell_found = .FALSE.

      DO i_cell = 1, nimg

         cell_i(1:3) = index_to_cell(i_cell, 1:3)

         IF (cell_i(1) == cell(1) .AND. cell_i(2) == cell(2) .AND. cell_i(3) == cell(3)) THEN
            cell_found = .TRUE.
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE is_cell_in_index_to_cell

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_cfm_mo_coeff_kp_and_eigenval_scf_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cfm_mo_coeff_kp_and_eigenval_scf_kp'

      INTEGER                                            :: handle, ikp, ispin, nkp_bs_and_DOS, re_im
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_scf
      REAL(KIND=dp)                                      :: CBM, VBM
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      TYPE(cp_cfm_type)                                  :: cfm_ks, cfm_mos, cfm_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(kpoint_type), POINTER                         :: kpoints_scf
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      matrix_ks_kp=matrix_ks, &
                      matrix_s_kp=matrix_s, &
                      kpoints=kpoints_scf)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints_scf, sab_nl=sab_nl, cell_to_index=cell_to_index_scf)

      CALL cp_cfm_create(cfm_ks, bs_env%cfm_work_mo%matrix_struct)
      CALL cp_cfm_create(cfm_s, bs_env%cfm_work_mo%matrix_struct)
      CALL cp_cfm_create(cfm_mos, bs_env%cfm_work_mo%matrix_struct)

      ! nkp_bs_and_DOS contains desymmetrized k-point mesh from SCF and k-points from GW bandstructure
      nkp_bs_and_DOS = bs_env%nkp_bs_and_DOS

      ALLOCATE (bs_env%eigenval_G0W0(bs_env%n_ao, nkp_bs_and_DOS, bs_env%n_spin))
      ALLOCATE (bs_env%eigenval_HF(bs_env%n_ao, nkp_bs_and_DOS, bs_env%n_spin))
      ALLOCATE (bs_env%fm_mo_coeff_kp(nkp_bs_and_DOS, bs_env%n_spin, 2))
      ALLOCATE (bs_env%fm_ks_kp(nkp_bs_and_DOS, bs_env%n_spin, 2))
      ALLOCATE (bs_env%fm_s_kp(nkp_bs_and_DOS, bs_env%n_spin, 2))
      DO ikp = 1, nkp_bs_and_DOS
         DO ispin = 1, bs_env%n_spin
            DO re_im = 1, 2
               CALL cp_fm_create(bs_env%fm_mo_coeff_kp(ikp, ispin, re_im), &
                                 bs_env%cfm_work_mo%matrix_struct)
               CALL cp_fm_create(bs_env%fm_ks_kp(ikp, ispin, re_im), &
                                 bs_env%cfm_work_mo%matrix_struct)
               CALL cp_fm_create(bs_env%fm_s_kp(ikp, ispin, re_im), &
                                 bs_env%cfm_work_mo%matrix_struct)
            END DO
         END DO
      END DO

      DO ispin = 1, bs_env%n_spin
         DO ikp = 1, nkp_bs_and_DOS

            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, ikp)

            ! h^KS^R -> h^KS(k)
            CALL rsmat_to_kp(matrix_ks, ispin, xkp, cell_to_index_scf, sab_nl, bs_env, cfm_ks)

            ! S^R -> S(k)
            CALL rsmat_to_kp(matrix_s, 1, xkp, cell_to_index_scf, sab_nl, bs_env, cfm_s)

            ! JW comment: one might remove h^KS(k) again later
            ! we store the complex KS matrix as fm matrix because the infrastructure for fm is
            ! much nicer compared to cfm
            CALL cp_cfm_to_fm(cfm_ks, &
                              bs_env%fm_ks_kp(ikp, ispin, 1), &
                              bs_env%fm_ks_kp(ikp, ispin, 2))
            CALL cp_cfm_to_fm(cfm_s, &
                              bs_env%fm_s_kp(ikp, ispin, 1), &
                              bs_env%fm_s_kp(ikp, ispin, 2))

            ! Diagonalize KS-matrix via Rothaan-Hall equation:
            ! H^KS(k) C(k) = S(k) C(k) ε(k)
            CALL cp_cfm_geeig_canon(cfm_ks, cfm_s, cfm_mos, &
                                    bs_env%eigenval_scf(:, ikp, ispin), &
                                    bs_env%cfm_work_mo, bs_env%eps_eigval_mat_s)

            ! we store the complex MO coeff as fm matrix because the infrastructure for fm is
            ! much nicer compared to cfm
            CALL cp_cfm_to_fm(cfm_mos, &
                              bs_env%fm_mo_coeff_kp(ikp, ispin, 1), &
                              bs_env%fm_mo_coeff_kp(ikp, ispin, 2))

         END DO

         VBM = MAXVAL(bs_env%eigenval_scf(bs_env%n_occ(ispin), :, ispin))
         CBM = MINVAL(bs_env%eigenval_scf(bs_env%n_occ(ispin) + 1, :, ispin))

         bs_env%e_fermi(ispin) = 0.5_dp*(VBM + CBM)

      END DO

      CALL cp_cfm_release(cfm_ks)
      CALL cp_cfm_release(cfm_s)
      CALL cp_cfm_release(cfm_mos)

      CALL timestop(handle)

   END SUBROUTINE compute_cfm_mo_coeff_kp_and_eigenval_scf_kp

! **************************************************************************************************
!> \brief ...
!> \param matrix_ks ...
!> \param ispin ...
!> \param xkp ...
!> \param cell_to_index_scf ...
!> \param sab_nl ...
!> \param bs_env ...
!> \param cfm_ks ...
! **************************************************************************************************
   SUBROUTINE rsmat_to_kp(matrix_ks, ispin, xkp, cell_to_index_scf, sab_nl, bs_env, cfm_ks)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      INTEGER                                            :: ispin
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_scf
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_ks

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rsmat_to_kp'

      INTEGER                                            :: handle
      TYPE(dbcsr_type), POINTER                          :: cmat, nsmat, rmat

      CALL timeset(routineN, handle)

      ALLOCATE (rmat, cmat, nsmat)
      CALL dbcsr_create(rmat, template=matrix_ks(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmat, template=matrix_ks(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(nsmat, template=matrix_ks(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(rmat, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmat, sab_nl)

      CALL dbcsr_set(rmat, 0.0_dp)
      CALL dbcsr_set(cmat, 0.0_dp)
      CALL rskp_transform(rmatrix=rmat, cmatrix=cmat, rsmat=matrix_ks, ispin=ispin, &
                          xkp=xkp, cell_to_index=cell_to_index_scf, sab_nl=sab_nl)
      CALL dbcsr_desymmetrize(rmat, nsmat)
      CALL copy_dbcsr_to_fm(nsmat, bs_env%fm_work_mo(1))
      CALL dbcsr_desymmetrize(cmat, nsmat)
      CALL copy_dbcsr_to_fm(nsmat, bs_env%fm_work_mo(2))
      CALL cp_fm_to_cfm(bs_env%fm_work_mo(1), bs_env%fm_work_mo(2), cfm_ks)

      CALL dbcsr_deallocate_matrix(rmat)
      CALL dbcsr_deallocate_matrix(cmat)
      CALL dbcsr_deallocate_matrix(nsmat)

      CALL timestop(handle)

   END SUBROUTINE rsmat_to_kp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_matrices_small_cell_full_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_matrices_small_cell_full_kp'

      INTEGER                                            :: handle, i_spin, i_t, img, n_spin, &
                                                            nimages_scf, num_time_freq_points
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nimages_scf = bs_env%nimages_scf_desymm
      num_time_freq_points = bs_env%num_time_freq_points
      n_spin = bs_env%n_spin

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      ALLOCATE (bs_env%fm_G_S(nimages_scf))
      ALLOCATE (bs_env%fm_Sigma_x_R(nimages_scf))
      ALLOCATE (bs_env%fm_chi_R_t(nimages_scf, num_time_freq_points))
      ALLOCATE (bs_env%fm_MWM_R_t(nimages_scf, num_time_freq_points))
      ALLOCATE (bs_env%fm_Sigma_c_R_neg_tau(nimages_scf, num_time_freq_points, n_spin))
      ALLOCATE (bs_env%fm_Sigma_c_R_pos_tau(nimages_scf, num_time_freq_points, n_spin))
      DO img = 1, nimages_scf
         CALL cp_fm_create(bs_env%fm_G_S(img), bs_env%fm_work_mo(1)%matrix_struct)
         CALL cp_fm_create(bs_env%fm_Sigma_x_R(img), bs_env%fm_work_mo(1)%matrix_struct)
         DO i_t = 1, num_time_freq_points
            CALL cp_fm_create(bs_env%fm_chi_R_t(img, i_t), bs_env%fm_RI_RI%matrix_struct)
            CALL cp_fm_create(bs_env%fm_MWM_R_t(img, i_t), bs_env%fm_RI_RI%matrix_struct)
            CALL cp_fm_set_all(bs_env%fm_MWM_R_t(img, i_t), 0.0_dp)
            DO i_spin = 1, n_spin
               CALL cp_fm_create(bs_env%fm_Sigma_c_R_neg_tau(img, i_t, i_spin), &
                                 bs_env%fm_work_mo(1)%matrix_struct)
               CALL cp_fm_create(bs_env%fm_Sigma_c_R_pos_tau(img, i_t, i_spin), &
                                 bs_env%fm_work_mo(1)%matrix_struct)
               CALL cp_fm_set_all(bs_env%fm_Sigma_c_R_neg_tau(img, i_t, i_spin), 0.0_dp)
               CALL cp_fm_set_all(bs_env%fm_Sigma_c_R_pos_tau(img, i_t, i_spin), 0.0_dp)
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_matrices_small_cell_full_kp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE trafo_V_xc_R_to_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'trafo_V_xc_R_to_kp'

      INTEGER                                            :: handle, ikp, img, ispin, n_ao
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_scf
      TYPE(cp_cfm_type)                                  :: cfm_mo_coeff, cfm_tmp, cfm_V_xc
      TYPE(cp_fm_type)                                   :: fm_V_xc_re
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(kpoint_type), POINTER                         :: kpoints_scf
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao

      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, kpoints=kpoints_scf)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints_scf, sab_nl=sab_nl, cell_to_index=cell_to_index_scf)

      CALL cp_cfm_create(cfm_V_xc, bs_env%cfm_work_mo%matrix_struct)
      CALL cp_cfm_create(cfm_mo_coeff, bs_env%cfm_work_mo%matrix_struct)
      CALL cp_cfm_create(cfm_tmp, bs_env%cfm_work_mo%matrix_struct)
      CALL cp_fm_create(fm_V_xc_re, bs_env%cfm_work_mo%matrix_struct)

      DO img = 1, bs_env%nimages_scf
         DO ispin = 1, bs_env%n_spin
            ! JW kind of hack because the format of matrix_ks remains dubious...
            CALL dbcsr_set(matrix_ks(ispin, img)%matrix, 0.0_dp)
            CALL copy_fm_to_dbcsr(bs_env%fm_V_xc_R(img, ispin), matrix_ks(ispin, img)%matrix)
         END DO
      END DO

      ALLOCATE (bs_env%v_xc_n(n_ao, bs_env%nkp_bs_and_DOS, bs_env%n_spin))

      DO ispin = 1, bs_env%n_spin
         DO ikp = 1, bs_env%nkp_bs_and_DOS

            ! v^xc^R -> v^xc(k)  (matrix_ks stores v^xc^R, see SUBROUTINE compute_V_xc)
            CALL rsmat_to_kp(matrix_ks, ispin, bs_env%kpoints_DOS%xkp(1:3, ikp), &
                             cell_to_index_scf, sab_nl, bs_env, cfm_V_xc)

            ! get C_µn(k)
            CALL cp_fm_to_cfm(bs_env%fm_mo_coeff_kp(ikp, ispin, 1), &
                              bs_env%fm_mo_coeff_kp(ikp, ispin, 2), cfm_mo_coeff)

            ! v^xc_nm(k_i) = sum_µν C^*_µn(k_i) v^xc_µν(k_i) C_νn(k_i)
            CALL parallel_gemm('N', 'N', n_ao, n_ao, n_ao, z_one, cfm_V_xc, cfm_mo_coeff, &
                               z_zero, cfm_tmp)
            CALL parallel_gemm('C', 'N', n_ao, n_ao, n_ao, z_one, cfm_mo_coeff, cfm_tmp, &
                               z_zero, cfm_V_xc)

            ! get v^xc_nn(k_i) which is a real quantity as v^xc is Hermitian
            CALL cp_cfm_to_fm(cfm_V_xc, fm_V_xc_re)
            CALL cp_fm_get_diag(fm_V_xc_re, bs_env%v_xc_n(:, ikp, ispin))

         END DO

      END DO

      ! just rebuild the overwritten KS matrix again
      CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)

      CALL cp_cfm_release(cfm_V_xc)
      CALL cp_cfm_release(cfm_mo_coeff)
      CALL cp_cfm_release(cfm_tmp)
      CALL cp_fm_release(fm_V_xc_re)

      CALL timestop(handle)

   END SUBROUTINE trafo_V_xc_R_to_kp

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param Sigma_c_n_time ...
!> \param Sigma_c_n_freq ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE time_to_freq(bs_env, Sigma_c_n_time, Sigma_c_n_freq, ispin)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: Sigma_c_n_time, Sigma_c_n_freq
      INTEGER                                            :: ispin

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'time_to_freq'

      INTEGER                                            :: handle, i_t, j_w, n_occ
      REAL(KIND=dp)                                      :: freq_j, time_i, w_cos_ij, w_sin_ij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Sigma_c_n_cos_time, Sigma_c_n_sin_time

      CALL timeset(routineN, handle)

      ALLOCATE (Sigma_c_n_cos_time(bs_env%n_ao, bs_env%num_time_freq_points))
      ALLOCATE (Sigma_c_n_sin_time(bs_env%n_ao, bs_env%num_time_freq_points))

      Sigma_c_n_cos_time(:, :) = 0.5_dp*(Sigma_c_n_time(:, :, 1) + Sigma_c_n_time(:, :, 2))
      Sigma_c_n_sin_time(:, :) = 0.5_dp*(Sigma_c_n_time(:, :, 1) - Sigma_c_n_time(:, :, 2))

      Sigma_c_n_freq(:, :, :) = 0.0_dp

      DO i_t = 1, bs_env%num_time_freq_points

         DO j_w = 1, bs_env%num_time_freq_points

            freq_j = bs_env%imag_freq_points(j_w)
            time_i = bs_env%imag_time_points(i_t)
            ! integration weights for cosine and sine transform
            w_cos_ij = bs_env%weights_cos_t_to_w(j_w, i_t)*COS(freq_j*time_i)
            w_sin_ij = bs_env%weights_sin_t_to_w(j_w, i_t)*SIN(freq_j*time_i)

            ! 1. Re(Σ^c_nn(k_i,iω)) from cosine transform
            Sigma_c_n_freq(:, j_w, 1) = Sigma_c_n_freq(:, j_w, 1) + &
                                        w_cos_ij*Sigma_c_n_cos_time(:, i_t)

            ! 2. Im(Σ^c_nn(k_i,iω)) from sine transform
            Sigma_c_n_freq(:, j_w, 2) = Sigma_c_n_freq(:, j_w, 2) + &
                                        w_sin_ij*Sigma_c_n_sin_time(:, i_t)

         END DO

      END DO

      ! for occupied levels, we need the correlation self-energy for negative omega.
      ! Therefore, weight_sin should be computed with -omega, which results in an
      ! additional minus for the imaginary part:
      n_occ = bs_env%n_occ(ispin)
      Sigma_c_n_freq(1:n_occ, :, 2) = -Sigma_c_n_freq(1:n_occ, :, 2)

      CALL timestop(handle)

   END SUBROUTINE time_to_freq

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param Sigma_c_ikp_n_freq ...
!> \param Sigma_x_ikp_n ...
!> \param V_xc_ikp_n ...
!> \param eigenval_scf ...
!> \param ikp ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE analyt_conti_and_print(bs_env, Sigma_c_ikp_n_freq, Sigma_x_ikp_n, V_xc_ikp_n, &
                                     eigenval_scf, ikp, ispin)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: Sigma_c_ikp_n_freq
      REAL(KIND=dp), DIMENSION(:)                        :: Sigma_x_ikp_n, V_xc_ikp_n, eigenval_scf
      INTEGER                                            :: ikp, ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'analyt_conti_and_print'

      CHARACTER(len=3)                                   :: occ_vir
      CHARACTER(len=default_string_length)               :: fname
      INTEGER                                            :: handle, i_mo, ikp_for_filename, iunit, &
                                                            n_mo, nkp
      LOGICAL                                            :: is_bandstruc_kpoint, print_DOS_kpoints, &
                                                            print_ikp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dummy, Sigma_c_ikp_n_qp

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao
      ALLOCATE (dummy(n_mo), Sigma_c_ikp_n_qp(n_mo))
      Sigma_c_ikp_n_qp(:) = 0.0_dp

      DO i_mo = 1, n_mo

         ! parallelization
         IF (MODULO(i_mo, bs_env%para_env%num_pe) /= bs_env%para_env%mepos) CYCLE

         CALL continuation_pade(Sigma_c_ikp_n_qp, &
                                bs_env%imag_freq_points_fit, dummy, dummy, &
                                Sigma_c_ikp_n_freq(:, 1:bs_env%num_freq_points_fit, 1)*z_one + &
                                Sigma_c_ikp_n_freq(:, 1:bs_env%num_freq_points_fit, 2)*gaussi, &
                                Sigma_x_ikp_n(:) - V_xc_ikp_n(:), &
                                eigenval_scf(:), &
                                eigenval_scf(:), &
                                i_mo, bs_env%n_occ(ispin), bs_env%nparam_pade, &
                                bs_env%num_freq_points_fit, &
                                ri_rpa_g0w0_crossing_newton, bs_env%n_occ(ispin), &
                                0.0_dp, .TRUE., .FALSE., 1, e_fermi_ext=bs_env%e_fermi(ispin))
      END DO

      CALL bs_env%para_env%sum(Sigma_c_ikp_n_qp)

      CALL correct_obvious_fitting_fails(Sigma_c_ikp_n_qp, ispin, bs_env)

      bs_env%eigenval_G0W0(:, ikp, ispin) = eigenval_scf(:) + &
                                            Sigma_c_ikp_n_qp(:) + &
                                            Sigma_x_ikp_n(:) - &
                                            V_xc_ikp_n(:)

      bs_env%eigenval_HF(:, ikp, ispin) = eigenval_scf(:) + Sigma_x_ikp_n(:) - V_xc_ikp_n(:)

      ! only print eigenvalues of DOS k-points in case no bandstructure path has been given
      print_DOS_kpoints = (bs_env%nkp_only_bs .LE. 0)
      ! in kpoints_DOS, the last nkp_only_bs are bandstructure k-points
      is_bandstruc_kpoint = (ikp > bs_env%nkp_only_DOS)
      print_ikp = print_DOS_kpoints .OR. is_bandstruc_kpoint

      IF (bs_env%para_env%is_source() .AND. print_ikp) THEN

         IF (print_DOS_kpoints) THEN
            nkp = bs_env%nkp_only_DOS
            ikp_for_filename = ikp
         ELSE
            nkp = bs_env%nkp_only_bs
            ikp_for_filename = ikp - bs_env%nkp_only_DOS
         END IF

         CALL get_fname(fname, bs_env, ikp_for_filename, nkp, "SCF_and_G0W0", ispin=ispin)

         CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", file_action="WRITE")

         WRITE (iunit, "(A)") " "

         WRITE (iunit, "(A10,3F10.4)") "kpoint: ", bs_env%kpoints_DOS%xkp(:, ikp)
         WRITE (iunit, "(A)") " "
         WRITE (iunit, "(A5,A24,2A17,A16,A18)") "n", "ϵ_nk^DFT (eV)", "Σ^c_nk (eV)", &
            "Σ^x_nk (eV)", "v_nk^xc (eV)", "ϵ_nk^G0W0 (eV)"
         WRITE (iunit, "(A)") " "

         DO i_mo = 1, n_mo
            IF (i_mo .LE. bs_env%n_occ(ispin)) occ_vir = 'occ'
            IF (i_mo > bs_env%n_occ(ispin)) occ_vir = 'vir'
            WRITE (iunit, "(I5,3A,4F16.3,F17.3)") i_mo, ' (', occ_vir, ') ', &
               eigenval_scf(i_mo)*evolt, &
               Sigma_c_ikp_n_qp(i_mo)*evolt, &
               Sigma_x_ikp_n(i_mo)*evolt, &
               V_xc_ikp_n(i_mo)*evolt, &
               bs_env%eigenval_G0W0(i_mo, ikp, ispin)*evolt
         END DO

         CALL close_file(iunit)

      END IF

      CALL timestop(handle)

   END SUBROUTINE analyt_conti_and_print

! **************************************************************************************************
!> \brief ...
!> \param Sigma_c_ikp_n_qp ...
!> \param ispin ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE correct_obvious_fitting_fails(Sigma_c_ikp_n_qp, ispin, bs_env)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_c_ikp_n_qp
      INTEGER                                            :: ispin
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'correct_obvious_fitting_fails'

      INTEGER                                            :: handle, homo, i_mo, j_mo, &
                                                            n_levels_scissor, n_mo
      LOGICAL                                            :: is_occ, is_vir
      REAL(KIND=dp)                                      :: sum_Sigma_c

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao
      homo = bs_env%n_occ(ispin)

      DO i_mo = 1, n_mo

         ! if |𝚺^c| > 13 eV, we use a scissors shift
         IF (ABS(Sigma_c_ikp_n_qp(i_mo)) > 13.0_dp/evolt) THEN

            is_occ = (i_mo .LE. homo)
            is_vir = (i_mo > homo)

            n_levels_scissor = 0
            sum_Sigma_c = 0.0_dp

            ! compute scissor
            DO j_mo = 1, n_mo

               ! only compute scissor from other GW levels close in energy
               IF (is_occ .AND. j_mo > homo) CYCLE
               IF (is_vir .AND. j_mo .LE. homo) CYCLE
               IF (ABS(i_mo - j_mo) > 10) CYCLE
               IF (i_mo == j_mo) CYCLE

               n_levels_scissor = n_levels_scissor + 1
               sum_Sigma_c = sum_Sigma_c + Sigma_c_ikp_n_qp(j_mo)

            END DO

            ! overwrite the self-energy with scissor shift
            Sigma_c_ikp_n_qp(i_mo) = sum_Sigma_c/REAL(n_levels_scissor, KIND=dp)

         END IF

      END DO ! i_mo

      CALL timestop(handle)

   END SUBROUTINE correct_obvious_fitting_fails

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_VBM_CBM_bandgaps(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_VBM_CBM_bandgaps'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL get_VBM_CBM_bandgaps_single(bs_env%band_edges_scf, bs_env%eigenval_scf, bs_env)
      CALL get_VBM_CBM_bandgaps_single(bs_env%band_edges_G0W0, bs_env%eigenval_G0W0, bs_env)
      CALL get_VBM_CBM_bandgaps_single(bs_env%band_edges_HF, bs_env%eigenval_HF, bs_env)

      CALL timestop(handle)

   END SUBROUTINE get_VBM_CBM_bandgaps

! **************************************************************************************************
!> \brief ...
!> \param band_edges ...
!> \param ev ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_VBM_CBM_bandgaps_single(band_edges, ev, bs_env)
      TYPE(band_edges_type)                              :: band_edges
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: ev
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_VBM_CBM_bandgaps_single'

      INTEGER                                            :: handle, homo, homo_1, homo_2, ikp, &
                                                            ispin, lumo, lumo_1, lumo_2, n_mo
      REAL(KIND=dp)                                      :: E_DBG_at_ikp

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao

      band_edges%DBG = 1000.0_dp

      SELECT CASE (bs_env%n_spin)
      CASE (1)
         homo = bs_env%n_occ(1)
         lumo = homo + 1
         band_edges%VBM = MAXVAL(ev(1:homo, :, 1))
         band_edges%CBM = MINVAL(ev(homo + 1:n_mo, :, 1))
      CASE (2)
         homo_1 = bs_env%n_occ(1)
         lumo_1 = homo_1 + 1
         homo_2 = bs_env%n_occ(2)
         lumo_2 = homo_2 + 1
         band_edges%VBM = MAX(MAXVAL(ev(1:homo_1, :, 1)), MAXVAL(ev(1:homo_2, :, 2)))
         band_edges%CBM = MIN(MINVAL(ev(homo_1 + 1:n_mo, :, 1)), MINVAL(ev(homo_2 + 1:n_mo, :, 2)))
      CASE DEFAULT
         CPABORT("Error with number of spins.")
      END SELECT

      band_edges%IDBG = band_edges%CBM - band_edges%VBM

      DO ispin = 1, bs_env%n_spin

         homo = bs_env%n_occ(ispin)

         DO ikp = 1, bs_env%nkp_bs_and_DOS

            E_DBG_at_ikp = -MAXVAL(ev(1:homo, ikp, ispin)) + MINVAL(ev(homo + 1:n_mo, ikp, ispin))

            IF (E_DBG_at_ikp < band_edges%DBG) band_edges%DBG = E_DBG_at_ikp

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE get_VBM_CBM_bandgaps_single

END MODULE gw_utils

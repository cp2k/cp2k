!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
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
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_geeig_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_p_type
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_create, dbt_destroy, dbt_filter, dbt_iterator_blocks_left, &
        dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, dbt_iterator_type, &
        dbt_mp_environ_pgrid, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE input_constants,                 ONLY: do_potential_truncated,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup,&
                                              cp_libint_static_init
   USE machine,                         ONLY: m_memory,&
                                              m_walltime
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
   USE mp2_ri_2c,                       ONLY: setup_trunc_coulomb_pot_for_exchange_self_energy
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
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
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: create_and_init_bs_env_for_gw, de_init_bs_env, get_i_j_atoms, &
             kpoint_init_cell_index_simple, compute_xkp

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

      CALL read_gw_input_parameters(bs_env, bs_sec)

      CALL print_header_and_input_parameters(bs_env)

      CALL get_RI_basis_and_basis_function_indices(qs_env, bs_env)

      CALL setup_kpoints_chi_eps_W(qs_env, bs_env, bs_env%kpoints_chi_eps_W)

      CALL set_heuristic_parameters(bs_env, qs_env)

      CALL set_parallelization_parameters(qs_env, bs_env)

      CALL allocate_and_fill_matrices_and_arrays(qs_env, bs_env)

      CALL cp_libint_static_init()

      CALL create_tensors(qs_env, bs_env)

      CALL set_sparsity_parallelization_parameters(bs_env)

      CALL check_for_restart_files(qs_env, bs_env)

      CALL compute_fm_V_xc_Gamma(qs_env, bs_env)

      CALL setup_time_and_frequency_minimax_grid(bs_env)

      ! free memory in qs_env; only if one is not doing local VBM, CBM, gap because
      ! we need real-space grid operations in pw_env, task_list there
      ! Recommendation: perform GW calculation without calculating local VBM, CBM, gap
      !                 to deallocate qs_env (and safe memory). Then use GW restart files
      !                 in a subsequent calculation to calculate local VBM, CBM, gap
      IF (.NOT. bs_env%do_local_VBM_CBM_gap) THEN
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
      CALL section_vals_val_get(gw_sec, "MEMORY_PER_PROC", i_val=bs_env%input_memory_per_proc_GB)
      CALL section_vals_val_get(gw_sec, "APPROX_KP_EXTRAPOL", l_val=bs_env%approx_kp_extrapol)

      CALL timestop(handle)

   END SUBROUTINE read_gw_input_parameters

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_RI_basis_and_basis_function_indices(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_RI_basis_and_basis_function_indices'

      INTEGER                                            :: handle, iatom, ikind, &
                                                            max_AO_bf_per_atom, n_ao_test, n_atom, &
                                                            n_kind, n_RI, nsgf, u
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
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

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,A)") " "
         WRITE (UNIT=u, FMT="(T2,2A,T75,I8)") "Number of auxiliary Gaussian basis functions ", &
            "for χ, ε, W", n_RI
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_RI_basis_and_basis_function_indices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_chi_eps_W(qs_env, bs_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_chi_eps_W'

      INTEGER                                            :: handle, i_dim, nkp, nkp_extra, nkp_orig, &
                                                            u
      INTEGER, DIMENSION(3)                              :: nkp_grid, nkp_grid_extra, periodic

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      kpoints%kp_scheme = "GENERAL"

      periodic(1:3) = bs_env%periodic(1:3)

      DO i_dim = 1, 3

         CPASSERT(periodic(i_dim) == 0 .OR. periodic(i_dim) == 1)

         IF (periodic(i_dim) == 1) nkp_grid(i_dim) = 4
         IF (periodic(i_dim) == 0) nkp_grid(i_dim) = 1

         IF (periodic(i_dim) == 1) THEN
            ! only even k-meshes in periodic direction implemented
            CPASSERT(MODULO(nkp_grid(i_dim), 2) == 0)
         END IF
         IF (periodic(i_dim) == 0) THEN
            ! single k-kpoint in non-periodic direction needed
            CPASSERT(nkp_grid(i_dim) == 1)
         END IF

         IF (periodic(i_dim) == 1) nkp_grid_extra(i_dim) = nkp_grid(i_dim) + 2
         IF (periodic(i_dim) == 0) nkp_grid_extra(i_dim) = 1

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

      ! k-point weights are chosen to automatically extrapolate the k-point mesh
      kpoints%wkp(1:nkp_orig) = 1.0_dp/REAL(nkp_orig, KIND=dp) &
                                /(1.0_dp - SQRT(REAL(nkp_extra, KIND=dp)/REAL(nkp_orig, KIND=dp)))

      ! exclude nkp_extra = 0 which is the case for molecules
      IF (nkp_extra > 0) THEN
         kpoints%wkp(nkp_orig + 1:nkp) = 1.0_dp/REAL(nkp_extra, KIND=dp) &
                                         /(1.0_dp - SQRT(REAL(nkp_orig, KIND=dp) &
                                                         /REAL(nkp_extra, KIND=dp)))
      END IF

      CALL compute_xkp(kpoints%xkp, 1, nkp_orig, nkp_grid)
      CALL compute_xkp(kpoints%xkp, nkp_orig + 1, nkp, nkp_grid_extra)

      IF (bs_env%approx_kp_extrapol) THEN
         bs_env%wkp_no_extra = 1.0_dp/REAL(nkp_orig, KIND=dp)
      END IF

      CALL kpoint_init_cell_index_simple(kpoints, qs_env)

      ! heuristic parameter: how many k-points for χ, ε, and W are used simultaneously
      ! (less simultaneous k-points: less memory, but more computational effort because of
      !  recomputation of V(k))
      bs_env%nkp_chi_eps_W_batch = 4

      bs_env%num_chi_eps_W_batches = (bs_env%nkp_chi_eps_W_orig_plus_extra - 1)/ &
                                     bs_env%nkp_chi_eps_W_batch + 1

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,A)") " "
         WRITE (UNIT=u, FMT="(T2,1A,T71,3I4)") "K-point mesh 1 for χ, ε, W", nkp_grid(1:3)
         WRITE (UNIT=u, FMT="(T2,2A,T71,3I4)") "K-point mesh 2 for χ, ε, W ", &
            "(for k-point extrapolation of W)", nkp_grid_extra(1:3)
         WRITE (UNIT=u, FMT="(T2,A,T79,L)") "Approximate the k-point extrapolation", &
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
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_matrices_and_arrays(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_fill_matrices_and_arrays'

      INTEGER                                            :: handle, i_t, num_time_freq_points
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_tensor
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_RI_global
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)

      fm_struct => bs_env%fm_ks_Gamma(1)%matrix_struct

      ALLOCATE (bs_env%fm_Gocc(bs_env%num_time_freq_points))
      ALLOCATE (bs_env%fm_Gvir(bs_env%num_time_freq_points))
      DO i_t = 1, bs_env%num_time_freq_points
         CALL cp_fm_create(bs_env%fm_Gocc(i_t), fm_struct)
         CALL cp_fm_create(bs_env%fm_Gvir(i_t), fm_struct)
      END DO
      CALL cp_fm_create(bs_env%fm_h_G0W0_Gamma, fm_struct)

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

      ALLOCATE (bs_env%eigenval_G0W0(bs_env%n_ao, bs_env%nkp_DOS))

      num_time_freq_points = bs_env%num_time_freq_points

      ALLOCATE (bs_env%imag_freq_points(num_time_freq_points))
      ALLOCATE (bs_env%imag_time_points(num_time_freq_points))
      ALLOCATE (bs_env%weights_cos_t_to_w(num_time_freq_points, num_time_freq_points))
      ALLOCATE (bs_env%weights_cos_w_to_t(num_time_freq_points, num_time_freq_points))
      ALLOCATE (bs_env%weights_sin_t_to_w(num_time_freq_points, num_time_freq_points))

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

   END SUBROUTINE allocate_and_fill_matrices_and_arrays

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE create_tensors(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_tensors'

      INTEGER                                            :: handle, n_atom_step, natom, nkind, &
                                                            RI_atom
      INTEGER(int_8)                                     :: mem, non_zero_elements_sum, nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_AO, sizes_RI
      REAL(dp)                                           :: max_dist_AO_atoms, occ, occupation_sum
      TYPE(dbt_type)                                     :: t_3c_global
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_global_array
      TYPE(neighbor_list_3c_type)                        :: nl_3c_global
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! code inspired by mp2_integrals.F

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      natom=natom, nkind=nkind)

      ! set up basis
      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (bs_env%basis_set_RI(nkind), bs_env%basis_set_AO(nkind))
      CALL basis_set_list_setup(bs_env%basis_set_RI, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=bs_env%basis_set_RI)
      CALL basis_set_list_setup(bs_env%basis_set_AO, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=bs_env%basis_set_AO)

      CALL init_interaction_radii(bs_env)

      ! split blocks does not improve load balancing/efficienfy for tensor contraction, so we go
      ! with the standard atomic blocks
      CALL create_3c_t(bs_env%t_RI_AO__AO, bs_env%para_env_tensor, "(RI AO | AO)", [1, 2], [3], &
                       sizes_RI, sizes_AO, create_nl_3c=.TRUE., nl_3c=bs_env%nl_3c, qs_env=qs_env)
      CALL create_3c_t(bs_env%t_RI__AO_AO, bs_env%para_env_tensor, "(RI | AO AO)", [1], [2, 3], &
                       sizes_RI, sizes_AO)

      CALL create_2c_t(bs_env, sizes_RI, sizes_AO)

      ! check the sparsity of 3c integral tensor (µν|P); calculate maximum distance between
      ! AO atoms µ, ν where at least a single integral (µν|P) is larger than the filter threshold
      CALL create_3c_t(t_3c_global, bs_env%para_env, "(RI AO | AO)", [1, 2], [3], &
                       sizes_RI, sizes_AO, create_nl_3c=.TRUE., nl_3c=nl_3c_global, qs_env=qs_env)

      DEALLOCATE (sizes_RI, sizes_AO)

      CALL m_memory(mem)
      CALL bs_env%para_env%max(mem)

      ALLOCATE (t_3c_global_array(1, 1))
      CALL dbt_create(t_3c_global, t_3c_global_array(1, 1))

      CALL bs_env%para_env%sync()
      bs_env%t1 = m_walltime()

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
                                 int_eps=bs_env%eps_3c_int, &
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

      bs_env%t2 = m_walltime()

      bs_env%occupation_3c_int = occupation_sum
      bs_env%max_dist_AO_atoms = max_dist_AO_atoms

      CALL dbt_destroy(t_3c_global)
      CALL dbt_destroy(t_3c_global_array(1, 1))
      DEALLOCATE (t_3c_global_array)

      CALL neighbor_list_3c_destroy(nl_3c_global)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') ''
         WRITE (bs_env%unit_nr, '(T2,A,F27.1,A)') &
            'Computed 3-center integrals (µν|P), execution time', bs_env%t2 - bs_env%t1, ' s'
         WRITE (bs_env%unit_nr, '(T2,A,F50.3,A)') 'Fraction of non-zero (µν|P)', &
            occupation_sum*100, ' %'
         WRITE (bs_env%unit_nr, '(T2,A,F33.1,A)') 'Max. distance between µ,ν in non-zero (µν|P)', &
            max_dist_AO_atoms*angstrom, ' A'
         WRITE (bs_env%unit_nr, '(T2,2A,I20,A)') 'Required memory if storing all 3-center ', &
            'integrals (µν|P)', INT(REAL(non_zero_elements_sum, KIND=dp)*8.0E-9_dp), ' GB'
      END IF

      CALL timestop(handle)

   END SUBROUTINE create_tensors

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
         CALL init_interaction_radii_orb_basis(orb_basis, bs_env%eps_3c_int)

         ri_basis => bs_env%basis_set_RI(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(ri_basis, bs_env%eps_3c_int)

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
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set, para_env)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set, para_env=para_env)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_3c_int, max_dist_AO_atoms, num_cells, index_to_cell, hmat, particle_set, cell) &
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

      CALL timeset(routineN, handle)

      ! heuristic parameter to prevent out of memory
      bs_env%safety_factor_memory = 0.10_dp

      ! choose atomic range for λ ("i_atom"), ν ("j_atom") in
      ! M_λνP(iτ) = sum_µ (µν|P) G^occ_µλ(i|τ|,k=0)
      ! N_νλQ(iτ) = sum_σ (σλ|Q) G^vir_σν(i|τ|,k=0)
      ! such that M and N fit into the memory
      n_atom_per_ivl = INT(SQRT(bs_env%safety_factor_memory*bs_env%input_memory_per_proc &
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

      ! choose atomic range for µ and σ ("inner loop (IL) atom") in
      ! M_λνP(iτ) = sum_µ (µν|P) G^occ_µλ(i|τ|,k=0)
      ! N_νλQ(iτ) = sum_σ (σλ|Q) G^vir_σν(i|τ|,k=0)
      ! such that M and N fit into the memory
      n_atom_per_IL_ivl = MIN(INT(bs_env%safety_factor_memory*bs_env%input_memory_per_proc &
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

      CHARACTER(LEN=default_path_length)                 :: f_chi, f_Goc, f_Gvi, f_S_n, f_S_p, &
                                                            f_S_x, f_W_t, prefix, project_name
      INTEGER                                            :: handle, i_t_or_w, num_time_freq_points
      LOGICAL                                            :: chi_exists, G_occ_exists, G_vir_exists, &
                                                            Sigma_neg_time_exists, &
                                                            Sigma_pos_time_exists, W_time_exists
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input, print_key

      CALL timeset(routineN, handle)

      num_time_freq_points = bs_env%num_time_freq_points

      ALLOCATE (bs_env%read_chi_and_G_occ_vir(num_time_freq_points))
      ALLOCATE (bs_env%calc_chi_and_G_occ_vir(num_time_freq_points))
      ALLOCATE (bs_env%Sigma_c_exists(num_time_freq_points))

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
            WRITE (f_chi, '(3A,I1,A)') TRIM(prefix), bs_env%chi_name, "_00", i_t_or_w, ".matrix"
            WRITE (f_Goc, '(3A,I1,A)') TRIM(prefix), bs_env%Gocc_name, "_00", i_t_or_w, ".matrix"
            WRITE (f_Gvi, '(3A,I1,A)') TRIM(prefix), bs_env%Gvir_name, "_00", i_t_or_w, ".matrix"
            WRITE (f_W_t, '(3A,I1,A)') TRIM(prefix), bs_env%W_time_name, "_00", i_t_or_w, ".matrix"
            WRITE (f_S_p, '(3A,I1,A)') TRIM(prefix), bs_env%Sigma_p_name, "_00", i_t_or_w, ".matrix"
            WRITE (f_S_n, '(3A,I1,A)') TRIM(prefix), bs_env%Sigma_n_name, "_00", i_t_or_w, ".matrix"
         ELSE IF (i_t_or_w < 100) THEN
            WRITE (f_chi, '(3A,I2,A)') TRIM(prefix), bs_env%chi_name, "_0", i_t_or_w, ".matrix"
            WRITE (f_Goc, '(3A,I2,A)') TRIM(prefix), bs_env%Gocc_name, "_0", i_t_or_w, ".matrix"
            WRITE (f_Gvi, '(3A,I2,A)') TRIM(prefix), bs_env%Gvir_name, "_0", i_t_or_w, ".matrix"
            WRITE (f_W_t, '(3A,I2,A)') TRIM(prefix), bs_env%W_time_name, "_0", i_t_or_w, ".matrix"
            WRITE (f_S_p, '(3A,I2,A)') TRIM(prefix), bs_env%Sigma_p_name, "_0", i_t_or_w, ".matrix"
            WRITE (f_S_n, '(3A,I2,A)') TRIM(prefix), bs_env%Sigma_n_name, "_0", i_t_or_w, ".matrix"
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         INQUIRE (file=TRIM(f_chi), exist=chi_exists)
         INQUIRE (file=TRIM(f_Goc), exist=G_occ_exists)
         INQUIRE (file=TRIM(f_Gvi), exist=G_vir_exists)
         INQUIRE (file=TRIM(f_W_t), exist=W_time_exists)
         INQUIRE (file=TRIM(f_S_p), exist=Sigma_pos_time_exists)
         INQUIRE (file=TRIM(f_S_n), exist=Sigma_neg_time_exists)

         bs_env%read_chi_and_G_occ_vir(i_t_or_w) = chi_exists .AND. G_occ_exists .AND. G_vir_exists
         bs_env%calc_chi_and_G_occ_vir(i_t_or_w) = .NOT. bs_env%read_chi_and_G_occ_vir(i_t_or_w)

         bs_env%all_W_exist = bs_env%all_W_exist .AND. W_time_exists

         bs_env%Sigma_c_exists(i_t_or_w) = Sigma_pos_time_exists .AND. Sigma_neg_time_exists

         IF (bs_env%Sigma_c_exists(i_t_or_w)) THEN
            bs_env%read_chi_and_G_occ_vir(i_t_or_w) = .FALSE.
            bs_env%calc_chi_and_G_occ_vir(i_t_or_w) = .FALSE.
         END IF

      END DO

      WRITE (f_S_x, '(3A)') TRIM(prefix), bs_env%Sigma_x_name, "_000.matrix"
      INQUIRE (file=TRIM(f_S_x), exist=bs_env%Sigma_x_exists)

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
      ! use all processors for the group (in principle, number could be changed, but performance
      ! seems to be best for a single group with all MPI processes per group)
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

      bs_env%input_memory_per_proc = INT(bs_env%input_memory_per_proc_GB*1.0E9_dp, KIND=int_8)

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,I47)') 'Group size for tensor operations', bs_env%group_size_tensor
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,I37,A)') 'Input: Available memory per MPI process', &
            bs_env%input_memory_per_proc_GB, ' GB'
         WRITE (u, '(T2,A,F35.1,A)') 'Used memory per MPI process before GW run', &
            REAL(mem, KIND=dp)/1.E9_dp, ' GB'
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

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! use the same threshold for computing 3-center integrals (µν|P) as for filtering
      ! tensor operations
      bs_env%eps_3c_int = bs_env%eps_filter

      ! Determines number of cells used for summing the cells R in the Coulomb matrix,
      ! V_PQ(k) = \sum_R <P,cell=0 | 1/r | Q,cell=R>. SIZE_LATTICE_SUM_V 4 gives
      ! excellent convergence
      bs_env%size_lattice_sum_V = 4

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

      ! only use interval ω in [0, 27.211 eV] (1 Hartree = 27.211 eV) for virt, and ω in
      ! [-27.211 eV, 0] for occ for use in analytic continuation of
      ! self-energy Σ^c_n(iω,k) -> Σ^c_n(ϵ,k)
      bs_env%freq_max_fit = 1.0_dp

      ! use a 16-parameter Padé fit
      bs_env%nparam_pade = 16

      ! minimum block size for tensor operations, taken from MP2/RPA input
      bs_env%min_block_size = 5

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
      CALL setup_trunc_coulomb_pot_for_exchange_self_energy(qs_env, bs_env%trunc_coulomb, &
                                                            rel_cutoff_trunc_coulomb_ri_x=0.5_dp)

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
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE diagonalize_ks_matrix(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_ks_matrix'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! use work matrices because the matrices are overwritten in cp_fm_geeig_canon
      CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(1), bs_env%fm_work_mo(1))
      CALL cp_fm_to_fm(bs_env%fm_s_Gamma, bs_env%fm_work_mo(2))

      ! diagonalize the Kohn-Sham matrix to obtain MO coefficients and SCF eigenvalues
      ! (at the Gamma-point)
      CALL cp_fm_geeig_canon(bs_env%fm_work_mo(1), &
                             bs_env%fm_work_mo(2), &
                             bs_env%fm_mo_coeff_Gamma(1), &
                             bs_env%eigenval_scf_Gamma, &
                             bs_env%fm_work_mo(3), &
                             bs_env%eps_eigval_mat_s)

      bs_env%e_fermi(1) = 0.5_dp*(bs_env%eigenval_scf_Gamma(bs_env%n_occ(1)) + &
                                  bs_env%eigenval_scf_Gamma(bs_env%n_occ(1) + 1))

      CALL timestop(handle)

   END SUBROUTINE diagonalize_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_fm_V_xc_Gamma(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_V_xc_Gamma'

      INTEGER                                            :: handle, ispin, myfun, nimages
      REAL(KIND=dp)                                      :: energy_ex, energy_exc, energy_total
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_ks_without_v_xc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input, energy=energy, dft_control=dft_control)

      ! previously, dft_control%nimages set to # neighbor cells, revert for Γ-only KS matrix
      nimages = dft_control%nimages
      dft_control%nimages = 1

      ! we need to reset XC functional, therefore, get XC input
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=xc_none)

      ! save the energy before the energy gets updated
      energy_total = energy%total
      energy_exc = energy%exc
      energy_ex = energy%ex

      NULLIFY (mat_ks_without_v_xc)
      CALL dbcsr_allocate_matrix_set(mat_ks_without_v_xc, bs_env%n_spin)

      DO ispin = 1, bs_env%n_spin
         ALLOCATE (mat_ks_without_v_xc(ispin)%matrix)
         CALL dbcsr_create(mat_ks_without_v_xc(ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
      END DO

      ! calculate KS-matrix without XC
      CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                        ext_ks_matrix=mat_ks_without_v_xc)

      ! transfer dbcsr matrix to fm
      CALL cp_fm_create(bs_env%fm_V_xc_Gamma(1), bs_env%fm_s_Gamma%matrix_struct)
      CALL copy_dbcsr_to_fm(mat_ks_without_v_xc(1)%matrix, bs_env%fm_V_xc_Gamma(1))

      ! finally compute the xc potential in the AO basis
      CALL cp_fm_scale_and_add(alpha=-1.0_dp, matrix_a=bs_env%fm_V_xc_Gamma(1), &
                               beta=1.0_dp, matrix_b=bs_env%fm_ks_Gamma(1))

      ! set back the energy
      energy%total = energy_total
      energy%exc = energy_exc
      energy%ex = energy_ex

      ! set back nimages
      dft_control%nimages = nimages

      ! set the DFT functional and HF fraction back
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)

      CALL dbcsr_deallocate_matrix_set(mat_ks_without_v_xc)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_V_xc_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_time_and_frequency_minimax_grid(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_time_and_frequency_minimax_grid'

      INTEGER                                            :: handle, homo, i_w, ierr, j_w, n_mo, &
                                                            num_time_freq_points, u
      REAL(KIND=dp)                                      :: E_max, E_min, E_range, max_error_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: points_and_weights

      CALL timeset(routineN, handle)

      homo = bs_env%n_occ(1)
      n_mo = bs_env%n_ao
      num_time_freq_points = bs_env%num_time_freq_points

      ! minimum and maximum difference between eigenvalues of unoccupied and an occupied MOs
      E_min = bs_env%eigenval_scf_Gamma(homo + 1) - bs_env%eigenval_scf_Gamma(homo)
      E_max = bs_env%eigenval_scf_Gamma(n_mo) - bs_env%eigenval_scf_Gamma(1)

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
      DO i_w = 1, bs_env%num_time_freq_points
         IF (bs_env%imag_freq_points(i_w) < bs_env%freq_max_fit) THEN
            bs_env%num_freq_points_fit = bs_env%num_freq_points_fit + 1
         END IF
      END DO

      ! iω values for the analytic continuation Σ^c_n(iω,k) -> Σ^c_n(ϵ,k)
      ALLOCATE (bs_env%imag_freq_points_fit(bs_env%num_freq_points_fit))
      j_w = 0
      DO i_w = 1, bs_env%num_time_freq_points
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

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,F44.2)') 'SCF direct band gap at Γ-point (eV)', E_min*evolt
         WRITE (u, '(T2,A,F42.2)') 'Max. SCF eigval diff. at Γ-point (eV)', E_max*evolt
         WRITE (u, '(T2,A,F55.2)') 'E-Range for minimax grid', E_range
         WRITE (u, '(T2,A,I27)') 'Number of Padé parameters for analytic continuation:', &
            bs_env%nparam_pade
         WRITE (u, '(T2,A)') ''
      END IF
      CALL timestop(handle)

   END SUBROUTINE setup_time_and_frequency_minimax_grid

END MODULE gw_utils

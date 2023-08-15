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
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_diag,                      ONLY: cp_fm_geeig_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_to_fm
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_type_no_symmetry
   USE dbt_api,                         ONLY: dbt_mp_environ_pgrid,&
                                              dbt_pgrid_create,&
                                              dbt_pgrid_destroy,&
                                              dbt_pgrid_type,&
                                              dbt_type
   USE gw_types,                        ONLY: gw_type
   USE input_constants,                 ONLY: do_potential_truncated
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup,&
                                              cp_libint_static_init
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
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tensors,                      ONLY: build_3c_neighbor_lists,&
                                              neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              distribution_3d_create,&
                                              distribution_3d_type
   USE scf_control_types,               ONLY: scf_control_type
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: create_and_init_gw_env, de_init_gw_env

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
!> \param gw_sec ...
! **************************************************************************************************
   SUBROUTINE create_and_init_gw_env(qs_env, gw_env, gw_sec)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(section_vals_type), POINTER                   :: gw_sec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_and_init_gw_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (gw_env)

      CALL read_input_parameters(gw_env, gw_sec)

      CALL print_header_and_input_parameters(gw_env)

      CALL get_parameters_from_qs_env(qs_env, gw_env)

      CALL setup_kpoints_chi_eps_W(qs_env, gw_env, gw_env%kpoints_chi_eps_W)

      CALL setup_kpoints_Sigma(qs_env, gw_env, gw_env%kpoints_Sigma)

      CALL set_heuristic_parameters(gw_env, qs_env)

      CALL set_parallelization_parameters(qs_env, gw_env)

      CALL allocate_and_fill_matrices_and_arrays(qs_env, gw_env)

      CALL cp_libint_static_init()

      CALL create_tensors(qs_env, gw_env)

      CALL diagonalize_ks_matrix(gw_env)

      CALL setup_time_and_frequency_minimax_grid(gw_env)

      CALL timestop(handle)

   END SUBROUTINE create_and_init_gw_env

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE de_init_gw_env(gw_env)
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'de_init_gw_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      ! deallocate quantities here which:
      ! 1. cannot be deallocated in gw_env_release due to circular dependencies
      ! 2. consume a lot of memory and should not be kept until the quantity is
      !    deallocated in gw_env_release

      IF (ASSOCIATED(gw_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(gw_env%nl_3c)

      CALL cp_libint_static_cleanup()

      CALL timestop(handle)

   END SUBROUTINE de_init_gw_env

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param gw_sec ...
! **************************************************************************************************
   SUBROUTINE read_input_parameters(gw_env, gw_sec)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(section_vals_type), POINTER                   :: gw_sec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_input_parameters'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(gw_sec, "NUM_TIME_FREQ_POINTS", i_val=gw_env%num_time_freq_points)
      CALL section_vals_val_get(gw_sec, "EPS_FILTER", r_val=gw_env%eps_filter)
      CALL section_vals_val_get(gw_sec, "EPS_3C_INT", r_val=gw_env%eps_3c_int)
      CALL section_vals_val_get(gw_sec, "KPOINTS_CHI_EPS_W", i_vals=gw_env%nkp_grid_chi_eps_W_input)
      CALL section_vals_val_get(gw_sec, "KPOINTS_SIGMA", i_vals=gw_env%nkp_grid_Sigma_input)
      CALL section_vals_val_get(gw_sec, "SIZE_LATTICE_SUM_V", i_val=gw_env%size_lattice_sum_V)

      CALL timestop(handle)

   END SUBROUTINE read_input_parameters

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE get_parameters_from_qs_env(qs_env, gw_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_parameters_from_qs_env'

      INTEGER                                            :: handle, homo, iatom, ikind, n_ao, &
                                                            n_atom, n_kind, n_RI, nsgf, u
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      scf_control=scf_control, &
                      mos=mos)

      gw_env%n_spin = dft_control%nspins
      IF (gw_env%n_spin == 1) gw_env%spin_degeneracy = 2.0_dp
      IF (gw_env%n_spin == 2) gw_env%spin_degeneracy = 1.0_dp

      ! so far, only closed-shell
      CPASSERT(gw_env%n_spin == 1)

      CALL get_mo_set(mo_set=mos(1), nao=n_ao, homo=homo)
      gw_env%n_ao = n_ao
      gw_env%n_occ(1:2) = homo
      gw_env%n_vir(1:2) = n_ao - homo

      IF (gw_env%n_spin == 2) THEN
         CALL get_mo_set(mo_set=mos(2), homo=homo)
         gw_env%n_occ(2) = homo
         gw_env%n_vir(2) = n_ao - homo
      END IF

      gw_env%eps_eigval_mat_s = scf_control%eps_eigval

      ! determine RI basis set size
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      n_kind = SIZE(qs_kind_set)
      n_atom = SIZE(particle_set)
      gw_env%n_atom = n_atom

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO ikind = 1, n_kind
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, &
                          basis_type="RI_AUX")
         CPASSERT(ASSOCIATED(basis_set_a))
      END DO

      n_RI = 0
      DO iatom = 1, n_atom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="RI_AUX")
         n_RI = n_RI + nsgf
      END DO

      gw_env%n_RI = n_RI

      u = gw_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,A)") " "
         WRITE (UNIT=u, FMT="(T2,2A,T77,I4)") "Number of occupied molecular orbitals (MOs) ", &
            "= Number of occupied bands", homo
         WRITE (UNIT=u, FMT="(T2,2A,T77,I4)") "Number of unoccupied (= virtual) MOs ", &
            "= Number of unoccupied bands", n_ao - homo
         WRITE (UNIT=u, FMT="(T2,A,T77,I4)") "Number of Gaussian basis functions for MOs", n_ao
         WRITE (UNIT=u, FMT="(T2,2A,T79,I4)") "Number of auxiliary Gaussian basis functions ", &
            "for χ, ε, W", n_RI
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_parameters_from_qs_env

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_chi_eps_W(qs_env, gw_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_chi_eps_W'

      INTEGER                                            :: handle, i_dim, nkp, nkp_extra, nkp_orig, &
                                                            u
      INTEGER, DIMENSION(3)                              :: nkp_grid, nkp_grid_extra, periodic
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      kpoints%kp_scheme = "GENERAL"

      NULLIFY (cell)
      CALL get_qs_env(qs_env=qs_env, cell=cell)
      CALL get_cell(cell=cell, periodic=periodic)

      DO i_dim = 1, 3

         CPASSERT(periodic(i_dim) == 0 .OR. periodic(i_dim) == 1)

         IF (gw_env%nkp_grid_chi_eps_W_input(i_dim) < 0) THEN
            IF (periodic(i_dim) == 1) nkp_grid(i_dim) = 4
            IF (periodic(i_dim) == 0) nkp_grid(i_dim) = 1
         ELSE
            nkp_grid(i_dim) = gw_env%nkp_grid_chi_eps_W_input(i_dim)
         END IF

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

      gw_env%nkp_grid_chi_eps_W_orig(1:3) = nkp_grid(1:3)
      gw_env%nkp_grid_chi_eps_W_extra(1:3) = nkp_grid_extra(1:3)
      gw_env%nkp_chi_eps_W_orig = nkp_orig
      gw_env%nkp_chi_eps_W_extra = nkp_extra

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))

      ! k-point weights are chosen to automatically extrapolate the k-point mesh
      kpoints%wkp(1:nkp_orig) = 1.0_dp/REAL(nkp_orig, KIND=dp) &
                                /(1.0_dp - SQRT(REAL(nkp_extra, KIND=dp)/REAL(nkp_orig, KIND=dp)))
      kpoints%wkp(nkp_orig + 1:nkp) = 1.0_dp/REAL(nkp_extra, KIND=dp) &
                                      /(1.0_dp - SQRT(REAL(nkp_orig, KIND=dp) &
                                                      /REAL(nkp_extra, KIND=dp)))

      CALL compute_xkp(kpoints%xkp, 1, nkp_orig, nkp_grid)
      CALL compute_xkp(kpoints%xkp, nkp_orig + 1, nkp, nkp_grid_extra)

      CALL kpoint_init_cell_index_simple(kpoints, qs_env)

      ! JW 2do: Parallelization over k-points (diagonalization of χ(iω,k) and further steps)
      gw_env%ikp_local_start = 1
      gw_env%ikp_local_end = kpoints%nkp

      u = gw_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,A)") " "
         WRITE (UNIT=u, FMT="(T2,1A,T71,3I4)") "K-point mesh 1 for χ, ε, W", nkp_grid(1:3)
         WRITE (UNIT=u, FMT="(T2,2A,T71,3I4)") "K-point mesh 2 for χ, ε, W ", &
            "(for k-point extrapolation of W)", nkp_grid_extra(1:3)
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_chi_eps_W

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_Sigma(qs_env, gw_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_Sigma'

      INTEGER                                            :: handle, i_dim, nkp
      INTEGER, DIMENSION(3)                              :: nkp_grid, periodic
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      kpoints%kp_scheme = "GENERAL"

      NULLIFY (cell)
      CALL get_qs_env(qs_env=qs_env, cell=cell)
      CALL get_cell(cell=cell, periodic=periodic)

      DO i_dim = 1, 3

         CPASSERT(periodic(i_dim) == 0 .OR. periodic(i_dim) == 1)

         IF (gw_env%nkp_grid_Sigma_input(i_dim) < 0) THEN
            IF (periodic(i_dim) == 1) nkp_grid(i_dim) = 8
            IF (periodic(i_dim) == 0) nkp_grid(i_dim) = 1
         ELSE
            nkp_grid(i_dim) = gw_env%nkp_grid_Sigma_input(i_dim)
         END IF

      END DO

      nkp = nkp_grid(1)*nkp_grid(2)*nkp_grid(3)

      kpoints%nkp_grid(1:3) = nkp_grid(1:3)
      kpoints%nkp = nkp

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))

      ! k-point weights are chosen to automatically extrapolate the k-point mesh
      kpoints%wkp(1:nkp) = 1.0_dp/REAL(nkp, KIND=dp)

      CALL compute_xkp(kpoints%xkp, 1, nkp, nkp_grid)

      CALL kpoint_init_cell_index_simple(kpoints, qs_env)

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_Sigma

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
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_matrices_and_arrays(qs_env, gw_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_fill_matrices_and_arrays'

      INTEGER                                            :: handle, i_t, i_work, ispin, &
                                                            num_time_freq_points
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_kp_diag, &
                                                            blacs_env_tensor
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=gw_env%n_ao, &
                               ncol_global=gw_env%n_ao, para_env=para_env)

      DO i_work = 1, SIZE(gw_env%fm_work_mo)
         CALL cp_fm_create(gw_env%fm_work_mo(i_work), fm_struct)
      END DO

      CALL cp_fm_create(gw_env%fm_s_Gamma, fm_struct)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, gw_env%fm_s_Gamma)

      DO ispin = 1, gw_env%n_spin
         CALL cp_fm_create(gw_env%fm_ks_Gamma(ispin), fm_struct)
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, gw_env%fm_ks_Gamma(ispin))
         CALL cp_fm_create(gw_env%fm_mo_coeff_Gamma(ispin), fm_struct)
      END DO

      ALLOCATE (gw_env%fm_Gocc(gw_env%num_time_freq_points))
      ALLOCATE (gw_env%fm_Gvir(gw_env%num_time_freq_points))
      DO i_t = 1, gw_env%num_time_freq_points
         CALL cp_fm_create(gw_env%fm_Gocc(i_t), fm_struct)
         CALL cp_fm_create(gw_env%fm_Gvir(i_t), fm_struct)
      END DO

      CALL cp_fm_struct_release(fm_struct)

! JW 2del THIS DOES NOT WORK UNFORUNATELY BECAUSE IT GIVES A CORRUPTED LIST AT THE END OF CP2K RUN
!      NULLIFY(fm_struct_RI_global)
!      CALL cp_fm_struct_create(fm_struct_RI_global, context=blacs_env, nrow_global=gw_env%n_RI, &
!                               ncol_global=gw_env%n_RI, para_env=para_env)
!
!!      CALL cp_fm_create(gw_env%fm_RI_RI, fm_struct_RI_global)
!!      CALL cp_cfm_create(gw_env%cfm_RI_RI, fm_struct_RI_global)
!
!!      ALLOCATE (gw_env%fm_V_kp(gw_env%kpoints_chi_eps_W%nkp, 2))
!!      DO ikp = 1, gw_env%kpoints_chi_eps_W%nkp
!!       DO re_im = 1,2
!!         CALL cp_fm_create(gw_env%fm_V_kp(ikp, re_im), fm_struct_RI_global)
!!        END DO
!!      END DO
!
!      CALL cp_fm_struct_release(fm_struct_RI_global)

      NULLIFY (gw_env%mat_ao_ao%matrix)
      ALLOCATE (gw_env%mat_ao_ao%matrix)
      CALL dbcsr_create(gw_env%mat_ao_ao%matrix, template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ALLOCATE (gw_env%eigenvalues_scf_Gamma(gw_env%n_ao))
      ALLOCATE (gw_env%eigenvalues_scf(gw_env%n_ao, gw_env%kpoints_Sigma%nkp))
      ALLOCATE (gw_env%eigenvalues_G0W0(gw_env%n_ao, gw_env%kpoints_Sigma%nkp))

      num_time_freq_points = gw_env%num_time_freq_points

      ALLOCATE (gw_env%imag_freq_points(num_time_freq_points))
      ALLOCATE (gw_env%imag_time_points(num_time_freq_points))
      ALLOCATE (gw_env%weights_cos_t_to_w(num_time_freq_points, num_time_freq_points))
      ALLOCATE (gw_env%weights_cos_w_to_t(num_time_freq_points, num_time_freq_points))
      ALLOCATE (gw_env%weights_sin_t_to_w(num_time_freq_points, num_time_freq_points))

      ! create blacs_env for subgroups of tensor operations
      NULLIFY (blacs_env_tensor)
      CALL cp_blacs_env_create(blacs_env=blacs_env_tensor, para_env=gw_env%para_env_tensor)

      NULLIFY (blacs_env_kp_diag)
      CALL cp_blacs_env_create(blacs_env=blacs_env_kp_diag, para_env=gw_env%para_env_kp_diag)

!! JW 2DO later cfm matrices in kp diag subgroups
!      NULLIFY (fm_struct_RI_kp_diag)
!      CALL cp_fm_struct_create(fm_struct_RI_kp_diag, context=blacs_env_kp_diag, &
!                               nrow_global=gw_env%n_RI, ncol_global=gw_env%n_RI, &
!                               para_env=gw_env%para_env_tensor)
!
!      PRINT *, 'check 4'
!
!      CALL cp_fm_struct_release(fm_struct_RI_kp_diag)

      ! JW: allocate dbcsr matrices in the tensor subgroup; actually, one only needs a small
      ! subset of blocks in the tensor subgroup, however, all atomic blocks are allocated.
      ! Later in case of poor performance/memory footprint, one might think of creating a
      ! dbcsr matrix with only the blocks that are needed in the tensor subgroup
      CALL create_mat_munu(gw_env%mat_ao_ao_tensor, &
                           qs_env, gw_env%eps_atom_grid_2d_mat, &
                           blacs_env_tensor, do_ri_aux_basis=.FALSE.)

      CALL create_mat_munu(gw_env%mat_RI_RI_tensor, &
                           qs_env, gw_env%eps_atom_grid_2d_mat, &
                           blacs_env_tensor, do_ri_aux_basis=.TRUE.)

      CALL create_mat_munu(gw_env%mat_RI_RI, &
                           qs_env, gw_env%eps_atom_grid_2d_mat, &
                           blacs_env, do_ri_aux_basis=.TRUE.)

      CALL cp_blacs_env_release(blacs_env_tensor)
      CALL cp_blacs_env_release(blacs_env_kp_diag)

      CALL timestop(handle)

   END SUBROUTINE allocate_and_fill_matrices_and_arrays

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE create_tensors(qs_env, gw_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_tensors'

      INTEGER                                            :: handle, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_AO, sizes_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! JW: code inspired by mp2_integrals.F

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      natom=natom, nkind=nkind)

      ! set up basis
      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (gw_env%basis_set_RI(nkind), gw_env%basis_set_AO(nkind))
      CALL basis_set_list_setup(gw_env%basis_set_RI, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=gw_env%basis_set_RI)
      CALL basis_set_list_setup(gw_env%basis_set_AO, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=gw_env%basis_set_AO)

      CALL init_interaction_radii(gw_env)

      ! split blocks to improve load balancing for tensor contraction
      ! JW TODO: that needs to be checked whether performance can be really improved by the split
!      min_block_size = gw_env%min_block_size
!      CALL pgf_block_sizes(atomic_kind_set, gw_env%basis_set_AO, min_block_size, sizes_AO_split)
!      CALL pgf_block_sizes(atomic_kind_set, gw_env%basis_set_RI, min_block_size, sizes_RI_split)

      CALL create_3c_t(gw_env%t_RI_AO__AO, gw_env%para_env_tensor, "(RI AO | AO)", [1, 2], [3], &
                       sizes_RI, sizes_AO, create_dist=.TRUE., qs_env=qs_env)
      CALL create_3c_t(gw_env%t_RI__AO_AO, gw_env%para_env_tensor, "(RI | AO AO)", [1], [2, 3], &
                       sizes_RI, sizes_AO)

      CALL create_2c_t(gw_env, sizes_RI, sizes_AO)

!      DO ibasis = 1, SIZE(gw_env%basis_set_AO)
!
!         orb_basis => gw_env%basis_set_AO(ibasis)%gto_basis_set
!         CALL init_interaction_radii_orb_basis(orb_basis, gw_env%eps_3c_int)
!
!         ri_basis => gw_env%basis_set_RI(ibasis)%gto_basis_set
!         CALL init_interaction_radii_orb_basis(ri_basis, gw_env%eps_3c_int)
!
!      END DO

      DEALLOCATE (sizes_RI, sizes_AO)

      CALL timestop(handle)

   END SUBROUTINE create_tensors

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param sizes_RI ...
!> \param sizes_AO ...
! **************************************************************************************************
   SUBROUTINE create_2c_t(gw_env, sizes_RI, sizes_AO)
      TYPE(gw_type), POINTER                             :: gw_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI, sizes_AO

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_2c_t'

      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      TYPE(dbt_pgrid_type)                               :: pgrid_2d

      CALL timeset(routineN, handle)

      ! JW: inspired from rpa_im_time.F / hfx_types.F

      pdims_2d = 0
      CALL dbt_pgrid_create(gw_env%para_env_tensor, pdims_2d, pgrid_2d)

      CALL create_2c_tensor(gw_env%t_G, dist_1, dist_2, pgrid_2d, sizes_AO, sizes_AO, &
                            name="(AO | AO)")
      DEALLOCATE (dist_1, dist_2)
      CALL create_2c_tensor(gw_env%t_chi, dist_1, dist_2, pgrid_2d, sizes_RI, sizes_RI, &
                            name="(RI | RI)")
      DEALLOCATE (dist_1, dist_2)
      CALL create_2c_tensor(gw_env%t_W, dist_1, dist_2, pgrid_2d, sizes_RI, sizes_RI, &
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
!> \param create_dist ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE create_3c_t(tensor, para_env, tensor_name, map1, map2, sizes_RI, sizes_AO, &
                          create_dist, qs_env)
      TYPE(dbt_type)                                     :: tensor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=12)                                  :: tensor_name
      INTEGER, DIMENSION(:)                              :: map1, map2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI, sizes_AO
      LOGICAL, OPTIONAL                                  :: create_dist
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_3c_t'

      INTEGER                                            :: handle, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_AO_1, dist_AO_2, dist_RI
      INTEGER, DIMENSION(3)                              :: pcoord, pdims, pdims_3d
      LOGICAL                                            :: my_create_dist
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

      IF (PRESENT(create_dist)) THEN
         my_create_dist = create_dist
      ELSE
         my_create_dist = .FALSE.
      END IF

      IF (my_create_dist) THEN
         CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set)
         CALL dbt_mp_environ_pgrid(pgrid_3d, pdims, pcoord)
         CALL mp_comm_t3c_2%create(pgrid_3d%mp_comm_2d, 3, pdims)
         CALL distribution_3d_create(dist_3d, dist_RI, dist_AO_1, dist_AO_2, &
                                     nkind, particle_set, mp_comm_t3c_2, own_comm=.TRUE.)

         CALL build_3c_neighbor_lists(qs_env%gw_env%nl_3c, &
                                      qs_env%gw_env%basis_set_RI, &
                                      qs_env%gw_env%basis_set_AO, &
                                      qs_env%gw_env%basis_set_AO, &
                                      dist_3d, qs_env%gw_env%ri_metric, &
                                      "GW_3c_nl", qs_env, own_dist=.TRUE.)
      END IF

      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
      CALL dbt_pgrid_destroy(pgrid_3d)

      CALL timestop(handle)

   END SUBROUTINE create_3c_t

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE init_interaction_radii(gw_env)
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_interaction_radii'

      INTEGER                                            :: handle, ibasis
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis

      CALL timeset(routineN, handle)

      DO ibasis = 1, SIZE(gw_env%basis_set_AO)

         orb_basis => gw_env%basis_set_AO(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(orb_basis, gw_env%eps_3c_int)

         ri_basis => gw_env%basis_set_RI(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(ri_basis, gw_env%eps_3c_int)

      END DO

      CALL timestop(handle)

   END SUBROUTINE init_interaction_radii

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE set_parallelization_parameters(qs_env, gw_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_parallelization_parameters'

      INTEGER :: color_sub, gcd_max, handle, i_atoms_per_group, i_group, ipcol, ipcol_loop, ipe, &
         iprow, iprow_loop, j_atoms_per_group, jpe, npcol, nprow, num_pe, num_t_groups
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env)

      num_pe = para_env%num_pe
      ! JW NEEDS TO BE OPTIMIZED, POTENTIALLY FROM AVAILABLE MEMORY
!      gw_env%tensor_group_size = MAX(num_pe/100, 1)
!      gw_env%tensor_group_size = num_pe
      gw_env%tensor_group_size = 1

      ! JW NEEDS TO BE OPTIMIZED, POTENTIALLY FROM AVAILABLE MEMORY
      gw_env%kp_diag_group_size = num_pe

      ! all MPI processes in the same subgroup; means that gw_env%para_env
      ! is identical to para_env from qs_env
      color_sub = 0
      ALLOCATE (gw_env%para_env)
      CALL gw_env%para_env%from_split(para_env, color_sub)

      color_sub = para_env%mepos/gw_env%tensor_group_size
      ! if the tensor group size does not divide the number of MPI processes
      ! without rest, make sure that there is no "small" last group
      IF (MODULO(num_pe, gw_env%tensor_group_size) .NE. 0 .AND. &
          color_sub == num_pe/gw_env%tensor_group_size) THEN
         color_sub = num_pe/gw_env%tensor_group_size - 1
      END IF
      gw_env%tensor_group_color = color_sub

      ALLOCATE (gw_env%para_env_tensor)
      CALL gw_env%para_env_tensor%from_split(para_env, color_sub)

      num_t_groups = para_env%num_pe/gw_env%tensor_group_size
      gw_env%num_tensor_groups = num_t_groups

      ! create a square mesh of tensor groups for iatom and jatom; code from blacs_env_create
      gcd_max = -1
      DO ipe = 1, CEILING(SQRT(REAL(num_t_groups, dp)))
         jpe = num_t_groups/ipe
         IF (ipe*jpe .NE. num_t_groups) CYCLE
         IF (gcd(ipe, jpe) >= gcd_max) THEN
            nprow = ipe
            npcol = jpe
            gcd_max = gcd(ipe, jpe)
         END IF
      END DO

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

      IF (MODULO(gw_env%n_atom, nprow) == 0) THEN
         i_atoms_per_group = gw_env%n_atom/nprow
      ELSE
         i_atoms_per_group = gw_env%n_atom/nprow + 1
      END IF

      IF (MODULO(gw_env%n_atom, npcol) == 0) THEN
         j_atoms_per_group = gw_env%n_atom/npcol
      ELSE
         j_atoms_per_group = gw_env%n_atom/npcol + 1
      END IF

      gw_env%atom_i_start = iprow*i_atoms_per_group + 1
      gw_env%atom_i_end = MIN((iprow + 1)*i_atoms_per_group, gw_env%n_atom)

      gw_env%atom_j_start = ipcol*j_atoms_per_group + 1
      gw_env%atom_j_end = MIN((ipcol + 1)*j_atoms_per_group, gw_env%n_atom)

      ! JW to think: Maybe, we need to change the number of atoms per i and per j cycle
      !              to optimize the performance
      gw_env%n_atom_ij = 1

      color_sub = para_env%mepos/gw_env%kp_diag_group_size
      IF (MODULO(num_pe, gw_env%kp_diag_group_size) .NE. 0 .AND. &
          color_sub == num_pe/gw_env%tensor_group_size) THEN
         color_sub = num_pe/gw_env%tensor_group_size - 1
      END IF
      gw_env%kp_diag_group_color = color_sub

      ALLOCATE (gw_env%para_env_kp_diag)
      CALL gw_env%para_env_kp_diag%from_split(para_env, color_sub)

      CALL timestop(handle)

   END SUBROUTINE set_parallelization_parameters

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE set_heuristic_parameters(gw_env, qs_env)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_heuristic_parameters'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! for generating numerically stable minimax Fourier integration weights
      gw_env%num_points_per_magnitude = 200
      gw_env%regularization_minimax = 1.0E-6_dp
      gw_env%stabilize_exp = 70.0_dp
      gw_env%eps_atom_grid_2d_mat = 1.0E-50_dp

      ! only use interval ω∈[0,20eV] for virt, or ω∈[-20eV,0] for occ for use in
      ! analytic continuation of self-energy Σ^c_n(iω,k) -> Σ^c_n(ϵ,k)
      gw_env%freq_max_fit = 20.0_dp/evolt

      ! use a 16-parameter Padé fit
      gw_env%nparam_pade = 16

      ! minimum block size for tensor operations, taken from MP2/RPA input
      gw_env%min_block_size = 5

      ! resolution of the identity with the truncated Coulomb metric, cutoff radius 3 Angström
      gw_env%ri_metric%potential_type = do_potential_truncated
      gw_env%ri_metric%omega = 0.0_dp
      ! cutoff radius = 3 Angström
      gw_env%ri_metric%cutoff_radius = 3.0_dp/angstrom
      gw_env%ri_metric%filename = "t_c_g.dat"

      gw_env%eps_eigval_mat_RI = 1.0E-3_dp
      gw_env%regularization_RI = 1.0E-2_dp

      ! truncated Coulomb operator for exchange self-energy
      ! (see details in Guidon, VandeVondele, Hutter, JCTC 5, 3010 (2009) and references therein)
      CALL setup_trunc_coulomb_pot_for_exchange_self_energy(qs_env, gw_env%trunc_coulomb, &
                                                            rel_cutoff_trunc_coulomb_ri_x=0.5_dp)

      CALL timestop(handle)

   END SUBROUTINE set_heuristic_parameters

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE print_header_and_input_parameters(gw_env)

      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_header_and_input_parameters'

      INTEGER                                            :: handle, u

      CALL timeset(routineN, handle)

      gw_env%unit_nr = cp_logger_get_default_io_unit()

      u = gw_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, *) ' '
         WRITE (u, '(T2,A)') '-------------------------------------------------------------------------------'
         WRITE (u, '(T2,A)') '-                                                                             -'
         WRITE (u, '(T2,A)') '-                              GW CALCULATION                                 -'
         WRITE (u, '(T2,A)') '-                                                                             -'
         WRITE (u, '(T2,A)') '-------------------------------------------------------------------------------'
         WRITE (u, '(T2,A)') ' '
         WRITE (u, '(T2,A,I45)') 'Input: Number of time/freq. points', gw_env%num_time_freq_points
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_header_and_input_parameters

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE diagonalize_ks_matrix(gw_env)
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_ks_matrix'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! use work matrices because the matrices are overwritten in cp_fm_geeig_canon
      CALL cp_fm_to_fm(gw_env%fm_ks_Gamma(1), gw_env%fm_work_mo(1))
      CALL cp_fm_to_fm(gw_env%fm_s_Gamma, gw_env%fm_work_mo(2))

      ! diagonalize the Kohn-Sham matrix to obtain MO coefficients and SCF eigenvalues
      ! (at the Gamma-point)
      CALL cp_fm_geeig_canon(gw_env%fm_work_mo(1), &
                             gw_env%fm_work_mo(2), &
                             gw_env%fm_mo_coeff_Gamma(1), &
                             gw_env%eigenvalues_scf_Gamma, &
                             gw_env%fm_work_mo(3), &
                             gw_env%eps_eigval_mat_s)

      gw_env%e_fermi(1) = 0.5_dp*(gw_env%eigenvalues_scf_Gamma(gw_env%n_occ(1)) + &
                                  gw_env%eigenvalues_scf_Gamma(gw_env%n_occ(1) + 1))

      CALL timestop(handle)

   END SUBROUTINE diagonalize_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE setup_time_and_frequency_minimax_grid(gw_env)
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_time_and_frequency_minimax_grid'

      INTEGER                                            :: handle, homo, i_w, ierr, j_w, n_mo, &
                                                            num_time_freq_points, u
      REAL(KIND=dp)                                      :: E_max, E_min, E_range, max_error_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: points_and_weights

      CALL timeset(routineN, handle)

      homo = gw_env%n_occ(1)
      n_mo = gw_env%n_ao
      num_time_freq_points = gw_env%num_time_freq_points

      ! minimum and maximum difference between eigenvalues of unoccupied and an occupied MOs
      E_min = gw_env%eigenvalues_scf_Gamma(homo + 1) - gw_env%eigenvalues_scf_Gamma(homo)
      E_max = gw_env%eigenvalues_scf_Gamma(n_mo) - gw_env%eigenvalues_scf_Gamma(1)

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
      gw_env%imag_freq_points(:) = points_and_weights(1:num_time_freq_points)*E_min

      ! determine number of fit points in the interval [0,ω_max] for virt, or [-ω_max,0] for occ
      gw_env%num_freq_points_fit = 0
      DO i_w = 1, gw_env%num_time_freq_points
         IF (gw_env%imag_freq_points(i_w) < gw_env%freq_max_fit) THEN
            gw_env%num_freq_points_fit = gw_env%num_freq_points_fit + 1
         END IF
      END DO

      ! iω values for the analytic continuation Σ^c_n(iω,k) -> Σ^c_n(ϵ,k)
      ALLOCATE (gw_env%imag_freq_points_fit(gw_env%num_freq_points_fit))
      j_w = 0
      DO i_w = 1, gw_env%num_time_freq_points
         IF (gw_env%imag_freq_points(i_w) < gw_env%freq_max_fit) THEN
            j_w = j_w + 1
            gw_env%imag_freq_points_fit(j_w) = gw_env%imag_freq_points(i_w)
         END IF
      END DO

      ! reset the number of Padé parameters if smaller than the number of
      ! imaginary-frequency points for the fit
      IF (gw_env%num_freq_points_fit < gw_env%nparam_pade) THEN
         gw_env%nparam_pade = gw_env%num_freq_points_fit
      END IF

      ! time points
      IF (num_time_freq_points .LE. 20) THEN
         CALL get_exp_minimax_coeff(num_time_freq_points, E_range, points_and_weights)
      ELSE
         CALL get_exp_minimax_coeff_gw(num_time_freq_points, E_range, points_and_weights)
      END IF

      gw_env%imag_time_points(:) = points_and_weights(1:num_time_freq_points)/(2.0_dp*E_min)

      DEALLOCATE (points_and_weights)

      ! in minimax grids, Fourier transforms t -> w and w -> t are split using
      ! e^(iwt) = cos(wt) + i sin(wt); we thus calculate weights for trafos with a cos and
      ! sine prefactor; details in Azizi, Wilhelm, Golze, Giantomassi, Panades-Barrueta,
      ! Rinke, Draxl, Gonze et al., 2 publications

      ! cosine transform weights imaginary time to imaginary frequency
      CALL get_l_sq_wghts_cos_tf_t_to_w(num_time_freq_points, &
                                        gw_env%imag_time_points, &
                                        gw_env%weights_cos_t_to_w, &
                                        gw_env%imag_freq_points, &
                                        E_min, E_max, max_error_min, &
                                        gw_env%num_points_per_magnitude, &
                                        gw_env%regularization_minimax)

      ! cosine transform weights imaginary frequency to imaginary time
      CALL get_l_sq_wghts_cos_tf_w_to_t(num_time_freq_points, &
                                        gw_env%imag_time_points, &
                                        gw_env%weights_cos_w_to_t, &
                                        gw_env%imag_freq_points, &
                                        E_min, E_max, max_error_min, &
                                        gw_env%num_points_per_magnitude, &
                                        gw_env%regularization_minimax)

      ! sine transform weights imaginary time to imaginary frequency
      CALL get_l_sq_wghts_sin_tf_t_to_w(num_time_freq_points, &
                                        gw_env%imag_time_points, &
                                        gw_env%weights_sin_t_to_w, &
                                        gw_env%imag_freq_points, &
                                        E_min, E_max, max_error_min, &
                                        gw_env%num_points_per_magnitude, &
                                        gw_env%regularization_minimax)

      u = gw_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,F44.2)') 'SCF direct band gap at Γ-point (eV)', E_min*evolt
         WRITE (u, '(T2,A,F42.2)') 'Max. SCF eigval diff. at Γ-point (eV)', E_max*evolt
         WRITE (u, '(T2,A,F55.2)') 'E-Range for minimax grid', E_range
         WRITE (u, '(T2,A,I27)') 'Number of Padé parameters for analytic continuation:', &
            gw_env%nparam_pade
         WRITE (u, '(T2,A)') ''
      END IF
      CALL timestop(handle)

   END SUBROUTINE setup_time_and_frequency_minimax_grid

END MODULE gw_utils

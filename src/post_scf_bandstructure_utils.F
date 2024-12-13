!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 07.2023
! **************************************************************************************************
MODULE post_scf_bandstructure_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig,&
                                              cp_cfm_geeig_canon,&
                                              cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
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
   USE cp_fm_diag,                      ONLY: cp_fm_geeig_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_parser_methods,               ONLY: read_float_object
   USE input_constants,                 ONLY: int_ldos_z,&
                                              large_cell_Gamma,&
                                              small_cell_full_kp
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              max_line_length
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE post_scf_bandstructure_types,    ONLY: band_edges_type,&
                                              post_scf_bandstructure_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE rpa_gw_im_time_util,             ONLY: compute_weight_re_im,&
                                              get_atom_index_from_basis_function_index
   USE scf_control_types,               ONLY: scf_control_type
   USE soc_pseudopotential_methods,     ONLY: V_SOC_xyz_from_pseudopotential,&
                                              remove_soc_outside_energy_window_mo
   USE soc_pseudopotential_utils,       ONLY: add_cfm_submat,&
                                              add_dbcsr_submat,&
                                              cfm_add_on_diag,&
                                              create_cfm_double,&
                                              get_cfm_submat
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: create_and_init_bs_env, &
             dos_pdos_ldos, cfm_ikp_from_fm_Gamma, MIC_contribution_from_ikp, &
             compute_xkp, kpoint_init_cell_index_simple, rsmat_to_kp, soc, &
             get_VBM_CBM_bandgaps, get_all_VBM_CBM_bandgaps

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'post_scf_bandstructure_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param post_scf_bandstructure_section ...
! **************************************************************************************************
   SUBROUTINE create_and_init_bs_env(qs_env, bs_env, post_scf_bandstructure_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_and_init_bs_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (bs_env)

      CALL print_header(bs_env)

      CALL read_bandstructure_input_parameters(bs_env, post_scf_bandstructure_section)

      CALL get_parameters_from_qs_env(qs_env, bs_env)

      CALL set_heuristic_parameters(bs_env)

      SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
      CASE (large_cell_Gamma)

         CALL setup_kpoints_DOS_large_cell_Gamma(qs_env, bs_env, bs_env%kpoints_DOS)

         CALL allocate_and_fill_fm_ks_fm_s(qs_env, bs_env)

         CALL diagonalize_ks_matrix(bs_env)

         CALL check_positive_definite_overlap_mat(bs_env, qs_env)

      CASE (small_cell_full_kp)

         CALL setup_kpoints_scf_desymm(qs_env, bs_env, bs_env%kpoints_scf_desymm, .TRUE.)
         CALL setup_kpoints_scf_desymm(qs_env, bs_env, bs_env%kpoints_scf_desymm_2, .FALSE.)

         CALL setup_kpoints_DOS_small_cell_full_kp(bs_env, bs_env%kpoints_DOS)

         CALL allocate_and_fill_fm_ks_fm_s(qs_env, bs_env)

         CALL compute_cfm_mo_coeff_kp_and_eigenval_scf_kp(qs_env, bs_env)

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE create_and_init_bs_env

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param bs_sec ...
! **************************************************************************************************
   SUBROUTINE read_bandstructure_input_parameters(bs_env, bs_sec)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: bs_sec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_bandstructure_input_parameters'

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: string_ptr
      CHARACTER(LEN=max_line_length)                     :: error_msg
      INTEGER                                            :: handle, i, ikp
      TYPE(section_vals_type), POINTER                   :: gw_sec, kp_bs_sec, ldos_sec, soc_sec

      CALL timeset(routineN, handle)

      NULLIFY (gw_sec)
      gw_sec => section_vals_get_subs_vals(bs_sec, "GW")
      CALL section_vals_get(gw_sec, explicit=bs_env%do_gw)

      NULLIFY (soc_sec)
      soc_sec => section_vals_get_subs_vals(bs_sec, "SOC")
      CALL section_vals_get(soc_sec, explicit=bs_env%do_soc)

      CALL section_vals_val_get(soc_sec, "ENERGY_WINDOW", r_val=bs_env%energy_window_soc)

      CALL section_vals_val_get(bs_sec, "DOS%KPOINTS", i_vals=bs_env%nkp_grid_DOS_input)
      CALL section_vals_val_get(bs_sec, "DOS%ENERGY_WINDOW", r_val=bs_env%energy_window_DOS)
      CALL section_vals_val_get(bs_sec, "DOS%ENERGY_STEP", r_val=bs_env%energy_step_DOS)
      CALL section_vals_val_get(bs_sec, "DOS%BROADENING", r_val=bs_env%broadening_DOS)

      NULLIFY (ldos_sec)
      ldos_sec => section_vals_get_subs_vals(bs_sec, "DOS%LDOS")
      CALL section_vals_get(ldos_sec, explicit=bs_env%do_ldos)

      CALL section_vals_val_get(ldos_sec, "INTEGRATION", i_val=bs_env%int_ldos_xyz)
      CALL section_vals_val_get(ldos_sec, "BIN_MESH", i_vals=bs_env%bin_mesh)

      NULLIFY (kp_bs_sec)
      kp_bs_sec => section_vals_get_subs_vals(bs_sec, "BANDSTRUCTURE_PATH")
      CALL section_vals_val_get(kp_bs_sec, "NPOINTS", i_val=bs_env%input_kp_bs_npoints)
      CALL section_vals_val_get(kp_bs_sec, "SPECIAL_POINT", n_rep_val=bs_env%input_kp_bs_n_sp_pts)

      ! read special points for band structure
      ALLOCATE (bs_env%xkp_special(3, bs_env%input_kp_bs_n_sp_pts))
      DO ikp = 1, bs_env%input_kp_bs_n_sp_pts
         CALL section_vals_val_get(kp_bs_sec, "SPECIAL_POINT", i_rep_val=ikp, c_vals=string_ptr)
         CPASSERT(SIZE(string_ptr(:), 1) == 4)
         DO i = 1, 3
            CALL read_float_object(string_ptr(i + 1), bs_env%xkp_special(i, ikp), error_msg)
            IF (LEN_TRIM(error_msg) > 0) CPABORT(TRIM(error_msg))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE read_bandstructure_input_parameters

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_header(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'print_header'

      INTEGER                                            :: handle, u

      CALL timeset(routineN, handle)

      bs_env%unit_nr = cp_logger_get_default_io_unit()

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, *) ' '
         WRITE (u, '(T2,2A)') '-------------------------------------------------', &
            '------------------------------'
         WRITE (u, '(T2,2A)') '-                                                ', &
            '                             -'
         WRITE (u, '(T2,2A)') '-                          BANDSTRUCTURE CALCULATION', &
            '                          -'
         WRITE (u, '(T2,2A)') '-                                                ', &
            '                             -'
         WRITE (u, '(T2,2A)') '--------------------------------------------------', &
            '-----------------------------'
         WRITE (u, '(T2,A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_header

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_DOS_large_cell_Gamma(qs_env, bs_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_DOS_large_cell_Gamma'

      INTEGER                                            :: handle, i_dim, i_kp_in_line, &
                                                            i_special_kp, ikk, n_kp_in_line, &
                                                            n_special_kp, nkp, nkp_only_bs, &
                                                            nkp_only_DOS, u
      INTEGER, DIMENSION(3)                              :: nkp_grid, periodic

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      kpoints%kp_scheme = "GENERAL"

      n_special_kp = bs_env%input_kp_bs_n_sp_pts
      n_kp_in_line = bs_env%input_kp_bs_npoints

      periodic(1:3) = bs_env%periodic(1:3)

      DO i_dim = 1, 3

         CPASSERT(periodic(i_dim) == 0 .OR. periodic(i_dim) == 1)

         IF (bs_env%nkp_grid_DOS_input(i_dim) < 0) THEN
            IF (periodic(i_dim) == 1) nkp_grid(i_dim) = 2
            IF (periodic(i_dim) == 0) nkp_grid(i_dim) = 1
         ELSE
            nkp_grid(i_dim) = bs_env%nkp_grid_DOS_input(i_dim)
         END IF

      END DO

      ! use the k <-> -k symmetry to reduce the number of kpoints
      IF (nkp_grid(1) > 1) THEN
         nkp_only_DOS = (nkp_grid(1) + 1)/2*nkp_grid(2)*nkp_grid(3)
      ELSE IF (nkp_grid(2) > 1) THEN
         nkp_only_DOS = nkp_grid(1)*(nkp_grid(2) + 1)/2*nkp_grid(3)
      ELSE IF (nkp_grid(3) > 1) THEN
         nkp_only_DOS = nkp_grid(1)*nkp_grid(2)*(nkp_grid(3) + 1)/2
      ELSE
         nkp_only_DOS = 1
      END IF

      ! we will compute the GW QP levels for all k's in the bandstructure path but also
      ! for all k-points from the SCF (e.g. for DOS or for self-consistent GW)
      IF (n_special_kp > 0) THEN
         nkp_only_bs = n_kp_in_line*(n_special_kp - 1) + 1
      ELSE
         nkp_only_bs = 0
      END IF

      nkp = nkp_only_DOS + nkp_only_bs

      kpoints%nkp_grid(1:3) = nkp_grid(1:3)
      kpoints%nkp = nkp

      bs_env%nkp_bs_and_DOS = nkp
      bs_env%nkp_only_bs = nkp_only_bs
      bs_env%nkp_only_DOS = nkp_only_DOS

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))
      kpoints%wkp(1:nkp_only_DOS) = 1.0_dp/REAL(nkp_only_DOS, KIND=dp)

      CALL compute_xkp(kpoints%xkp, 1, nkp_only_DOS, nkp_grid)

      IF (n_special_kp > 0) THEN
         kpoints%xkp(1:3, nkp_only_DOS + 1) = bs_env%xkp_special(1:3, 1)
         ikk = nkp_only_DOS + 1
         DO i_special_kp = 2, n_special_kp
            DO i_kp_in_line = 1, n_kp_in_line
               ikk = ikk + 1
               kpoints%xkp(1:3, ikk) = bs_env%xkp_special(1:3, i_special_kp - 1) + &
                                       REAL(i_kp_in_line, KIND=dp)/REAL(n_kp_in_line, KIND=dp)* &
                                       (bs_env%xkp_special(1:3, i_special_kp) - &
                                        bs_env%xkp_special(1:3, i_special_kp - 1))
               kpoints%wkp(ikk) = 0.0_dp
            END DO
         END DO
      END IF

      CALL kpoint_init_cell_index_simple(kpoints, qs_env)

      u = bs_env%unit_nr

      IF (u > 0) THEN
         IF (nkp_only_bs > 0) THEN
            WRITE (u, FMT="(T2,1A,T77,I4)") &
               "Number of special k-points for the bandstructure", n_special_kp
            WRITE (u, FMT="(T2,1A,T77,I4)") "Number of k-points for the bandstructure", nkp
            WRITE (u, FMT="(T2,1A,T69,3I4)") &
               "K-point mesh for the density of states (DOS)", nkp_grid(1:3)
         ELSE
            WRITE (u, FMT="(T2,1A,T69,3I4)") &
               "K-point mesh for the density of states (DOS) and the self-energy", nkp_grid(1:3)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_DOS_large_cell_Gamma

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

      ! we need in periodic directions at least 4 k-points in the SCF
      DO i_dim = 1, 3
         IF (bs_env%periodic(i_dim) == 1) THEN
            CPASSERT(nkp_grid(i_dim) .GE. 4)
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
      bs_env%cell_grid_scf_desymm(1:3) = cell_grid(1:3)

      IF (ASSOCIATED(kpoints%index_to_cell)) DEALLOCATE (kpoints%index_to_cell)
      IF (ASSOCIATED(kpoints%cell_to_index)) DEALLOCATE (kpoints%cell_to_index)

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
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_DOS_small_cell_full_kp(bs_env, kpoints)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_DOS_small_cell_full_kp'

      INTEGER                                            :: handle, i_kp_in_line, i_special_kp, ikk, &
                                                            n_kp_in_line, n_special_kp, nkp, &
                                                            nkp_only_bs, nkp_scf_desymm, u

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      n_special_kp = bs_env%input_kp_bs_n_sp_pts
      n_kp_in_line = bs_env%input_kp_bs_npoints
      nkp_scf_desymm = bs_env%nkp_scf_desymm

      ! we will compute the GW QP levels for all k's in the bandstructure path but also
      ! for all k-points from the SCF (e.g. for DOS or for self-consistent GW)
      IF (n_special_kp > 0) THEN
         nkp_only_bs = n_kp_in_line*(n_special_kp - 1) + 1
      ELSE
         nkp_only_bs = 0
      END IF
      nkp = nkp_only_bs + nkp_scf_desymm

      ALLOCATE (kpoints%xkp(3, nkp))
      ALLOCATE (kpoints%wkp(nkp))

      kpoints%nkp = nkp

      bs_env%nkp_bs_and_DOS = nkp
      bs_env%nkp_only_bs = nkp_only_bs
      bs_env%nkp_only_DOS = nkp_scf_desymm

      kpoints%xkp(1:3, 1:nkp_scf_desymm) = bs_env%kpoints_scf_desymm%xkp(1:3, 1:nkp_scf_desymm)
      kpoints%wkp(1:nkp_scf_desymm) = 1.0_dp/REAL(nkp_scf_desymm, KIND=dp)

      IF (n_special_kp > 0) THEN
         kpoints%xkp(1:3, nkp_scf_desymm + 1) = bs_env%xkp_special(1:3, 1)
         ikk = nkp_scf_desymm + 1
         DO i_special_kp = 2, n_special_kp
            DO i_kp_in_line = 1, n_kp_in_line
               ikk = ikk + 1
               kpoints%xkp(1:3, ikk) = bs_env%xkp_special(1:3, i_special_kp - 1) + &
                                       REAL(i_kp_in_line, KIND=dp)/REAL(n_kp_in_line, KIND=dp)* &
                                       (bs_env%xkp_special(1:3, i_special_kp) - &
                                        bs_env%xkp_special(1:3, i_special_kp - 1))
               kpoints%wkp(ikk) = 0.0_dp
            END DO
         END DO
      END IF

      IF (ASSOCIATED(kpoints%index_to_cell)) DEALLOCATE (kpoints%index_to_cell)

      ALLOCATE (kpoints%index_to_cell(bs_env%nimages_scf_desymm, 3))
      kpoints%index_to_cell(:, :) = bs_env%kpoints_scf_desymm%index_to_cell(:, :)

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (u, FMT="(T2,1A,T77,I4)") "Number of special k-points for the bandstructure", &
            n_special_kp
         WRITE (u, FMT="(T2,1A,T77,I4)") "Number of k-points for the bandstructure", nkp
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_DOS_small_cell_full_kp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_cfm_mo_coeff_kp_and_eigenval_scf_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cfm_mo_coeff_kp_and_eigenval_scf_kp'

      INTEGER                                            :: handle, ikp, ispin, nkp_bs_and_DOS
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
      ALLOCATE (bs_env%cfm_mo_coeff_kp(nkp_bs_and_DOS, bs_env%n_spin))
      ALLOCATE (bs_env%cfm_ks_kp(nkp_bs_and_DOS, bs_env%n_spin))
      ALLOCATE (bs_env%cfm_s_kp(nkp_bs_and_DOS))
      DO ikp = 1, nkp_bs_and_DOS
      DO ispin = 1, bs_env%n_spin
         CALL cp_cfm_create(bs_env%cfm_mo_coeff_kp(ikp, ispin), bs_env%cfm_work_mo%matrix_struct)
         CALL cp_cfm_create(bs_env%cfm_ks_kp(ikp, ispin), bs_env%cfm_work_mo%matrix_struct)
      END DO
      CALL cp_cfm_create(bs_env%cfm_s_kp(ikp), bs_env%cfm_work_mo%matrix_struct)
      END DO

      DO ispin = 1, bs_env%n_spin
         DO ikp = 1, nkp_bs_and_DOS

            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, ikp)

            ! h^KS^R -> h^KS(k)
            CALL rsmat_to_kp(matrix_ks, ispin, xkp, cell_to_index_scf, sab_nl, bs_env, cfm_ks)

            ! S^R -> S(k)
            CALL rsmat_to_kp(matrix_s, 1, xkp, cell_to_index_scf, sab_nl, bs_env, cfm_s)

            ! we store the complex KS matrix as fm matrix because the infrastructure for fm is
            ! much nicer compared to cfm
            CALL cp_cfm_to_cfm(cfm_ks, bs_env%cfm_ks_kp(ikp, ispin))
            CALL cp_cfm_to_cfm(cfm_s, bs_env%cfm_s_kp(ikp))

            ! Diagonalize KS-matrix via Rothaan-Hall equation:
            ! H^KS(k) C(k) = S(k) C(k) ε(k)
            CALL cp_cfm_geeig_canon(cfm_ks, cfm_s, cfm_mos, &
                                    bs_env%eigenval_scf(:, ikp, ispin), &
                                    bs_env%cfm_work_mo, bs_env%eps_eigval_mat_s)

            ! we store the complex MO coeff as fm matrix because the infrastructure for fm is
            ! much nicer compared to cfm
            CALL cp_cfm_to_cfm(cfm_mos, bs_env%cfm_mo_coeff_kp(ikp, ispin))

         END DO

         VBM = MAXVAL(bs_env%eigenval_scf(bs_env%n_occ(ispin), :, ispin))
         CBM = MINVAL(bs_env%eigenval_scf(bs_env%n_occ(ispin) + 1, :, ispin))

         bs_env%e_fermi(ispin) = 0.5_dp*(VBM + CBM)

      END DO

      CALL get_VBM_CBM_bandgaps(bs_env%band_edges_scf, bs_env%eigenval_scf, bs_env)

      CALL cp_cfm_release(cfm_ks)
      CALL cp_cfm_release(cfm_s)
      CALL cp_cfm_release(cfm_mos)

      CALL timestop(handle)

   END SUBROUTINE compute_cfm_mo_coeff_kp_and_eigenval_scf_kp

! **************************************************************************************************
!> \brief ...
!> \param mat_rs ...
!> \param ispin ...
!> \param xkp ...
!> \param cell_to_index_scf ...
!> \param sab_nl ...
!> \param bs_env ...
!> \param cfm_kp ...
!> \param imag_rs_mat ...
! **************************************************************************************************
   SUBROUTINE rsmat_to_kp(mat_rs, ispin, xkp, cell_to_index_scf, sab_nl, bs_env, cfm_kp, imag_rs_mat)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_rs
      INTEGER                                            :: ispin
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_scf
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_kp
      LOGICAL, OPTIONAL                                  :: imag_rs_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rsmat_to_kp'

      INTEGER                                            :: handle
      LOGICAL                                            :: imag_rs_mat_private
      TYPE(dbcsr_type), POINTER                          :: cmat, nsmat, rmat

      CALL timeset(routineN, handle)

      ALLOCATE (rmat, cmat, nsmat)

      imag_rs_mat_private = .FALSE.
      IF (PRESENT(imag_rs_mat)) imag_rs_mat_private = imag_rs_mat

      IF (imag_rs_mat_private) THEN
         CALL dbcsr_create(rmat, template=mat_rs(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_create(cmat, template=mat_rs(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      ELSE
         CALL dbcsr_create(rmat, template=mat_rs(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(cmat, template=mat_rs(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      END IF
      CALL dbcsr_create(nsmat, template=mat_rs(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(rmat, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmat, sab_nl)

      CALL dbcsr_set(rmat, 0.0_dp)
      CALL dbcsr_set(cmat, 0.0_dp)
      CALL rskp_transform(rmatrix=rmat, cmatrix=cmat, rsmat=mat_rs, ispin=ispin, &
                          xkp=xkp, cell_to_index=cell_to_index_scf, sab_nl=sab_nl)

      CALL dbcsr_desymmetrize(rmat, nsmat)
      CALL copy_dbcsr_to_fm(nsmat, bs_env%fm_work_mo(1))
      CALL dbcsr_desymmetrize(cmat, nsmat)
      CALL copy_dbcsr_to_fm(nsmat, bs_env%fm_work_mo(2))
      CALL cp_fm_to_cfm(bs_env%fm_work_mo(1), bs_env%fm_work_mo(2), cfm_kp)

      CALL dbcsr_deallocate_matrix(rmat)
      CALL dbcsr_deallocate_matrix(cmat)
      CALL dbcsr_deallocate_matrix(nsmat)

      CALL timestop(handle)

   END SUBROUTINE rsmat_to_kp

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE diagonalize_ks_matrix(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_ks_matrix'

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: CBM, VBM

      CALL timeset(routineN, handle)

      ALLOCATE (bs_env%eigenval_scf_Gamma(bs_env%n_ao, bs_env%n_spin))

      DO ispin = 1, bs_env%n_spin

         ! use work matrices because the matrices are overwritten in cp_fm_geeig_canon
         CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(ispin), bs_env%fm_work_mo(1))
         CALL cp_fm_to_fm(bs_env%fm_s_Gamma, bs_env%fm_work_mo(2))

         ! diagonalize the Kohn-Sham matrix to obtain MO coefficients and SCF eigenvalues
         ! (at the Gamma-point)
         CALL cp_fm_geeig_canon(bs_env%fm_work_mo(1), &
                                bs_env%fm_work_mo(2), &
                                bs_env%fm_mo_coeff_Gamma(ispin), &
                                bs_env%eigenval_scf_Gamma(:, ispin), &
                                bs_env%fm_work_mo(3), &
                                bs_env%eps_eigval_mat_s)

         VBM = bs_env%eigenval_scf_Gamma(bs_env%n_occ(ispin), ispin)
         CBM = bs_env%eigenval_scf_Gamma(bs_env%n_occ(ispin) + 1, ispin)

         bs_env%band_edges_scf_Gamma(ispin)%VBM = VBM
         bs_env%band_edges_scf_Gamma(ispin)%CBM = CBM
         bs_env%e_fermi(ispin) = 0.5_dp*(VBM + CBM)

      END DO

      CALL timestop(handle)

   END SUBROUTINE diagonalize_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE check_positive_definite_overlap_mat(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_positive_definite_overlap_mat'

      INTEGER                                            :: handle, ikp, info, u
      TYPE(cp_cfm_type)                                  :: cfm_s_ikp

      CALL timeset(routineN, handle)

      DO ikp = 1, bs_env%kpoints_DOS%nkp

         ! get S_µν(k_i) from S_µν(k=0)
         CALL cfm_ikp_from_fm_Gamma(cfm_s_ikp, bs_env%fm_s_Gamma, &
                                    ikp, qs_env, bs_env%kpoints_DOS, "ORB")

         ! check whether S_µν(k_i) is positive definite
         CALL cp_cfm_cholesky_decompose(matrix=cfm_s_ikp, n=bs_env%n_ao, info_out=info)

         ! check if Cholesky decomposition failed (Cholesky decomposition only works for
         ! positive definite matrices
         IF (info .NE. 0) THEN
            u = bs_env%unit_nr

            IF (u > 0) THEN
               WRITE (u, FMT="(T2,A)") ""
               WRITE (u, FMT="(T2,A)") "ERROR: The Cholesky decomposition "// &
                  "of the k-point overlap matrix failed. This is"
               WRITE (u, FMT="(T2,A)") "because the algorithm is "// &
                  "only correct in the limit of large cells. The cell of "
               WRITE (u, FMT="(T2,A)") "the calculation is too small. "// &
                  "Use MULTIPLE_UNIT_CELL to create a larger cell "
               WRITE (u, FMT="(T2,A)") "and to prevent this error."
            END IF

            CALL bs_env%para_env%sync()
            CPABORT("Please see information on the error above.")

         END IF ! Cholesky decomposition failed

      END DO ! ikp

      CALL cp_cfm_release(cfm_s_ikp)

      CALL timestop(handle)

   END SUBROUTINE check_positive_definite_overlap_mat

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_parameters_from_qs_env(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_parameters_from_qs_env'

      INTEGER                                            :: color_sub, handle, homo, n_ao, n_atom, u
      INTEGER, DIMENSION(3)                              :: periodic
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      scf_control=scf_control, &
                      mos=mos)

      bs_env%n_spin = dft_control%nspins
      IF (bs_env%n_spin == 1) bs_env%spin_degeneracy = 2.0_dp
      IF (bs_env%n_spin == 2) bs_env%spin_degeneracy = 1.0_dp

      CALL get_mo_set(mo_set=mos(1), nao=n_ao, homo=homo)
      bs_env%n_ao = n_ao
      bs_env%n_occ(1:2) = homo
      bs_env%n_vir(1:2) = n_ao - homo

      IF (bs_env%n_spin == 2) THEN
         CALL get_mo_set(mo_set=mos(2), homo=homo)
         bs_env%n_occ(2) = homo
         bs_env%n_vir(2) = n_ao - homo
      END IF

      bs_env%eps_eigval_mat_s = scf_control%eps_eigval

      ! get para_env from qs_env (bs_env%para_env is identical to para_env in qs_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      color_sub = 0
      ALLOCATE (bs_env%para_env)
      CALL bs_env%para_env%from_split(para_env, color_sub)

      CALL get_qs_env(qs_env, particle_set=particle_set)

      n_atom = SIZE(particle_set)
      bs_env%n_atom = n_atom

      CALL get_qs_env(qs_env=qs_env, cell=cell)
      CALL get_cell(cell=cell, periodic=periodic, h=hmat)
      bs_env%periodic(1:3) = periodic(1:3)
      bs_env%hmat(1:3, 1:3) = hmat
      bs_env%nimages_scf = dft_control%nimages
      IF (dft_control%nimages == 1) THEN
         bs_env%small_cell_full_kp_or_large_cell_Gamma = large_cell_Gamma
      ELSE IF (dft_control%nimages > 1) THEN
         bs_env%small_cell_full_kp_or_large_cell_Gamma = small_cell_full_kp
      ELSE
         CPABORT("Wrong number of cells from DFT calculation.")
      END IF

      u = bs_env%unit_nr

      ! Marek : Get and save the rtp method
      CALL get_qs_env(qs_env=qs_env, input=input)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%RTBSE%_SECTION_PARAMETERS_", i_val=bs_env%rtp_method)

      IF (u > 0) THEN
         WRITE (u, FMT="(T2,2A,T73,I8)") "Number of occupied molecular orbitals (MOs) ", &
            "= Number of occupied bands", homo
         WRITE (u, FMT="(T2,2A,T73,I8)") "Number of unoccupied (= virtual) MOs ", &
            "= Number of unoccupied bands", n_ao - homo
         WRITE (u, FMT="(T2,A,T73,I8)") "Number of Gaussian basis functions for MOs", n_ao
         IF (bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp) THEN
            WRITE (u, FMT="(T2,2A,T73,I8)") "Number of cells considered in the DFT ", &
               "calculation", bs_env%nimages_scf
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_parameters_from_qs_env

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE set_heuristic_parameters(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_heuristic_parameters'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      bs_env%n_bins_max_for_printing = 5000

      CALL timestop(handle)

   END SUBROUTINE set_heuristic_parameters

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_fm_ks_fm_s(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_fill_fm_ks_fm_s'

      INTEGER                                            :: handle, i_work, ispin
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      matrix_ks_kp=matrix_ks, &
                      matrix_s_kp=matrix_s)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=bs_env%n_ao, &
                               ncol_global=bs_env%n_ao, para_env=para_env)

      DO i_work = 1, SIZE(bs_env%fm_work_mo)
         CALL cp_fm_create(bs_env%fm_work_mo(i_work), fm_struct)
      END DO

      CALL cp_cfm_create(bs_env%cfm_work_mo, fm_struct)
      CALL cp_cfm_create(bs_env%cfm_work_mo_2, fm_struct)

      CALL cp_fm_create(bs_env%fm_s_Gamma, fm_struct)
      CALL copy_dbcsr_to_fm(matrix_s(1, 1)%matrix, bs_env%fm_s_Gamma)

      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(bs_env%fm_ks_Gamma(ispin), fm_struct)
         CALL copy_dbcsr_to_fm(matrix_ks(ispin, 1)%matrix, bs_env%fm_ks_Gamma(ispin))
         CALL cp_fm_create(bs_env%fm_mo_coeff_Gamma(ispin), fm_struct)
      END DO

      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (bs_env%mat_ao_ao%matrix)
      ALLOCATE (bs_env%mat_ao_ao%matrix)
      CALL dbcsr_create(bs_env%mat_ao_ao%matrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ALLOCATE (bs_env%eigenval_scf(bs_env%n_ao, bs_env%nkp_bs_and_DOS, bs_env%n_spin))

      CALL timestop(handle)

   END SUBROUTINE allocate_and_fill_fm_ks_fm_s

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE dos_pdos_ldos(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dos_pdos_ldos'

      INTEGER                                            :: handle, homo, homo_1, homo_2, &
                                                            homo_spinor, ikp, ikp_for_file, ispin, &
                                                            n_ao, n_E, nkind, nkp
      LOGICAL                                            :: is_bandstruc_kpoint, print_DOS_kpoints, &
                                                            print_ikp
      REAL(KIND=dp)                                      :: broadening, E_max, E_max_G0W0, E_min, &
                                                            E_min_G0W0, E_total_window, &
                                                            energy_step_DOS, energy_window_DOS, t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: DOS_G0W0, DOS_G0W0_SOC, DOS_scf, DOS_scf_SOC, &
         eigenval, eigenval_spinor, eigenval_spinor_G0W0, eigenval_spinor_no_SOC
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: PDOS_G0W0, PDOS_G0W0_SOC, PDOS_scf, &
                                                            PDOS_scf_SOC, proj_mo_on_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_G0W0_2d, LDOS_scf_2d, &
                                                            LDOS_scf_2d_SOC
      TYPE(band_edges_type)                              :: band_edges_G0W0, band_edges_G0W0_SOC, &
                                                            band_edges_scf, band_edges_scf_guess, &
                                                            band_edges_scf_SOC
      TYPE(cp_cfm_type) :: cfm_ks_ikp, cfm_ks_ikp_spinor, cfm_mos_ikp_spinor, cfm_s_ikp, &
         cfm_s_ikp_copy, cfm_s_ikp_spinor, cfm_s_ikp_spinor_copy, cfm_SOC_ikp_spinor, &
         cfm_spinor_wf_ikp, cfm_work_ikp, cfm_work_ikp_spinor
      TYPE(cp_cfm_type), DIMENSION(2)                    :: cfm_mos_ikp

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao

      energy_window_DOS = bs_env%energy_window_DOS
      energy_step_DOS = bs_env%energy_step_DOS
      broadening = bs_env%broadening_DOS

      ! if we have done GW or a full kpoint SCF, we already have the band edges
      IF (bs_env%do_gw .OR. &
          bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp) THEN
         band_edges_scf = bs_env%band_edges_scf
         band_edges_scf_guess = band_edges_scf
      ELSE

         IF (bs_env%n_spin == 1) THEN
            homo = bs_env%n_occ(1)
            band_edges_scf_guess%VBM = bs_env%eigenval_scf_Gamma(homo, 1)
            band_edges_scf_guess%CBM = bs_env%eigenval_scf_Gamma(homo + 1, 1)
         ELSE
            homo_1 = bs_env%n_occ(1)
            homo_2 = bs_env%n_occ(2)
            band_edges_scf_guess%VBM = MAX(bs_env%eigenval_scf_Gamma(homo_1, 1), &
                                           bs_env%eigenval_scf_Gamma(homo_2, 2))
            band_edges_scf_guess%CBM = MIN(bs_env%eigenval_scf_Gamma(homo_1 + 1, 1), &
                                           bs_env%eigenval_scf_Gamma(homo_2 + 1, 2))
         END IF

         ! initialization
         band_edges_scf%VBM = -1000.0_dp
         band_edges_scf%CBM = 1000.0_dp
         band_edges_scf%DBG = 1000.0_dp
      END IF

      E_min = band_edges_scf_guess%VBM - 0.5_dp*energy_window_DOS
      E_max = band_edges_scf_guess%CBM + 0.5_dp*energy_window_DOS

      IF (bs_env%do_gw) THEN
         band_edges_G0W0 = bs_env%band_edges_G0W0
         E_min_G0W0 = band_edges_G0W0%VBM - 0.5_dp*energy_window_DOS
         E_max_G0W0 = band_edges_G0W0%CBM + 0.5_dp*energy_window_DOS
         E_min = MIN(E_min, E_min_G0W0)
         E_max = MAX(E_max, E_max_G0W0)
      END IF

      E_total_window = E_max - E_min

      n_E = INT(E_total_window/energy_step_DOS)

      ALLOCATE (DOS_scf(n_E))
      DOS_scf(:) = 0.0_dp
      ALLOCATE (DOS_scf_SOC(n_E))
      DOS_scf_SOC(:) = 0.0_dp

      CALL get_qs_env(qs_env, nkind=nkind)

      ALLOCATE (PDOS_scf(n_E, nkind))
      PDOS_scf(:, :) = 0.0_dp
      ALLOCATE (PDOS_scf_SOC(n_E, nkind))
      PDOS_scf_SOC(:, :) = 0.0_dp

      ALLOCATE (proj_mo_on_kind(n_ao, nkind))
      proj_mo_on_kind(:, :) = 0.0_dp

      ALLOCATE (eigenval(n_ao))
      ALLOCATE (eigenval_spinor(2*n_ao))
      ALLOCATE (eigenval_spinor_no_SOC(2*n_ao))
      ALLOCATE (eigenval_spinor_G0W0(2*n_ao))

      IF (bs_env%do_gw) THEN

         ALLOCATE (DOS_G0W0(n_E))
         DOS_G0W0(:) = 0.0_dp
         ALLOCATE (DOS_G0W0_SOC(n_E))
         DOS_G0W0_SOC(:) = 0.0_dp

         ALLOCATE (PDOS_G0W0(n_E, nkind))
         PDOS_G0W0(:, :) = 0.0_dp
         ALLOCATE (PDOS_G0W0_SOC(n_E, nkind))
         PDOS_G0W0_SOC(:, :) = 0.0_dp

      END IF

      CALL cp_cfm_create(cfm_mos_ikp(1), bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_cfm_create(cfm_mos_ikp(2), bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_cfm_create(cfm_work_ikp, bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_cfm_create(cfm_s_ikp_copy, bs_env%fm_ks_Gamma(1)%matrix_struct)

      IF (bs_env%do_soc) THEN

         CALL cp_cfm_create(cfm_mos_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
         CALL cp_cfm_create(cfm_work_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
         CALL cp_cfm_create(cfm_s_ikp_spinor_copy, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
         CALL cp_cfm_create(cfm_ks_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
         CALL cp_cfm_create(cfm_SOC_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
         CALL cp_cfm_create(cfm_s_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
         CALL cp_cfm_create(cfm_spinor_wf_ikp, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)

         homo_spinor = bs_env%n_occ(1) + bs_env%n_occ(bs_env%n_spin)

         band_edges_scf_SOC%VBM = -1000.0_dp
         band_edges_scf_SOC%CBM = 1000.0_dp
         band_edges_scf_SOC%DBG = 1000.0_dp

         IF (bs_env%do_gw) THEN
            band_edges_G0W0_SOC%VBM = -1000.0_dp
            band_edges_G0W0_SOC%CBM = 1000.0_dp
            band_edges_G0W0_SOC%DBG = 1000.0_dp
         END IF

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(A)') ''
            WRITE (bs_env%unit_nr, '(T2,A,F43.1,A)') 'SOC requested, SOC energy window:', &
               bs_env%energy_window_soc*evolt, ' eV'
         END IF

      END IF

      IF (bs_env%do_ldos) THEN
         CPASSERT(bs_env%int_ldos_xyz == int_ldos_z)
      END IF

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ''
      END IF

      IF (bs_env%small_cell_full_kp_or_large_cell_Gamma == small_cell_full_kp) THEN
         CALL cp_cfm_create(cfm_ks_ikp, bs_env%cfm_ks_kp(1, 1)%matrix_struct)
         CALL cp_cfm_create(cfm_s_ikp, bs_env%cfm_ks_kp(1, 1)%matrix_struct)
      END IF

      DO ikp = 1, bs_env%nkp_bs_and_DOS

         t1 = m_walltime()

         DO ispin = 1, bs_env%n_spin

            SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
            CASE (large_cell_Gamma)

               ! 1. get H^KS_µν(k_i) from H^KS_µν(k=0)
               CALL cfm_ikp_from_fm_Gamma(cfm_ks_ikp, bs_env%fm_ks_Gamma(ispin), &
                                          ikp, qs_env, bs_env%kpoints_DOS, "ORB")

               ! 2. get S_µν(k_i) from S_µν(k=0)
               CALL cfm_ikp_from_fm_Gamma(cfm_s_ikp, bs_env%fm_s_Gamma, &
                                          ikp, qs_env, bs_env%kpoints_DOS, "ORB")
               CALL cp_cfm_to_cfm(cfm_s_ikp, cfm_s_ikp_copy)

               ! 3. Diagonalize (Roothaan-Hall): H_KS(k_i)*C(k_i) = S(k_i)*C(k_i)*ϵ(k_i)
               CALL cp_cfm_geeig(cfm_ks_ikp, cfm_s_ikp_copy, cfm_mos_ikp(ispin), &
                                 eigenval, cfm_work_ikp)

            CASE (small_cell_full_kp)

               ! 1. get H^KS_µν(k_i)
               CALL cp_cfm_to_cfm(bs_env%cfm_ks_kp(ikp, ispin), cfm_ks_ikp)

               ! 2. get S_µν(k_i)
               CALL cp_cfm_to_cfm(bs_env%cfm_s_kp(ikp), cfm_s_ikp)

               ! 3. get C_µn(k_i) and ϵ_n(k_i)
               CALL cp_cfm_to_cfm(bs_env%cfm_mo_coeff_kp(ikp, ispin), cfm_mos_ikp(ispin))
               eigenval(:) = bs_env%eigenval_scf(:, ikp, ispin)

            END SELECT

            ! 4. Projection p_nk^A of MO ψ_nk(r) on atom type A (inspired by Mulliken charge)
            !    p_nk^A = sum_µ^A,ν C*_µ^A,n(k) S_µ^A,ν(k) C_ν,n(k)
            CALL compute_proj_mo_on_kind(proj_mo_on_kind, qs_env, cfm_mos_ikp(ispin), cfm_s_ikp)

            ! 5. DOS and PDOS
            CALL add_to_DOS_PDOS(DOS_scf, PDOS_scf, eigenval, ikp, bs_env, n_E, E_min, &
                                 proj_mo_on_kind)
            IF (bs_env%do_gw) THEN
               CALL add_to_DOS_PDOS(DOS_G0W0, PDOS_G0W0, bs_env%eigenval_G0W0(:, ikp, ispin), &
                                    ikp, bs_env, n_E, E_min, proj_mo_on_kind)
            END IF

            IF (bs_env%do_ldos) THEN
               CALL add_to_LDOS_2d(LDOS_scf_2d, qs_env, ikp, bs_env, cfm_mos_ikp(ispin), &
                                   eigenval(:), band_edges_scf_guess)

               IF (bs_env%do_gw) THEN
                  CALL add_to_LDOS_2d(LDOS_G0W0_2d, qs_env, ikp, bs_env, cfm_mos_ikp(ispin), &
                                      bs_env%eigenval_G0W0(:, ikp, 1), band_edges_G0W0)
               END IF

            END IF

            homo = bs_env%n_occ(ispin)

            band_edges_scf%VBM = MAX(band_edges_scf%VBM, eigenval(homo))
            band_edges_scf%CBM = MIN(band_edges_scf%CBM, eigenval(homo + 1))
            band_edges_scf%DBG = MIN(band_edges_scf%DBG, eigenval(homo + 1) - eigenval(homo))

         END DO ! spin

         ! now the same with spin-orbit coupling
         IF (bs_env%do_soc) THEN

            ! only print eigenvalues of DOS k-points in case no bandstructure path has been given
            print_DOS_kpoints = (bs_env%nkp_only_bs .LE. 0)
            ! in kpoints_DOS, the last nkp_only_bs are bandstructure k-points
            is_bandstruc_kpoint = (ikp > bs_env%nkp_only_DOS)
            print_ikp = print_DOS_kpoints .OR. is_bandstruc_kpoint

            IF (print_DOS_kpoints) THEN
               nkp = bs_env%nkp_only_DOS
               ikp_for_file = ikp
            ELSE
               nkp = bs_env%nkp_only_bs
               ikp_for_file = ikp - bs_env%nkp_only_DOS
            END IF

            ! compute DFT+SOC eigenvalues; based on these, compute band edges, DOS and LDOS
            CALL SOC_ev(bs_env, qs_env, ikp, bs_env%eigenval_scf, band_edges_scf, &
                        E_min, cfm_mos_ikp, DOS_scf_SOC, PDOS_scf_SOC, &
                        band_edges_scf_SOC, eigenval_spinor, cfm_spinor_wf_ikp)

            IF (.NOT. bs_env%do_gw .AND. print_ikp) THEN
               CALL write_SOC_eigenvalues(eigenval_spinor, ikp_for_file, ikp, bs_env)
            END IF

            IF (bs_env%do_ldos) THEN
               CALL add_to_LDOS_2d(LDOS_scf_2d_SOC, qs_env, ikp, bs_env, cfm_spinor_wf_ikp, &
                                   eigenval_spinor, band_edges_scf_guess, .TRUE., cfm_work_ikp)
            END IF

            IF (bs_env%do_gw) THEN

               ! compute G0W0+SOC eigenvalues; based on these, compute band edges, DOS and LDOS
               CALL SOC_ev(bs_env, qs_env, ikp, bs_env%eigenval_G0W0, band_edges_G0W0, &
                           E_min, cfm_mos_ikp, DOS_G0W0_SOC, PDOS_G0W0_SOC, &
                           band_edges_G0W0_SOC, eigenval_spinor_G0W0, cfm_spinor_wf_ikp)

               IF (print_ikp) THEN
                  ! write SCF+SOC and G0W0+SOC eigenvalues to file
                  ! SCF_and_G0W0_band_structure_for_kpoint_<ikp>_+_SOC
                  CALL write_SOC_eigenvalues(eigenval_spinor, ikp_for_file, ikp, bs_env, &
                                             eigenval_spinor_G0W0)
               END IF

            END IF ! do_gw

         END IF ! do_soc

         IF (bs_env%unit_nr > 0 .AND. m_walltime() - t1 > 20.0_dp) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T43,I5,A,I3,A,F7.1,A)') &
               'Compute DOS, LDOS for k-point ', ikp, ' /', bs_env%nkp_bs_and_DOS, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! ikp_DOS

      band_edges_scf%IDBG = band_edges_scf%CBM - band_edges_scf%VBM
      IF (bs_env%do_soc) THEN
         band_edges_scf_SOC%IDBG = band_edges_scf_SOC%CBM - band_edges_scf_SOC%VBM
         IF (bs_env%do_gw) THEN
            band_edges_G0W0_SOC%IDBG = band_edges_G0W0_SOC%CBM - band_edges_G0W0_SOC%VBM
         END IF
      END IF

      CALL write_band_edges(band_edges_scf, "SCF", bs_env)
      CALL write_dos_pdos(DOS_scf, PDOS_scf, bs_env, qs_env, "SCF", E_min, band_edges_scf%VBM)
      IF (bs_env%do_ldos) THEN
         CALL print_LDOS_main(LDOS_scf_2d, bs_env, band_edges_scf, "SCF")
      END IF

      IF (bs_env%do_soc) THEN
         CALL write_band_edges(band_edges_scf_SOC, "SCF+SOC", bs_env)
         CALL write_dos_pdos(DOS_scf_SOC, PDOS_scf_SOC, bs_env, qs_env, "SCF_SOC", &
                             E_min, band_edges_scf_SOC%VBM)
         IF (bs_env%do_ldos) THEN
            ! argument band_edges_scf is actually correct because the non-SOC band edges
            ! have been used as reference in add_to_LDOS_2d
            CALL print_LDOS_main(LDOS_scf_2d_SOC, bs_env, band_edges_scf, &
                                 "SCF_SOC")
         END IF
      END IF

      IF (bs_env%do_gw) THEN
         CALL write_band_edges(band_edges_G0W0, "G0W0", bs_env)
         CALL write_band_edges(bs_env%band_edges_HF, "Hartree-Fock with SCF orbitals", bs_env)
         CALL write_dos_pdos(DOS_G0W0, PDOS_G0W0, bs_env, qs_env, "G0W0", E_min, &
                             band_edges_G0W0%VBM)
         IF (bs_env%do_ldos) THEN
            CALL print_LDOS_main(LDOS_G0W0_2d, bs_env, band_edges_G0W0, "G0W0")
         END IF
      END IF

      IF (bs_env%do_soc .AND. bs_env%do_gw) THEN
         CALL write_band_edges(band_edges_G0W0_SOC, "G0W0+SOC", bs_env)
         CALL write_dos_pdos(DOS_G0W0_SOC, PDOS_G0W0_SOC, bs_env, qs_env, "G0W0_SOC", E_min, &
                             band_edges_G0W0_SOC%VBM)
      END IF

      CALL cp_cfm_release(cfm_s_ikp)
      CALL cp_cfm_release(cfm_ks_ikp)
      CALL cp_cfm_release(cfm_mos_ikp(1))
      CALL cp_cfm_release(cfm_mos_ikp(2))
      CALL cp_cfm_release(cfm_work_ikp)
      CALL cp_cfm_release(cfm_s_ikp_copy)

      CALL cp_cfm_release(cfm_s_ikp_spinor)
      CALL cp_cfm_release(cfm_ks_ikp_spinor)
      CALL cp_cfm_release(cfm_SOC_ikp_spinor)
      CALL cp_cfm_release(cfm_mos_ikp_spinor)
      CALL cp_cfm_release(cfm_work_ikp_spinor)
      CALL cp_cfm_release(cfm_s_ikp_spinor_copy)
      CALL cp_cfm_release(cfm_spinor_wf_ikp)

      CALL timestop(handle)

   END SUBROUTINE dos_pdos_ldos

! **************************************************************************************************
!> \brief ...
!> \param LDOS_2d ...
!> \param bs_env ...
!> \param band_edges ...
!> \param scf_gw_soc ...
! **************************************************************************************************
   SUBROUTINE print_LDOS_main(LDOS_2d, bs_env, band_edges, scf_gw_soc)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_2d
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(band_edges_type)                              :: band_edges
      CHARACTER(LEN=*)                                   :: scf_gw_soc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'print_LDOS_main'

      INTEGER :: handle, i_x, i_x_bin, i_x_end, i_x_end_bin, i_x_end_glob, i_x_start, &
         i_x_start_bin, i_x_start_glob, i_y, i_y_bin, i_y_end, i_y_end_bin, i_y_end_glob, &
         i_y_start, i_y_start_bin, i_y_start_glob, n_E
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: n_sum_for_bins
      INTEGER, DIMENSION(2)                              :: bin_mesh
      LOGICAL                                            :: do_xy_bins
      REAL(KIND=dp)                                      :: E_min, energy_step, energy_window
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_2d_bins

      CALL timeset(routineN, handle)

      n_E = SIZE(LDOS_2d, 3)

      energy_window = bs_env%energy_window_DOS
      energy_step = bs_env%energy_step_DOS
      E_min = band_edges%VBM - 0.5_dp*energy_window

      bin_mesh(1:2) = bs_env%bin_mesh(1:2)
      do_xy_bins = (bin_mesh(1) > 0 .AND. bin_mesh(2) > 0)

      i_x_start = LBOUND(LDOS_2d, 1)
      i_x_end = UBOUND(LDOS_2d, 1)
      i_y_start = LBOUND(LDOS_2d, 2)
      i_y_end = UBOUND(LDOS_2d, 2)

      IF (do_xy_bins) THEN
         i_x_start_bin = 1
         i_x_end_bin = bin_mesh(1)
         i_y_start_bin = 1
         i_y_end_bin = bin_mesh(2)
      ELSE
         i_x_start_bin = i_x_start
         i_x_end_bin = i_x_end
         i_y_start_bin = i_y_start
         i_y_end_bin = i_y_end
      END IF

      ALLOCATE (LDOS_2d_bins(i_x_start_bin:i_x_end_bin, i_y_start_bin:i_y_end_bin, n_E))
      LDOS_2d_bins(:, :, :) = 0.0_dp

      IF (do_xy_bins) THEN

         i_x_start_glob = i_x_start
         i_x_end_glob = i_x_end
         i_y_start_glob = i_y_start
         i_y_end_glob = i_y_end

         CALL bs_env%para_env%min(i_x_start_glob)
         CALL bs_env%para_env%max(i_x_end_glob)
         CALL bs_env%para_env%min(i_y_start_glob)
         CALL bs_env%para_env%max(i_y_end_glob)

         ALLOCATE (n_sum_for_bins(bin_mesh(1), bin_mesh(2)))
         n_sum_for_bins(:, :) = 0

         ! transform interval [i_x_start, i_x_end] to [1, bin_mesh(1)] (and same for y)
         DO i_x = i_x_start, i_x_end
            DO i_y = i_y_start, i_y_end
               i_x_bin = bin_mesh(1)*(i_x - i_x_start_glob)/(i_x_end_glob - i_x_start_glob + 1) + 1
               i_y_bin = bin_mesh(2)*(i_y - i_y_start_glob)/(i_y_end_glob - i_y_start_glob + 1) + 1
               LDOS_2d_bins(i_x_bin, i_y_bin, :) = LDOS_2d_bins(i_x_bin, i_y_bin, :) + &
                                                   LDOS_2d(i_x, i_y, :)
               n_sum_for_bins(i_x_bin, i_y_bin) = n_sum_for_bins(i_x_bin, i_y_bin) + 1
            END DO
         END DO

         CALL bs_env%para_env%sum(LDOS_2d_bins)
         CALL bs_env%para_env%sum(n_sum_for_bins)

         ! divide by number of terms in the sum so we have the average LDOS(x,y,E)
         DO i_x_bin = 1, bin_mesh(1)
            DO i_y_bin = 1, bin_mesh(2)
               LDOS_2d_bins(i_x_bin, i_y_bin, :) = LDOS_2d_bins(i_x_bin, i_y_bin, :)/ &
                                                   REAL(n_sum_for_bins(i_x_bin, i_y_bin), KIND=dp)
            END DO
         END DO

      ELSE

         LDOS_2d_bins(:, :, :) = LDOS_2d(:, :, :)

      END IF

      IF (bin_mesh(1)*bin_mesh(2) < bs_env%n_bins_max_for_printing) THEN
         CALL print_LDOS_2d_bins(LDOS_2d_bins, bs_env, E_min, scf_gw_soc)
      ELSE
         CPWARN("The number of bins for the LDOS is too large. Decrease BIN_MESH.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_LDOS_main

! **************************************************************************************************
!> \brief ...
!> \param LDOS_2d_bins ...
!> \param bs_env ...
!> \param E_min ...
!> \param scf_gw_soc ...
! **************************************************************************************************
   SUBROUTINE print_LDOS_2d_bins(LDOS_2d_bins, bs_env, E_min, scf_gw_soc)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_2d_bins
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: E_min
      CHARACTER(LEN=*)                                   :: scf_gw_soc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_LDOS_2d_bins'

      CHARACTER(LEN=18)                                  :: print_format
      CHARACTER(LEN=4)                                   :: print_format_1, print_format_2
      CHARACTER(len=default_string_length)               :: fname
      INTEGER                                            :: handle, i_E, i_x, i_x_end, i_x_start, &
                                                            i_y, i_y_end, i_y_start, iunit, n_E, &
                                                            n_x, n_y
      REAL(KIND=dp)                                      :: energy
      REAL(KIND=dp), DIMENSION(3)                        :: coord, idx

      CALL timeset(routineN, handle)

      i_x_start = LBOUND(LDOS_2d_bins, 1)
      i_x_end = UBOUND(LDOS_2d_bins, 1)
      i_y_start = LBOUND(LDOS_2d_bins, 2)
      i_y_end = UBOUND(LDOS_2d_bins, 2)
      n_E = SIZE(LDOS_2d_bins, 3)

      n_x = i_x_end - i_x_start + 1
      n_y = i_y_end - i_y_start + 1

      IF (bs_env%para_env%is_source()) THEN

         DO i_x = i_x_start, i_x_end
            DO i_y = i_y_start, i_y_end

               idx(1) = (REAL(i_x, KIND=dp) - 0.5_dp)/REAL(n_x, KIND=dp)
               idx(2) = (REAL(i_y, KIND=dp) - 0.5_dp)/REAL(n_y, KIND=dp)
               idx(3) = 0.0_dp
               coord(1:3) = MATMUL(bs_env%hmat, idx)

               CALL get_print_format(coord(1), print_format_1)
               CALL get_print_format(coord(2), print_format_2)

               print_format = "(3A,"//print_format_1//",A,"//print_format_2//",A)"

               WRITE (fname, print_format) "LDOS_", scf_gw_soc, &
                  "_at_x_", coord(1)*angstrom, '_A_and_y_', coord(2)*angstrom, '_A'

               CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", &
                              file_action="WRITE")

               WRITE (iunit, "(2A)") "        Energy E (eV)    average LDOS(x,y,E) (1/(eV*Å^2), ", &
                  "integrated over z, averaged inside bin)"

               DO i_E = 1, n_E
                  energy = E_min + i_E*bs_env%energy_step_DOS
                  WRITE (iunit, "(2F17.3)") energy*evolt, &
                     LDOS_2d_bins(i_x, i_y, i_E)* &
                     bs_env%unit_ldos_int_z_inv_Ang2_eV
               END DO

               CALL close_file(iunit)

            END DO
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE print_LDOS_2d_bins

! **************************************************************************************************
!> \brief ...
!> \param coord ...
!> \param print_format ...
! **************************************************************************************************
   SUBROUTINE get_print_format(coord, print_format)
      REAL(KIND=dp)                                      :: coord
      CHARACTER(LEN=4)                                   :: print_format

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_print_format'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (coord < -10000/angstrom) THEN
         print_format = "F9.2"
      ELSE IF (coord < -1000/angstrom) THEN
         print_format = "F8.2"
      ELSE IF (coord < -100/angstrom) THEN
         print_format = "F7.2"
      ELSE IF (coord < -10/angstrom) THEN
         print_format = "F6.2"
      ELSE IF (coord < -1/angstrom) THEN
         print_format = "F5.2"
      ELSE IF (coord < 10/angstrom) THEN
         print_format = "F4.2"
      ELSE IF (coord < 100/angstrom) THEN
         print_format = "F5.2"
      ELSE IF (coord < 1000/angstrom) THEN
         print_format = "F6.2"
      ELSE IF (coord < 10000/angstrom) THEN
         print_format = "F7.2"
      ELSE
         print_format = "F8.2"
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_print_format

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param ikp ...
!> \param eigenval_no_SOC ...
!> \param band_edges_no_SOC ...
!> \param E_min ...
!> \param cfm_mos_ikp ...
!> \param DOS ...
!> \param PDOS ...
!> \param band_edges ...
!> \param eigenval_spinor ...
!> \param cfm_spinor_wf_ikp ...
! **************************************************************************************************
   SUBROUTINE SOC_ev(bs_env, qs_env, ikp, eigenval_no_SOC, band_edges_no_SOC, E_min, cfm_mos_ikp, &
                     DOS, PDOS, band_edges, eigenval_spinor, cfm_spinor_wf_ikp)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: ikp
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: eigenval_no_SOC
      TYPE(band_edges_type)                              :: band_edges_no_SOC
      REAL(KIND=dp)                                      :: E_min
      TYPE(cp_cfm_type), DIMENSION(2)                    :: cfm_mos_ikp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: DOS
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: PDOS
      TYPE(band_edges_type)                              :: band_edges
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_spinor
      TYPE(cp_cfm_type)                                  :: cfm_spinor_wf_ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'SOC_ev'

      INTEGER                                            :: handle, homo_spinor, n_ao, n_E, nkind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_spinor_no_SOC
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: proj_mo_on_kind_spinor
      TYPE(cp_cfm_type)                                  :: cfm_eigenvec_ikp_spinor, &
                                                            cfm_ks_ikp_spinor, cfm_mos_ikp_spinor, &
                                                            cfm_SOC_ikp_spinor, cfm_work_ikp_spinor

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      homo_spinor = bs_env%n_occ(1) + bs_env%n_occ(bs_env%n_spin)
      n_E = SIZE(DOS)
      nkind = SIZE(PDOS, 2)

      CALL cp_cfm_create(cfm_ks_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
      CALL cp_cfm_create(cfm_SOC_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
      CALL cp_cfm_create(cfm_mos_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
      CALL cp_cfm_create(cfm_work_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)
      CALL cp_cfm_create(cfm_eigenvec_ikp_spinor, bs_env%cfm_SOC_spinor_ao(1)%matrix_struct)

      ALLOCATE (eigenval_spinor_no_SOC(2*n_ao))
      ALLOCATE (proj_mo_on_kind_spinor(2*n_ao, nkind))
      ! PDOS not yet implemented -> projection is just zero -> PDOS is zero
      proj_mo_on_kind_spinor(:, :) = 0.0_dp

      ! 1. get V^SOC_µν,σσ'(k_i)
      SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
      CASE (large_cell_Gamma)

         ! 1. get V^SOC_µν,σσ'(k_i) from V^SOC_µν,σσ'(k=0)
         CALL cfm_ikp_from_cfm_spinor_Gamma(cfm_SOC_ikp_spinor, &
                                            bs_env%cfm_SOC_spinor_ao(1), &
                                            bs_env%fm_s_Gamma%matrix_struct, &
                                            ikp, qs_env, bs_env%kpoints_DOS, "ORB")

      CASE (small_cell_full_kp)

         ! 1. V^SOC_µν,σσ'(k_i) already there
         CALL cp_cfm_to_cfm(bs_env%cfm_SOC_spinor_ao(ikp), cfm_SOC_ikp_spinor)

      END SELECT

      ! 2. V^SOC_nn',σσ'(k_i) = sum_µν C^*_µn,σ(k_i) V^SOC_µν,σσ'(k_i) C_νn'(k_i),
      !    C_µn,σ(k_i): MO coefficiencts from diagonalizing KS-matrix h^KS_nn',σσ'(k_i)

      ! 2.1 build matrix C_µn,σ(k_i)
      CALL cp_cfm_set_all(cfm_mos_ikp_spinor, z_zero)
      CALL add_cfm_submat(cfm_mos_ikp_spinor, cfm_mos_ikp(1), 1, 1)
      CALL add_cfm_submat(cfm_mos_ikp_spinor, cfm_mos_ikp(bs_env%n_spin), n_ao + 1, n_ao + 1)

      ! 2.2 work_nν,σσ' = sum_µ C^*_µn,σ(k_i) V^SOC_µν,σσ'(k_i)
      CALL parallel_gemm('C', 'N', 2*n_ao, 2*n_ao, 2*n_ao, z_one, &
                         cfm_mos_ikp_spinor, cfm_SOC_ikp_spinor, &
                         z_zero, cfm_work_ikp_spinor)

      ! 2.3 V^SOC_nn',σσ'(k_i) = sum_ν work_nν,σσ' C_νn'(k_i)
      CALL parallel_gemm('N', 'N', 2*n_ao, 2*n_ao, 2*n_ao, z_one, &
                         cfm_work_ikp_spinor, cfm_mos_ikp_spinor, &
                         z_zero, cfm_ks_ikp_spinor)

      ! 3. remove SOC outside of energy window (otherwise, numerical problems arise
      !    because energetically low semicore states and energetically very high
      !    unbound states couple to the states around the Fermi level)
      eigenval_spinor_no_SOC(1:n_ao) = eigenval_no_SOC(1:n_ao, ikp, 1)
      eigenval_spinor_no_SOC(n_ao + 1:) = eigenval_no_SOC(1:n_ao, ikp, bs_env%n_spin)
      IF (bs_env%energy_window_soc > 0.0_dp) THEN
         CALL remove_soc_outside_energy_window_mo(cfm_ks_ikp_spinor, &
                                                  bs_env%energy_window_soc, &
                                                  eigenval_spinor_no_SOC, &
                                                  band_edges_no_SOC%VBM, &
                                                  band_edges_no_SOC%CBM)
      END IF

      ! 4. h^G0W0+SOC_nn',σσ'(k_i) = ε_nσ^G0W0(k_i) δ_nn' δ_σσ' + V^SOC_nn',σσ'(k_i)
      CALL cfm_add_on_diag(cfm_ks_ikp_spinor, eigenval_spinor_no_SOC)

      ! 5. diagonalize h^G0W0+SOC_nn',σσ'(k_i) to get eigenvalues
      CALL cp_cfm_heevd(cfm_ks_ikp_spinor, cfm_eigenvec_ikp_spinor, eigenval_spinor)

      ! 6. DOS from spinors, no PDOS
      CALL add_to_DOS_PDOS(DOS, PDOS, eigenval_spinor, &
                           ikp, bs_env, n_E, E_min, proj_mo_on_kind_spinor)

      ! 7. valence band max. (VBM), conduction band min. (CBM) and direct bandgap (DBG)
      band_edges%VBM = MAX(band_edges%VBM, eigenval_spinor(homo_spinor))
      band_edges%CBM = MIN(band_edges%CBM, eigenval_spinor(homo_spinor + 1))
      band_edges%DBG = MIN(band_edges%DBG, eigenval_spinor(homo_spinor + 1) &
                           - eigenval_spinor(homo_spinor))

      ! 8. spinor wavefunctions:
      CALL parallel_gemm('N', 'N', 2*n_ao, 2*n_ao, 2*n_ao, z_one, &
                         cfm_mos_ikp_spinor, cfm_eigenvec_ikp_spinor, &
                         z_zero, cfm_spinor_wf_ikp)

      CALL cp_cfm_release(cfm_ks_ikp_spinor)
      CALL cp_cfm_release(cfm_SOC_ikp_spinor)
      CALL cp_cfm_release(cfm_work_ikp_spinor)
      CALL cp_cfm_release(cfm_eigenvec_ikp_spinor)
      CALL cp_cfm_release(cfm_mos_ikp_spinor)

      CALL timestop(handle)

   END SUBROUTINE SOC_ev

! **************************************************************************************************
!> \brief ...
!> \param DOS ...
!> \param PDOS ...
!> \param eigenval ...
!> \param ikp ...
!> \param bs_env ...
!> \param n_E ...
!> \param E_min ...
!> \param proj_mo_on_kind ...
! **************************************************************************************************
   SUBROUTINE add_to_DOS_PDOS(DOS, PDOS, eigenval, ikp, bs_env, n_E, E_min, proj_mo_on_kind)

      REAL(KIND=dp), DIMENSION(:)                        :: DOS
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: PDOS
      REAL(KIND=dp), DIMENSION(:)                        :: eigenval
      INTEGER                                            :: ikp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: n_E
      REAL(KIND=dp)                                      :: E_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: proj_mo_on_kind

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_to_DOS_PDOS'

      INTEGER                                            :: handle, i_E, i_kind, i_mo, n_mo, nkind
      REAL(KIND=dp)                                      :: broadening, energy, energy_step_DOS, wkp

      CALL timeset(routineN, handle)

      energy_step_DOS = bs_env%energy_step_DOS
      broadening = bs_env%broadening_DOS

      n_mo = SIZE(eigenval)
      nkind = SIZE(proj_mo_on_kind, 2)

      ! normalize to closed-shell / open-shell
      wkp = bs_env%kpoints_DOS%wkp(ikp)*bs_env%spin_degeneracy
      DO i_E = 1, n_E
         energy = E_min + i_E*energy_step_DOS
         DO i_mo = 1, n_mo
            ! DOS
            DOS(i_E) = DOS(i_E) + wkp*Gaussian(energy - eigenval(i_mo), broadening)

            ! PDOS
            DO i_kind = 1, nkind
               IF (proj_mo_on_kind(i_mo, i_kind) > 0.0_dp) THEN
                  PDOS(i_E, i_kind) = PDOS(i_E, i_kind) + &
                                      proj_mo_on_kind(i_mo, i_kind)*wkp* &
                                      Gaussian(energy - eigenval(i_mo), broadening)
               END IF
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE add_to_DOS_PDOS

! **************************************************************************************************
!> \brief ...
!> \param LDOS_2d ...
!> \param qs_env ...
!> \param ikp ...
!> \param bs_env ...
!> \param cfm_mos_ikp ...
!> \param eigenval ...
!> \param band_edges ...
!> \param do_spinor ...
!> \param cfm_non_spinor ...
! **************************************************************************************************
   SUBROUTINE add_to_LDOS_2d(LDOS_2d, qs_env, ikp, bs_env, cfm_mos_ikp, eigenval, &
                             band_edges, do_spinor, cfm_non_spinor)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_2d
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: ikp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_mos_ikp
      REAL(KIND=dp), DIMENSION(:)                        :: eigenval
      TYPE(band_edges_type)                              :: band_edges
      LOGICAL, OPTIONAL                                  :: do_spinor
      TYPE(cp_cfm_type), OPTIONAL                        :: cfm_non_spinor

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_to_LDOS_2d'

      INTEGER :: handle, i_E, i_x_end, i_x_start, i_y_end, i_y_start, i_z, i_z_end, i_z_start, &
         j_col, j_mo, n_E, n_mo, n_z, ncol_local, nimages, z_end_global, z_start_global
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      LOGICAL                                            :: is_any_weight_non_zero, my_do_spinor
      REAL(KIND=dp)                                      :: broadening, E_max, E_min, &
                                                            E_total_window, energy, energy_step, &
                                                            energy_window, spin_degeneracy, weight
      TYPE(cp_cfm_type)                                  :: cfm_weighted_dm_ikp, cfm_work
      TYPE(cp_fm_type)                                   :: fm_non_spinor, fm_weighted_dm_MIC
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: weighted_dm_MIC
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type)                               :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: LDOS_3d
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      my_do_spinor = .FALSE.
      IF (PRESENT(do_spinor)) my_do_spinor = do_spinor

      CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env, dft_control=dft_control)

      ! previously, dft_control%nimages set to # neighbor cells, revert for Γ-only KS matrix
      nimages = dft_control%nimages
      dft_control%nimages = bs_env%nimages_scf

      energy_window = bs_env%energy_window_DOS
      energy_step = bs_env%energy_step_DOS
      broadening = bs_env%broadening_DOS

      E_min = band_edges%VBM - 0.5_dp*energy_window
      E_max = band_edges%CBM + 0.5_dp*energy_window
      E_total_window = E_max - E_min

      n_E = INT(E_total_window/energy_step)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL auxbas_pw_pool%create_pw(LDOS_3d)
      CALL auxbas_pw_pool%create_pw(rho_g)

      i_x_start = LBOUND(LDOS_3d%array, 1)
      i_x_end = UBOUND(LDOS_3d%array, 1)
      i_y_start = LBOUND(LDOS_3d%array, 2)
      i_y_end = UBOUND(LDOS_3d%array, 2)
      i_z_start = LBOUND(LDOS_3d%array, 3)
      i_z_end = UBOUND(LDOS_3d%array, 3)

      z_start_global = i_z_start
      z_end_global = i_z_end

      CALL bs_env%para_env%min(z_start_global)
      CALL bs_env%para_env%max(z_end_global)
      n_z = z_end_global - z_start_global + 1

      IF (ANY(ABS(bs_env%hmat(1:2, 3)) > 1.0E-6_dp) .OR. ANY(ABS(bs_env%hmat(3, 1:2)) > 1.0E-6_dp)) &
         CPABORT("Please choose a cell that has 90° angles to the z-direction.")
      ! for integration, we need the dz and the conversion from H -> eV and a_Bohr -> Å
      bs_env%unit_ldos_int_z_inv_Ang2_eV = bs_env%hmat(3, 3)/REAL(n_z, KIND=dp)/evolt/angstrom**2

      IF (ikp == 1) THEN
         ALLOCATE (LDOS_2d(i_x_start:i_x_end, i_y_start:i_y_end, n_E))
         LDOS_2d(:, :, :) = 0.0_dp
      END IF

      CALL cp_cfm_create(cfm_work, cfm_mos_ikp%matrix_struct)
      CALL cp_cfm_create(cfm_weighted_dm_ikp, cfm_mos_ikp%matrix_struct)
      CALL cp_fm_create(fm_weighted_dm_MIC, cfm_mos_ikp%matrix_struct)
      IF (my_do_spinor) THEN
         CALL cp_fm_create(fm_non_spinor, cfm_non_spinor%matrix_struct)
      END IF

      CALL cp_cfm_get_info(matrix=cfm_mos_ikp, &
                           ncol_global=n_mo, &
                           ncol_local=ncol_local, &
                           col_indices=col_indices)

      NULLIFY (weighted_dm_MIC)
      CALL dbcsr_allocate_matrix_set(weighted_dm_MIC, 1)
      ALLOCATE (weighted_dm_MIC(1)%matrix)
      CALL dbcsr_create(weighted_dm_MIC(1)%matrix, template=bs_env%mat_ao_ao%matrix, &
                        matrix_type=dbcsr_type_symmetric)

      DO i_E = 1, n_E

         energy = E_min + i_E*energy_step

         is_any_weight_non_zero = .FALSE.

         DO j_col = 1, ncol_local

            j_mo = col_indices(j_col)

            IF (my_do_spinor) THEN
               spin_degeneracy = 1.0_dp
            ELSE
               spin_degeneracy = bs_env%spin_degeneracy
            END IF

            weight = Gaussian(energy - eigenval(j_mo), broadening)*spin_degeneracy

            cfm_work%local_data(:, j_col) = cfm_mos_ikp%local_data(:, j_col)*weight

            IF (weight > 1.0E-5_dp) is_any_weight_non_zero = .TRUE.

         END DO

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(is_any_weight_non_zero)
         CALL bs_env%para_env%sync()

         ! cycle if there are no states at the energy i_E
         IF (is_any_weight_non_zero) THEN

            CALL parallel_gemm('N', 'C', n_mo, n_mo, n_mo, z_one, &
                               cfm_mos_ikp, cfm_work, z_zero, cfm_weighted_dm_ikp)

            IF (my_do_spinor) THEN

               ! contribution from up,up to fm_non_spinor
               CALL get_cfm_submat(cfm_non_spinor, cfm_weighted_dm_ikp, 1, 1)
               CALL cp_fm_set_all(fm_non_spinor, 0.0_dp)
               CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_non_spinor, &
                                              cfm_non_spinor, ikp, bs_env%kpoints_DOS, &
                                              "ORB", bs_env%kpoints_DOS%wkp(ikp))

               ! add contribution from down,down to fm_non_spinor
               CALL get_cfm_submat(cfm_non_spinor, cfm_weighted_dm_ikp, n_mo/2, n_mo/2)
               CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_non_spinor, &
                                              cfm_non_spinor, ikp, bs_env%kpoints_DOS, &
                                              "ORB", bs_env%kpoints_DOS%wkp(ikp))
               CALL copy_fm_to_dbcsr(fm_non_spinor, weighted_dm_MIC(1)%matrix, &
                                     keep_sparsity=.FALSE.)
            ELSE
               CALL cp_fm_set_all(fm_weighted_dm_MIC, 0.0_dp)
               CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_weighted_dm_MIC, &
                                              cfm_weighted_dm_ikp, ikp, bs_env%kpoints_DOS, &
                                              "ORB", bs_env%kpoints_DOS%wkp(ikp))
               CALL copy_fm_to_dbcsr(fm_weighted_dm_MIC, weighted_dm_MIC(1)%matrix, &
                                     keep_sparsity=.FALSE.)
            END IF

            LDOS_3d%array(:, :, :) = 0.0_dp

            CALL calculate_rho_elec(matrix_p_kp=weighted_dm_MIC, &
                                    rho=LDOS_3d, &
                                    rho_gspace=rho_g, &
                                    ks_env=ks_env)

            DO i_z = i_z_start, i_z_end
               LDOS_2d(:, :, i_E) = LDOS_2d(:, :, i_E) + LDOS_3d%array(:, :, i_z)
            END DO

         END IF

      END DO

      ! set back nimages
      dft_control%nimages = nimages

      CALL auxbas_pw_pool%give_back_pw(LDOS_3d)
      CALL auxbas_pw_pool%give_back_pw(rho_g)

      CALL cp_cfm_release(cfm_work)
      CALL cp_cfm_release(cfm_weighted_dm_ikp)

      CALL cp_fm_release(fm_weighted_dm_MIC)

      CALL dbcsr_deallocate_matrix_set(weighted_dm_MIC)

      IF (my_do_spinor) THEN
         CALL cp_fm_release(fm_non_spinor)
      END IF

      CALL timestop(handle)

   END SUBROUTINE add_to_LDOS_2d

! **************************************************************************************************
!> \brief ...
!> \param eigenval_spinor ...
!> \param ikp_for_file ...
!> \param ikp ...
!> \param bs_env ...
!> \param eigenval_spinor_G0W0 ...
! **************************************************************************************************
   SUBROUTINE write_SOC_eigenvalues(eigenval_spinor, ikp_for_file, ikp, bs_env, eigenval_spinor_G0W0)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_spinor
      INTEGER                                            :: ikp_for_file, ikp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), OPTIONAL :: eigenval_spinor_G0W0

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_SOC_eigenvalues'

      CHARACTER(len=3)                                   :: occ_vir
      CHARACTER(LEN=default_string_length)               :: fname
      INTEGER                                            :: handle, i_mo, iunit, n_occ_spinor

      CALL timeset(routineN, handle)

      fname = "bandstructure_SCF_and_G0W0_plus_SOC"

      IF (bs_env%para_env%is_source()) THEN

         IF (ikp_for_file == 1) THEN
            CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", &
                           file_action="WRITE")
         ELSE
            CALL open_file(TRIM(fname), unit_number=iunit, file_status="OLD", &
                           file_action="WRITE", file_position="APPEND")
         END IF

         WRITE (iunit, "(A)") " "
         WRITE (iunit, "(A10,I7,A25,3F10.4)") "kpoint: ", ikp_for_file, "coordinate: ", &
            bs_env%kpoints_DOS%xkp(:, ikp)
         WRITE (iunit, "(A)") " "

         IF (PRESENT(eigenval_spinor_G0W0)) THEN
            ! SCF+SOC and G0W0+SOC eigenvalues
            WRITE (iunit, "(A5,A12,2A22)") "n", "k", "ϵ_nk^DFT+SOC (eV)", "ϵ_nk^G0W0+SOC (eV)"
         ELSE
            ! SCF+SOC eigenvalues only
            WRITE (iunit, "(A5,A12,A22)") "n", "k", "ϵ_nk^DFT+SOC (eV)"
         END IF

         n_occ_spinor = bs_env%n_occ(1) + bs_env%n_occ(bs_env%n_spin)

         DO i_mo = 1, SIZE(eigenval_spinor)
            IF (i_mo .LE. n_occ_spinor) occ_vir = 'occ'
            IF (i_mo > n_occ_spinor) occ_vir = 'vir'
            IF (PRESENT(eigenval_spinor_G0W0)) THEN
               ! SCF+SOC and G0W0+SOC eigenvalues
               WRITE (iunit, "(I5,3A,I5,4F16.3,2F17.3)") i_mo, ' (', occ_vir, ') ', &
                  ikp_for_file, eigenval_spinor(i_mo)*evolt, eigenval_spinor_G0W0(i_mo)*evolt
            ELSE
               ! SCF+SOC eigenvalues only
               WRITE (iunit, "(I5,3A,I5,4F16.3,F17.3)") i_mo, ' (', occ_vir, ') ', &
                  ikp_for_file, eigenval_spinor(i_mo)*evolt
            END IF
         END DO

         CALL close_file(iunit)

      END IF

      CALL timestop(handle)

   END SUBROUTINE write_SOC_eigenvalues

! **************************************************************************************************
!> \brief ...
!> \param int_number ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION count_digits(int_number)

      INTEGER, INTENT(IN)                                :: int_number
      INTEGER                                            :: count_digits

      INTEGER                                            :: digitCount, tempInt

      digitCount = 0

      tempInt = int_number

      DO WHILE (tempInt /= 0)
         tempInt = tempInt/10
         digitCount = digitCount + 1
      END DO

      count_digits = digitCount

   END FUNCTION count_digits

! **************************************************************************************************
!> \brief ...
!> \param band_edges ...
!> \param scf_gw_soc ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE write_band_edges(band_edges, scf_gw_soc, bs_env)

      TYPE(band_edges_type)                              :: band_edges
      CHARACTER(LEN=*)                                   :: scf_gw_soc
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_band_edges'

      CHARACTER(LEN=17)                                  :: print_format
      INTEGER                                            :: handle, u

      CALL timeset(routineN, handle)

      ! print format
      print_format = "(T2,2A,T61,F20.3)"

      u = bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, print_format) scf_gw_soc, ' valence band maximum (eV):', band_edges%VBM*evolt
         WRITE (u, print_format) scf_gw_soc, ' conduction band minimum (eV):', band_edges%CBM*evolt
         WRITE (u, print_format) scf_gw_soc, ' indirect band gap (eV):', band_edges%IDBG*evolt
         WRITE (u, print_format) scf_gw_soc, ' direct band gap (eV):', band_edges%DBG*evolt
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_band_edges

! **************************************************************************************************
!> \brief ...
!> \param DOS ...
!> \param PDOS ...
!> \param bs_env ...
!> \param qs_env ...
!> \param scf_gw_soc ...
!> \param E_min ...
!> \param E_VBM ...
! **************************************************************************************************
   SUBROUTINE write_dos_pdos(DOS, PDOS, bs_env, qs_env, scf_gw_soc, E_min, E_VBM)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: DOS
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: PDOS
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*)                                   :: scf_gw_soc
      REAL(KIND=dp)                                      :: E_min, E_VBM

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_dos_pdos'

      CHARACTER(LEN=3), DIMENSION(100)                   :: elements
      CHARACTER(LEN=default_string_length)               :: atom_name, fname, output_string
      INTEGER                                            :: handle, i_E, i_kind, iatom, iunit, n_A, &
                                                            n_E, nkind
      REAL(KIND=dp)                                      :: energy
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      WRITE (fname, "(3A)") "DOS_PDOS_", scf_gw_soc, ".out"

      n_E = SIZE(PDOS, 1)
      nkind = SIZE(PDOS, 2)
      CALL get_qs_env(qs_env, particle_set=particle_set)

      IF (bs_env%para_env%is_source()) THEN

         CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", file_action="WRITE")

         n_A = 2 + nkind

         DO iatom = 1, bs_env%n_atom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 kind_number=i_kind, name=atom_name)
            elements(i_kind) = atom_name
         END DO

         WRITE (output_string, "(A,I1,A)") "(", n_A, "A)"

         WRITE (iunit, TRIM(output_string)) "Energy-E_F (eV)    DOS (1/eV)    PDOS (1/eV) ", &
            " of atom type ", elements(1:nkind)

         WRITE (output_string, "(A,I1,A)") "(", n_A, "F13.5)"

         DO i_E = 1, n_E
            ! energy is relative to valence band maximum => - E_VBM
            energy = E_min + i_E*bs_env%energy_step_DOS - E_VBM
            WRITE (iunit, TRIM(output_string)) energy*evolt, DOS(i_E)/evolt, PDOS(i_E, :)/evolt
         END DO

         CALL close_file(iunit)

      END IF

      CALL timestop(handle)

   END SUBROUTINE write_dos_pdos

! **************************************************************************************************
!> \brief ...
!> \param energy ...
!> \param broadening ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION Gaussian(energy, broadening)

      REAL(KIND=dp), INTENT(IN)                          :: energy, broadening
      REAL(KIND=dp)                                      :: Gaussian

      IF (ABS(energy) < 5*broadening) THEN
         Gaussian = 1.0_dp/broadening/SQRT(twopi)*EXP(-0.5_dp*energy**2/broadening**2)
      ELSE
         Gaussian = 0.0_dp
      END IF

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param proj_mo_on_kind ...
!> \param qs_env ...
!> \param cfm_mos ...
!> \param cfm_s ...
! **************************************************************************************************
   SUBROUTINE compute_proj_mo_on_kind(proj_mo_on_kind, qs_env, cfm_mos, cfm_s)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: proj_mo_on_kind
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_type)                                  :: cfm_mos, cfm_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_proj_mo_on_kind'

      INTEGER                                            :: handle, i_atom, i_global, i_kind, i_row, &
                                                            j_col, n_ao, n_mo, ncol_local, nkind, &
                                                            nrow_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf, kind_of
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_cfm_type)                                  :: cfm_proj, cfm_s_i_kind, cfm_work
      TYPE(cp_fm_type)                                   :: fm_proj_im, fm_proj_re

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, nkind=nkind)
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      CALL cp_cfm_get_info(matrix=cfm_mos, &
                           nrow_global=n_mo, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      n_ao = qs_env%bs_env%n_ao

      ALLOCATE (atom_from_bf(n_ao))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf, n_ao, "ORB")

      proj_mo_on_kind(:, :) = 0.0_dp

      CALL cp_cfm_create(cfm_s_i_kind, cfm_s%matrix_struct)
      CALL cp_cfm_create(cfm_work, cfm_s%matrix_struct)
      CALL cp_cfm_create(cfm_proj, cfm_s%matrix_struct)
      CALL cp_fm_create(fm_proj_re, cfm_s%matrix_struct)
      CALL cp_fm_create(fm_proj_im, cfm_s%matrix_struct)

      DO i_kind = 1, nkind

         CALL cp_cfm_to_cfm(cfm_s, cfm_s_i_kind)

         ! set entries in overlap matrix to zero which do not belong to atoms of i_kind
         DO i_row = 1, nrow_local
            DO j_col = 1, ncol_local

               i_global = row_indices(i_row)

               IF (i_global .LE. n_ao) THEN
                  i_atom = atom_from_bf(i_global)
               ELSE IF (i_global .LE. 2*n_ao) THEN
                  i_atom = atom_from_bf(i_global - n_ao)
               ELSE
                  CPABORT("Wrong indices.")
               END IF

               IF (i_kind .NE. kind_of(i_atom)) THEN
                  cfm_s_i_kind%local_data(i_row, j_col) = z_zero
               END IF

            END DO
         END DO

         CALL parallel_gemm('N', 'N', n_mo, n_mo, n_mo, z_one, &
                            cfm_s_i_kind, cfm_mos, z_zero, cfm_work)
         CALL parallel_gemm('C', 'N', n_mo, n_mo, n_mo, z_one, &
                            cfm_mos, cfm_work, z_zero, cfm_proj)

         CALL cp_cfm_to_fm(cfm_proj, fm_proj_re, fm_proj_im)

         CALL cp_fm_get_diag(fm_proj_im, proj_mo_on_kind(:, i_kind))
         CALL cp_fm_get_diag(fm_proj_re, proj_mo_on_kind(:, i_kind))

      END DO ! i_kind

      CALL cp_cfm_release(cfm_s_i_kind)
      CALL cp_cfm_release(cfm_work)
      CALL cp_cfm_release(cfm_proj)
      CALL cp_fm_release(fm_proj_re)
      CALL cp_fm_release(fm_proj_im)

      CALL timestop(handle)

   END SUBROUTINE compute_proj_mo_on_kind

! **************************************************************************************************
!> \brief ...
!> \param cfm_spinor_ikp ...
!> \param cfm_spinor_Gamma ...
!> \param fm_struct_non_spinor ...
!> \param ikp ...
!> \param qs_env ...
!> \param kpoints ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE cfm_ikp_from_cfm_spinor_Gamma(cfm_spinor_ikp, cfm_spinor_Gamma, fm_struct_non_spinor, &
                                            ikp, qs_env, kpoints, basis_type)
      TYPE(cp_cfm_type)                                  :: cfm_spinor_ikp, cfm_spinor_Gamma
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_non_spinor
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*)                                   :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cfm_ikp_from_cfm_spinor_Gamma'

      INTEGER                                            :: handle, i_block, i_offset, j_block, &
                                                            j_offset, n_ao
      TYPE(cp_cfm_type)                                  :: cfm_non_spinor_Gamma, cfm_non_spinor_ikp
      TYPE(cp_fm_type)                                   :: fm_non_spinor_Gamma_im, &
                                                            fm_non_spinor_Gamma_re

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_non_spinor_Gamma, fm_struct_non_spinor)
      CALL cp_cfm_create(cfm_non_spinor_ikp, fm_struct_non_spinor)
      CALL cp_fm_create(fm_non_spinor_Gamma_re, fm_struct_non_spinor)
      CALL cp_fm_create(fm_non_spinor_Gamma_im, fm_struct_non_spinor)

      CALL cp_cfm_get_info(cfm_non_spinor_Gamma, nrow_global=n_ao)

      CALL cp_cfm_set_all(cfm_spinor_ikp, z_zero)

      DO i_block = 0, 1
         DO j_block = 0, 1
            i_offset = i_block*n_ao + 1
            j_offset = j_block*n_ao + 1
            CALL get_cfm_submat(cfm_non_spinor_Gamma, cfm_spinor_Gamma, i_offset, j_offset)
            CALL cp_cfm_to_fm(cfm_non_spinor_Gamma, fm_non_spinor_Gamma_re, fm_non_spinor_Gamma_im)

            ! transform real part of Gamma-point matrix to ikp
            CALL cfm_ikp_from_fm_Gamma(cfm_non_spinor_ikp, fm_non_spinor_Gamma_re, &
                                       ikp, qs_env, kpoints, basis_type)
            CALL add_cfm_submat(cfm_spinor_ikp, cfm_non_spinor_ikp, i_offset, j_offset)

            ! transform imag part of Gamma-point matrix to ikp
            CALL cfm_ikp_from_fm_Gamma(cfm_non_spinor_ikp, fm_non_spinor_Gamma_im, &
                                       ikp, qs_env, kpoints, basis_type)
            CALL add_cfm_submat(cfm_spinor_ikp, cfm_non_spinor_ikp, i_offset, j_offset, gaussi)

         END DO
      END DO

      CALL cp_cfm_release(cfm_non_spinor_Gamma)
      CALL cp_cfm_release(cfm_non_spinor_ikp)
      CALL cp_fm_release(fm_non_spinor_Gamma_re)
      CALL cp_fm_release(fm_non_spinor_Gamma_im)

      CALL timestop(handle)

   END SUBROUTINE cfm_ikp_from_cfm_spinor_Gamma

! **************************************************************************************************
!> \brief ...
!> \param cfm_ikp ...
!> \param fm_Gamma ...
!> \param ikp ...
!> \param qs_env ...
!> \param kpoints ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE cfm_ikp_from_fm_Gamma(cfm_ikp, fm_Gamma, ikp, qs_env, kpoints, basis_type)
      TYPE(cp_cfm_type)                                  :: cfm_ikp
      TYPE(cp_fm_type)                                   :: fm_Gamma
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*)                                   :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cfm_ikp_from_fm_Gamma'

      INTEGER :: col_global, handle, i, i_atom, i_atom_old, i_cell, i_mic_cell, i_row, j, j_atom, &
         j_atom_old, j_cell, j_col, n_bf, ncol_local, nrow_local, num_cells, row_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL :: i_cell_is_the_minimum_image_cell
      REAL(KIND=dp)                                      :: abs_rab_cell_i, abs_rab_cell_j, arg
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, cell_vector_j, rab_cell_i, &
                                                            rab_cell_j
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      IF (.NOT. ASSOCIATED(cfm_ikp%local_data)) THEN
         CALL cp_cfm_create(cfm_ikp, fm_Gamma%matrix_struct)
      END IF
      CALL cp_cfm_set_all(cfm_ikp, z_zero)

      CALL cp_fm_get_info(matrix=fm_Gamma, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! get number of basis functions (bf) for different basis sets
      IF (basis_type == "ORB") THEN
         n_bf = qs_env%bs_env%n_ao
      ELSE IF (basis_type == "RI_AUX") THEN
         n_bf = qs_env%bs_env%n_RI
      ELSE
         CPABORT("Only ORB and RI_AUX basis implemented.")
      END IF

      ALLOCATE (atom_from_bf(n_bf))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf, n_bf, basis_type)

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      index_to_cell => kpoints%index_to_cell

      num_cells = SIZE(index_to_cell, 2)
      i_atom_old = 0
      j_atom_old = 0

      DO i_row = 1, nrow_local
         DO j_col = 1, ncol_local

            row_global = row_indices(i_row)
            col_global = col_indices(j_col)

            i_atom = atom_from_bf(row_global)
            j_atom = atom_from_bf(col_global)

            ! we only need to check for new MIC cell for new i_atom-j_atom pair
            IF (i_atom .NE. i_atom_old .OR. j_atom .NE. j_atom_old) THEN
               DO i_cell = 1, num_cells

                  ! only check nearest neigbors
                  IF (ANY(ABS(index_to_cell(1:3, i_cell)) > 1)) CYCLE

                  cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, i_cell), dp))

                  rab_cell_i(1:3) = pbc(particle_set(i_atom)%r(1:3), cell) - &
                                    (pbc(particle_set(j_atom)%r(1:3), cell) + cell_vector(1:3))
                  abs_rab_cell_i = SQRT(rab_cell_i(1)**2 + rab_cell_i(2)**2 + rab_cell_i(3)**2)

                  ! minimum image convention
                  i_cell_is_the_minimum_image_cell = .TRUE.
                  DO j_cell = 1, num_cells
                     cell_vector_j(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, j_cell), dp))
                     rab_cell_j(1:3) = pbc(particle_set(i_atom)%r(1:3), cell) - &
                                       (pbc(particle_set(j_atom)%r(1:3), cell) + cell_vector_j(1:3))
                     abs_rab_cell_j = SQRT(rab_cell_j(1)**2 + rab_cell_j(2)**2 + rab_cell_j(3)**2)

                     IF (abs_rab_cell_i > abs_rab_cell_j + 1.0E-6_dp) THEN
                        i_cell_is_the_minimum_image_cell = .FALSE.
                     END IF
                  END DO

                  IF (i_cell_is_the_minimum_image_cell) THEN
                     i_mic_cell = i_cell
                  END IF

               END DO ! i_cell
            END IF

            arg = REAL(index_to_cell(1, i_mic_cell), dp)*kpoints%xkp(1, ikp) + &
                  REAL(index_to_cell(2, i_mic_cell), dp)*kpoints%xkp(2, ikp) + &
                  REAL(index_to_cell(3, i_mic_cell), dp)*kpoints%xkp(3, ikp)

            i = i_row
            j = j_col

            cfm_ikp%local_data(i, j) = COS(twopi*arg)*fm_Gamma%local_data(i, j)*z_one + &
                                       SIN(twopi*arg)*fm_Gamma%local_data(i, j)*gaussi

            j_atom_old = j_atom
            i_atom_old = i_atom

         END DO ! j_col
      END DO ! i_row

      CALL timestop(handle)

   END SUBROUTINE cfm_ikp_from_fm_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_j ...
!> \param cfm_W_ikp_freq_j ...
!> \param ikp ...
!> \param kpoints ...
!> \param basis_type ...
!> \param wkp_ext ...
! **************************************************************************************************
   SUBROUTINE MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, &
                                        cfm_W_ikp_freq_j, ikp, kpoints, basis_type, wkp_ext)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      TYPE(cp_cfm_type)                                  :: cfm_W_ikp_freq_j
      INTEGER                                            :: ikp
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*)                                   :: basis_type
      REAL(KIND=dp), OPTIONAL                            :: wkp_ext

      CHARACTER(LEN=*), PARAMETER :: routineN = 'MIC_contribution_from_ikp'

      INTEGER                                            :: handle, i_bf, iatom, iatom_old, irow, &
                                                            j_bf, jatom, jatom_old, jcol, n_bf, &
                                                            ncol_local, nrow_local, num_cells
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp)                                      :: contribution, weight_im, weight_re, &
                                                            wkp_of_ikp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      ! get number of basis functions (bf) for different basis sets
      IF (basis_type == "ORB") THEN
         n_bf = qs_env%bs_env%n_ao
      ELSE IF (basis_type == "RI_AUX") THEN
         n_bf = qs_env%bs_env%n_RI
      ELSE
         CPABORT("Only ORB and RI_AUX basis implemented.")
      END IF

      ALLOCATE (atom_from_bf_index(n_bf))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf_index, n_bf, basis_type)

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      CALL cp_cfm_get_info(matrix=cfm_W_ikp_freq_j, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp)
      index_to_cell => kpoints%index_to_cell
      num_cells = SIZE(index_to_cell, 2)

      iatom_old = 0
      jatom_old = 0

      DO irow = 1, nrow_local
         DO jcol = 1, ncol_local

            i_bf = row_indices(irow)
            j_bf = col_indices(jcol)

            iatom = atom_from_bf_index(i_bf)
            jatom = atom_from_bf_index(j_bf)

            IF (PRESENT(wkp_ext)) THEN
               wkp_of_ikp = wkp_ext
            ELSE
               SELECT CASE (bs_env%l_RI(i_bf) + bs_env%l_RI(j_bf))
               CASE (0)
                  ! both RI functions are s-functions, k-extrapolation for 2D and 3D
                  wkp_of_ikp = wkp(ikp)
               CASE (1)
                  ! one function is an s-function, the other a p-function, k-extrapolation for 3D
                  wkp_of_ikp = bs_env%wkp_s_p(ikp)
               CASE DEFAULT
                  ! for any other matrix element of W, there is no need for extrapolation
                  wkp_of_ikp = bs_env%wkp_no_extra(ikp)
               END SELECT
            END IF

            IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN

               CALL compute_weight_re_im(weight_re, weight_im, &
                                         num_cells, iatom, jatom, xkp(1:3, ikp), wkp_of_ikp, &
                                         cell, index_to_cell, hmat, particle_set)

               iatom_old = iatom
               jatom_old = jatom

            END IF

            contribution = weight_re*REAL(cfm_W_ikp_freq_j%local_data(irow, jcol)) + &
                           weight_im*AIMAG(cfm_W_ikp_freq_j%local_data(irow, jcol))

            fm_W_MIC_freq_j%local_data(irow, jcol) = fm_W_MIC_freq_j%local_data(irow, jcol) &
                                                     + contribution

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE MIC_contribution_from_ikp

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
!> \param kpoints ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE kpoint_init_cell_index_simple(kpoints, qs_env)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_init_cell_index_simple'

      INTEGER                                            :: handle, nimages
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, para_env, sab_orb)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, dft_control=dft_control, sab_orb=sab_orb)
      nimages = dft_control%nimages
      CALL kpoint_init_cell_index(kpoints, sab_orb, para_env, dft_control)

      ! set back dft_control%nimages
      dft_control%nimages = nimages

      CALL timestop(handle)

   END SUBROUTINE kpoint_init_cell_index_simple

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE soc(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'soc'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! V^SOC_µν^(α),R = ħ/2 < ϕ_µ cell O | sum_ℓ ΔV_ℓ^SO(r,r') L^(α) | ϕ_ν cell R>, α = x,y,z
      ! see Hartwigsen, Goedecker, Hutter, Eq.(18), (19) (doi.org/10.1103/PhysRevB.58.3641)
      CALL V_SOC_xyz_from_pseudopotential(qs_env, bs_env%mat_V_SOC_xyz)

      ! Calculate H^SOC_µν,σσ'(k) = sum_α V^SOC_µν^(α)(k)*Pauli-matrix^(α)_σσ'
      ! see Hartwigsen, Goedecker, Hutter, Eq.(18) (doi.org/10.1103/PhysRevB.58.3641)
      SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)
      CASE (large_cell_Gamma)

         ! H^SOC_µν,σσ' = sum_α V^SOC_µν^(α)*Pauli-matrix^(α)_σσ'
         CALL H_KS_spinor_Gamma(bs_env)

      CASE (small_cell_full_kp)

         ! V^SOC_µν^(α),R -> V^SOC_µν^(α)(k); then calculate spinor H^SOC_µν,σσ'(k) (see above)
         CALL H_KS_spinor_kp(qs_env, bs_env)

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE soc

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE H_KS_spinor_Gamma(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'H_KS_spinor_Gamma'

      INTEGER                                            :: handle, nao, s
      TYPE(cp_fm_struct_type), POINTER                   :: str

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(bs_env%fm_ks_Gamma(1), nrow_global=nao)

      ALLOCATE (bs_env%cfm_SOC_spinor_ao(1))
      CALL create_cfm_double(bs_env%cfm_SOC_spinor_ao(1), fm_orig=bs_env%fm_ks_Gamma(1))
      CALL cp_cfm_set_all(bs_env%cfm_SOC_spinor_ao(1), z_zero)

      str => bs_env%fm_ks_Gamma(1)%matrix_struct

      s = nao + 1

      ! careful: inside add_dbcsr_submat, mat_V_SOC_xyz is multiplied by i because the real matrix
      !          mat_V_SOC_xyz is antisymmetric as V_SOC matrix is purely imaginary and Hermitian
      CALL add_dbcsr_submat(bs_env%cfm_SOC_spinor_ao(1), bs_env%mat_V_SOC_xyz(1, 1)%matrix, &
                            str, s, 1, z_one, .TRUE.)
      CALL add_dbcsr_submat(bs_env%cfm_SOC_spinor_ao(1), bs_env%mat_V_SOC_xyz(2, 1)%matrix, &
                            str, s, 1, gaussi, .TRUE.)
      CALL add_dbcsr_submat(bs_env%cfm_SOC_spinor_ao(1), bs_env%mat_V_SOC_xyz(3, 1)%matrix, &
                            str, 1, 1, z_one, .FALSE.)
      CALL add_dbcsr_submat(bs_env%cfm_SOC_spinor_ao(1), bs_env%mat_V_SOC_xyz(3, 1)%matrix, &
                            str, s, s, -z_one, .FALSE.)

      CALL timestop(handle)

   END SUBROUTINE H_KS_spinor_Gamma

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE H_KS_spinor_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'H_KS_spinor_kp'

      INTEGER                                            :: handle, i_dim, ikp, n_spin, &
                                                            nkp_bs_and_DOS, s
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_scf
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      TYPE(cp_cfm_type)                                  :: cfm_V_SOC_xyz_ikp
      TYPE(cp_fm_struct_type), POINTER                   :: str
      TYPE(kpoint_type), POINTER                         :: kpoints_scf
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      nkp_bs_and_DOS = bs_env%nkp_bs_and_DOS
      n_spin = bs_env%n_spin
      s = bs_env%n_ao + 1
      str => bs_env%cfm_ks_kp(1, 1)%matrix_struct

      CALL cp_cfm_create(cfm_V_SOC_xyz_ikp, bs_env%cfm_work_mo%matrix_struct)

      CALL alloc_cfm_double_array_1d(bs_env%cfm_SOC_spinor_ao, bs_env%cfm_ks_kp(1, 1), nkp_bs_and_DOS)

      CALL get_qs_env(qs_env, kpoints=kpoints_scf)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints_scf, sab_nl=sab_nl, cell_to_index=cell_to_index_scf)

      DO i_dim = 1, 3

         DO ikp = 1, nkp_bs_and_DOS

            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, ikp)

            CALL cp_cfm_set_all(cfm_V_SOC_xyz_ikp, z_zero)

            CALL rsmat_to_kp(bs_env%mat_V_SOC_xyz, i_dim, xkp, cell_to_index_scf, &
                             sab_nl, bs_env, cfm_V_SOC_xyz_ikp, imag_rs_mat=.TRUE.)

            ! multiply V_SOC with i because bs_env%mat_V_SOC_xyz stores imag. part (real part = 0)
            CALL cp_cfm_scale(gaussi, cfm_V_SOC_xyz_ikp)

            SELECT CASE (i_dim)
            CASE (1)
               ! add V^SOC_x * σ_x for σ_x = ( (0,1) (1,0) )
               CALL add_cfm_submat(bs_env%cfm_SOC_spinor_ao(ikp), cfm_V_SOC_xyz_ikp, 1, s)
               CALL add_cfm_submat(bs_env%cfm_SOC_spinor_ao(ikp), cfm_V_SOC_xyz_ikp, s, 1)
            CASE (2)
               ! add V^SOC_y * σ_y for σ_y = ( (0,-i) (i,0) )
               CALL cp_cfm_scale(gaussi, cfm_V_SOC_xyz_ikp)
               CALL add_cfm_submat(bs_env%cfm_SOC_spinor_ao(ikp), cfm_V_SOC_xyz_ikp, 1, s)
               CALL cp_cfm_scale(-z_one, cfm_V_SOC_xyz_ikp)
               CALL add_cfm_submat(bs_env%cfm_SOC_spinor_ao(ikp), cfm_V_SOC_xyz_ikp, s, 1)
            CASE (3)
               ! add V^SOC_z * σ_z for σ_z = ( (1,0) (0,1) )
               CALL add_cfm_submat(bs_env%cfm_SOC_spinor_ao(ikp), cfm_V_SOC_xyz_ikp, 1, 1)
               CALL cp_cfm_scale(-z_one, cfm_V_SOC_xyz_ikp)
               CALL add_cfm_submat(bs_env%cfm_SOC_spinor_ao(ikp), cfm_V_SOC_xyz_ikp, s, s)
            END SELECT

         END DO

      END DO ! ikp

      CALL cp_cfm_release(cfm_V_SOC_xyz_ikp)

      CALL timestop(handle)

   END SUBROUTINE H_KS_spinor_kp

! **************************************************************************************************
!> \brief ...
!> \param cfm_array ...
!> \param cfm_template ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE alloc_cfm_double_array_1d(cfm_array, cfm_template, n)
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_array
      TYPE(cp_cfm_type)                                  :: cfm_template
      INTEGER                                            :: n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_cfm_double_array_1d'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      ALLOCATE (cfm_array(n))
      DO i = 1, n
         CALL create_cfm_double(cfm_array(i), cfm_orig=cfm_template)
         CALL cp_cfm_set_all(cfm_array(i), z_zero)
      END DO

      CALL timestop(handle)

   END SUBROUTINE alloc_cfm_double_array_1d

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_all_VBM_CBM_bandgaps(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_all_VBM_CBM_bandgaps'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL get_VBM_CBM_bandgaps(bs_env%band_edges_scf, bs_env%eigenval_scf, bs_env)
      CALL get_VBM_CBM_bandgaps(bs_env%band_edges_G0W0, bs_env%eigenval_G0W0, bs_env)
      CALL get_VBM_CBM_bandgaps(bs_env%band_edges_HF, bs_env%eigenval_HF, bs_env)

      CALL timestop(handle)

   END SUBROUTINE get_all_VBM_CBM_bandgaps

! **************************************************************************************************
!> \brief ...
!> \param band_edges ...
!> \param ev ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_VBM_CBM_bandgaps(band_edges, ev, bs_env)
      TYPE(band_edges_type)                              :: band_edges
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: ev
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_VBM_CBM_bandgaps'

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

   END SUBROUTINE get_VBM_CBM_bandgaps

END MODULE post_scf_bandstructure_utils

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
MODULE post_scf_bandstructure_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_geeig_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_parser_methods,               ONLY: read_float_object
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_symmetric
   USE gw_utils,                        ONLY: compute_xkp,&
                                              kpoint_init_cell_index_simple
   USE input_constants,                 ONLY: int_ldos_z
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              max_line_length
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: complex_diag,&
                                              inv_3x3
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
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE rpa_gw_im_time_util,             ONLY: get_atom_index_from_basis_function_index
   USE scf_control_types,               ONLY: scf_control_type
   USE soc_pseudopotential_utils,       ONLY: add_cfm_submat,&
                                              get_cfm_submat
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: create_and_init_bs_env, &
             bandstructure_primitive_cell, bandstructure_primitive_cell_spinor, &
             dos_pdos_ldos, cfm_ikp_from_fm_Gamma

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

      CALL setup_kpoints_DOS(qs_env, bs_env, bs_env%kpoints_DOS)

      CALL allocate_and_fill_fm_ks_fm_s(qs_env, bs_env)

      CALL diagonalize_ks_matrix(bs_env)

      CALL check_positive_definite_overlap_mat(bs_env, qs_env)

      IF (bs_env%do_bs) THEN
         CALL setup_kpoints_bandstructure(bs_env, bs_env%kpoints_bandstructure)
         CALL setup_primitive_cell_for_bandstructure(qs_env, bs_env)
      END IF

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
      TYPE(section_vals_type), POINTER                   :: gw_sec, kp_bs_sec, lvcg_sec, soc_sec

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

      NULLIFY (lvcg_sec)
      lvcg_sec => section_vals_get_subs_vals(bs_sec, "LOCAL_VBM_CBM_GAP")
      CALL section_vals_get(lvcg_sec, explicit=bs_env%do_local_vbm_cbm_gap)

      CALL section_vals_val_get(lvcg_sec, "INTEGRATION", i_val=bs_env%int_ldos_xyz)
      CALL section_vals_val_get(lvcg_sec, "EPS_LDOS", r_vals=bs_env%eps_ldos)
      CALL section_vals_val_get(lvcg_sec, "ENERGY_WINDOW", r_vals=bs_env%energy_window_ldos)
      CALL section_vals_val_get(lvcg_sec, "ENERGY_STEP", r_val=bs_env%energy_step_ldos)
      CALL section_vals_val_get(lvcg_sec, "BIN_MESH", i_vals=bs_env%bin_mesh)

      NULLIFY (kp_bs_sec)
      kp_bs_sec => section_vals_get_subs_vals(bs_sec, "BANDSTRUCTURE_PATH")
      CALL section_vals_get(kp_bs_sec, explicit=bs_env%do_bs)
      CALL section_vals_val_get(kp_bs_sec, "NPOINTS", i_val=bs_env%input_kp_bs_npoints)
      CALL section_vals_val_get(kp_bs_sec, "SPECIAL_POINT", n_rep_val=bs_env%input_kp_bs_n_sp_pts)

      ! read special points for band structure
      ALLOCATE (bs_env%xkp_special(3, bs_env%input_kp_bs_n_sp_pts))
      ALLOCATE (bs_env%kp_special_name(bs_env%input_kp_bs_n_sp_pts))
      DO ikp = 1, bs_env%input_kp_bs_n_sp_pts
         CALL section_vals_val_get(kp_bs_sec, "SPECIAL_POINT", i_rep_val=ikp, c_vals=string_ptr)
         CPASSERT(SIZE(string_ptr(:), 1) == 4)
         bs_env%kp_special_name(ikp) = string_ptr(1)
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
   SUBROUTINE setup_kpoints_DOS(qs_env, bs_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'setup_kpoints_DOS'

      INTEGER                                            :: handle, i_dim, nkp, nkp_of_kp_grid, u
      INTEGER, DIMENSION(3)                              :: nkp_grid, periodic

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      kpoints%kp_scheme = "GENERAL"

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

      nkp_of_kp_grid = nkp_grid(1)*nkp_grid(2)*nkp_grid(3)

      IF (periodic(1) == 0 .AND. periodic(2) == 0 .AND. periodic(3) == 0) THEN
         nkp = 1
      ELSE
         ! we include the Γ-point for periodic systems
         nkp = nkp_of_kp_grid + 1
      END IF

      kpoints%nkp_grid(1:3) = nkp_grid(1:3)
      kpoints%nkp = nkp

      bs_env%nkp_DOS = nkp

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))
      kpoints%wkp(1:nkp_of_kp_grid) = 1.0_dp/REAL(nkp_of_kp_grid, KIND=dp)
      CALL compute_xkp(kpoints%xkp, 1, nkp_of_kp_grid, nkp_grid)
      ! Γ-point
      kpoints%xkp(1:3, nkp) = 0.0_dp
      kpoints%wkp(nkp) = 0.0_dp

      CALL kpoint_init_cell_index_simple(kpoints, qs_env)

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,1A,T69,3I4)") "K-point mesh for the density of states (DOS)", &
            nkp_grid(1:3)
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_DOS

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE setup_kpoints_bandstructure(bs_env, kpoints)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_kpoints_bandstructure'

      INTEGER                                            :: handle, i_kp_in_line, i_special_kp, ikk, &
                                                            n_kp_in_line, n_special_kp, nkp, u

      CALL timeset(routineN, handle)

      ! routine adapted from mp2_integrals.F
      NULLIFY (kpoints)
      CALL kpoint_create(kpoints)

      n_special_kp = bs_env%input_kp_bs_n_sp_pts
      n_kp_in_line = bs_env%input_kp_bs_npoints

      nkp = n_kp_in_line*(n_special_kp - 1) + 1

      IF (n_special_kp < 1) &
         CPABORT("Please specify special k-points in the Brillouin zone via SPECIAL_POINT.")
      IF (n_kp_in_line < 1) &
         CPABORT("Please specify the number of k-points between special k-points.")

      ALLOCATE (kpoints%xkp(3, nkp))

      kpoints%nkp = nkp
      kpoints%xkp(1:3, 1) = bs_env%xkp_special(1:3, 1)

      bs_env%nkp_bs = nkp

      ikk = 1

      DO i_special_kp = 2, n_special_kp
         DO i_kp_in_line = 1, n_kp_in_line
            ikk = ikk + 1
            kpoints%xkp(1:3, ikk) = bs_env%xkp_special(1:3, i_special_kp - 1) + &
                                    REAL(i_kp_in_line, KIND=dp)/REAL(n_kp_in_line, KIND=dp)* &
                                    (bs_env%xkp_special(1:3, i_special_kp) - &
                                     bs_env%xkp_special(1:3, i_special_kp - 1))
         END DO
      END DO

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,1A,T77,I4)") "Number of special k-points for the bandstructure", &
            n_special_kp
         WRITE (UNIT=u, FMT="(T2,1A,T77,I4)") "Number of k-points for the bandstructure", nkp
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_kpoints_bandstructure

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

      ALLOCATE (bs_env%eigenval_scf_Gamma(bs_env%n_ao))

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
               WRITE (UNIT=u, FMT="(T2,A)") ""
               WRITE (UNIT=u, FMT="(T2,A)") "ERROR: The Cholesky decomposition "// &
                  "of the k-point overlap matrix failed. This is"
               WRITE (UNIT=u, FMT="(T2,A)") "because the algorithm is "// &
                  "only correct in the limit of large cells. The cell of "
               WRITE (UNIT=u, FMT="(T2,A)") "the calculation is too small. "// &
                  "Use MULTIPLE_UNIT_CELL to create a larger cell "
               WRITE (UNIT=u, FMT="(T2,A)") "and to prevent this error."
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

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      scf_control=scf_control, &
                      mos=mos)

      bs_env%n_spin = dft_control%nspins
      IF (bs_env%n_spin == 1) bs_env%spin_degeneracy = 2.0_dp
      IF (bs_env%n_spin == 2) bs_env%spin_degeneracy = 1.0_dp

      ! so far, only closed-shell
      CPASSERT(bs_env%n_spin == 1)

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

      u = bs_env%unit_nr

      IF (u > 0) THEN
         WRITE (UNIT=u, FMT="(T2,2A,T73,I8)") "Number of occupied molecular orbitals (MOs) ", &
            "= Number of occupied bands", homo
         WRITE (UNIT=u, FMT="(T2,2A,T73,I8)") "Number of unoccupied (= virtual) MOs ", &
            "= Number of unoccupied bands", n_ao - homo
         WRITE (UNIT=u, FMT="(T2,A,T73,I8)") "Number of Gaussian basis functions for MOs", n_ao
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

      bs_env%n_bins_max_for_printing = 1000

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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=bs_env%n_ao, &
                               ncol_global=bs_env%n_ao, para_env=para_env)

      DO i_work = 1, SIZE(bs_env%fm_work_mo)
         CALL cp_fm_create(bs_env%fm_work_mo(i_work), fm_struct)
      END DO

      CALL cp_fm_create(bs_env%fm_s_Gamma, fm_struct)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, bs_env%fm_s_Gamma)

      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(bs_env%fm_ks_Gamma(ispin), fm_struct)
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, bs_env%fm_ks_Gamma(ispin))
         CALL cp_fm_create(bs_env%fm_mo_coeff_Gamma(ispin), fm_struct)
      END DO

      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (bs_env%mat_ao_ao%matrix)
      ALLOCATE (bs_env%mat_ao_ao%matrix)
      CALL dbcsr_create(bs_env%mat_ao_ao%matrix, template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ALLOCATE (bs_env%eigenval_scf(bs_env%n_ao, bs_env%nkp_DOS))

      CALL timestop(handle)

   END SUBROUTINE allocate_and_fill_fm_ks_fm_s

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE setup_primitive_cell_for_bandstructure(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_primitive_cell_for_bandstructure'

      INTEGER :: handle, i_atom, i_x_cell, index_j, j_atom, j_atom_prim_cell, j_y_cell, k_z_cell, &
         n_atom, n_atom_in_primitive_cell, n_max_check, n_max_x, n_max_y, n_max_z, &
         n_mult_unit_cell_x, n_mult_unit_cell_y, n_mult_unit_cell_z, n_primitive_cells, n_x, n_y, &
         n_z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      LOGICAL :: i_atom_has_image_in_every_subcell, i_atom_has_image_in_subcell_ijk
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: valid_multiple_unit_cell
      REAL(dp), DIMENSION(3) :: center_primitive_cell, coord_ijk, coord_sub_cell_ijk, &
         index_atom_i, index_ijk, index_sub_cell_ijk, offset, ra, ra_ijk, rab, rb
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set, cell=cell)
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)
      CALL get_cell(cell=cell, h=hmat)

      ! automatically check structure for smallest possible unit cell
      n_max_check = 20
      n_max_x = n_max_check*bs_env%periodic(1) + 1
      n_max_y = n_max_check*bs_env%periodic(2) + 1
      n_max_z = n_max_check*bs_env%periodic(3) + 1

      ALLOCATE (valid_multiple_unit_cell(n_max_x, n_max_y, n_max_z))
      valid_multiple_unit_cell(:, :, :) = .TRUE.

      n_atom = bs_env%n_atom

      DO i_atom = 1, n_atom

         IF (.NOT. MODULO(i_atom, bs_env%para_env%num_pe) == bs_env%para_env%mepos) CYCLE

         ra(1:3) = particle_set(i_atom)%r(1:3)

         DO n_x = 1, n_max_x
            DO n_y = 1, n_max_y
               DO n_z = 1, n_max_z

                  i_atom_has_image_in_every_subcell = .TRUE.

                  DO i_x_cell = 0, n_x - 1
                     DO j_y_cell = 0, n_y - 1
                        DO k_z_cell = 0, n_z - 1

                           i_atom_has_image_in_subcell_ijk = .FALSE.

                           DO j_atom = 1, n_atom

                              IF (kind_of(i_atom) .NE. kind_of(j_atom)) CYCLE

                              IF (i_atom_has_image_in_subcell_ijk) CYCLE

                              IF (.NOT. i_atom_has_image_in_every_subcell) CYCLE

                              index_sub_cell_ijk(1:3) = (/REAL(i_x_cell, dp)/REAL(n_x, dp), &
                                                          REAL(j_y_cell, dp)/REAL(n_y, dp), &
                                                          REAL(k_z_cell, dp)/REAL(n_z, dp)/)

                              coord_sub_cell_ijk(1:3) = MATMUL(hmat, index_sub_cell_ijk)

                              ra_ijk(1:3) = ra(1:3) + coord_sub_cell_ijk(1:3)

                              rb(1:3) = pbc(particle_set(j_atom)%r(1:3), cell)

                              rab(1:3) = rb(1:3) - pbc(ra_ijk(1:3), cell)

                              dab = SQRT((rab(1))**2 + (rab(2))**2 + (rab(3))**2)

                              IF (dab < 1.0E-5) i_atom_has_image_in_subcell_ijk = .TRUE.

                           END DO

                           IF (.NOT. i_atom_has_image_in_subcell_ijk) THEN
                              i_atom_has_image_in_every_subcell = .FALSE.
                           END IF

                        END DO
                     END DO
                  END DO

                  ! a valid multiple unit cell must be valid for all atoms
                  valid_multiple_unit_cell(n_x, n_y, n_z) = i_atom_has_image_in_every_subcell .AND. &
                                                            valid_multiple_unit_cell(n_x, n_y, n_z)

               END DO
            END DO
         END DO

      END DO

      CALL mpi_AND(valid_multiple_unit_cell, bs_env%para_env)

      n_mult_unit_cell_x = 1
      n_mult_unit_cell_y = 1
      n_mult_unit_cell_z = 1

      DO n_x = 1, n_max_x
         DO n_y = 1, n_max_y
            DO n_z = 1, n_max_z
               IF (valid_multiple_unit_cell(n_x, n_y, n_z)) THEN
                  n_mult_unit_cell_x = MAX(n_mult_unit_cell_x, n_x)
                  n_mult_unit_cell_y = MAX(n_mult_unit_cell_y, n_y)
                  n_mult_unit_cell_z = MAX(n_mult_unit_cell_z, n_z)
               END IF
            END DO
         END DO
      END DO

      bs_env%multiple_unit_cell(1) = n_mult_unit_cell_x
      bs_env%multiple_unit_cell(2) = n_mult_unit_cell_y
      bs_env%multiple_unit_cell(3) = n_mult_unit_cell_z

      IF (n_mult_unit_cell_x .NE. 1 .OR. &
          n_mult_unit_cell_y .NE. 1 .OR. &
          n_mult_unit_cell_z .NE. 1) THEN
         bs_env%calculate_bandstructure_of_primitive_cell = .TRUE.
      ELSE
         bs_env%calculate_bandstructure_of_primitive_cell = .FALSE.
      END IF

      n_atom_in_primitive_cell = n_atom/n_mult_unit_cell_x/n_mult_unit_cell_y/n_mult_unit_cell_z
      bs_env%n_atom_in_primitive_cell = n_atom_in_primitive_cell

      n_primitive_cells = n_atom/n_atom_in_primitive_cell
      bs_env%n_primitive_cells = n_primitive_cells

      bs_env%hmat_primitive_cell(1, 1:3) = hmat(1, 1:3)/REAL(n_mult_unit_cell_x)
      bs_env%hmat_primitive_cell(2, 1:3) = hmat(2, 1:3)/REAL(n_mult_unit_cell_y)
      bs_env%hmat_primitive_cell(3, 1:3) = hmat(3, 1:3)/REAL(n_mult_unit_cell_z)

      bs_env%hinv_primitive_cell = inv_3x3(bs_env%hmat_primitive_cell)

      bs_env%do_bs_primitive_cell = bs_env%do_bs .AND. n_atom_in_primitive_cell < 20 &
                                    .AND. n_primitive_cells > 1

      ALLOCATE (bs_env%atoms_i_primitive_cell(n_atom_in_primitive_cell))
      bs_env%atoms_i_primitive_cell(:) = 0

      ! just a small offset to avoid that atoms are precisely on egdes or faces
      offset(1:3) = MATMUL(bs_env%hmat_primitive_cell, (/0.001_dp, 0.001_dp, 0.001_dp/))
      center_primitive_cell(1:3) = pbc(particle_set(1)%r(1:3), cell) - offset(1:3)

      index_j = 0
      DO i_atom = 1, n_atom

         rb(1:3) = pbc(particle_set(i_atom)%r(1:3), cell) - center_primitive_cell(1:3)

         index_atom_i(1:3) = MATMUL(bs_env%hinv_primitive_cell, rb)

         IF (index_atom_i(1) > -0.5_dp .AND. index_atom_i(1) < 0.5_dp .AND. &
             index_atom_i(2) > -0.5_dp .AND. index_atom_i(2) < 0.5_dp .AND. &
             index_atom_i(3) > -0.5_dp .AND. index_atom_i(3) < 0.5_dp) THEN

            index_j = index_j + 1
            CPASSERT(index_j .LE. n_atom_in_primitive_cell)
            bs_env%atoms_i_primitive_cell(index_j) = i_atom

         END IF

      END DO

      ALLOCATE (bs_env%ref_atom_primitive_cell(n_atom))
      ALLOCATE (bs_env%cell_of_i_atom(n_atom, 3))

      DO i_atom = 1, n_atom

         ra(1:3) = pbc(particle_set(i_atom)%r(1:3), cell)

         DO j_atom_prim_cell = 1, n_atom_in_primitive_cell

            j_atom = bs_env%atoms_i_primitive_cell(j_atom_prim_cell)

            rb(1:3) = pbc(particle_set(j_atom)%r(1:3), cell)

            DO i_x_cell = -n_mult_unit_cell_x/2, n_mult_unit_cell_x/2
               DO j_y_cell = -n_mult_unit_cell_y/2, n_mult_unit_cell_y/2
                  DO k_z_cell = -n_mult_unit_cell_z/2, n_mult_unit_cell_z/2

                     index_ijk(1:3) = (/REAL(i_x_cell, dp), REAL(j_y_cell, dp), REAL(k_z_cell, dp)/)
                     coord_ijk(1:3) = MATMUL(bs_env%hmat_primitive_cell, index_ijk)

                     ra_ijk(1:3) = ra(1:3) + coord_ijk(1:3)

                     rab(1:3) = rb(1:3) - pbc(ra_ijk(1:3), cell)

                     dab = SQRT((rab(1))**2 + (rab(2))**2 + (rab(3))**2)

                     IF (dab < 1.0E-5) THEN
                        bs_env%ref_atom_primitive_cell(i_atom) = j_atom
                        bs_env%cell_of_i_atom(i_atom, 1) = i_x_cell
                        bs_env%cell_of_i_atom(i_atom, 2) = j_y_cell
                        bs_env%cell_of_i_atom(i_atom, 3) = k_z_cell
                     END IF

                  END DO
               END DO
            END DO
         END DO
      END DO
      IF (bs_env%unit_nr > 0 .AND. bs_env%calculate_bandstructure_of_primitive_cell) THEN
         WRITE (bs_env%unit_nr, '(T2,A,3I4)') &
            'Detected a multiple unit cell (will be used for band structure)    ', &
            bs_env%multiple_unit_cell
         WRITE (bs_env%unit_nr, '(T2,A,I28)') &
            'Number of occupied bands in the primitive unit cell', &
            bs_env%n_occ(1)/bs_env%n_primitive_cells
         WRITE (bs_env%unit_nr, '(T2,A,I26)') &
            'Number of unoccupied bands in the primitive unit cell', &
            bs_env%n_vir(1)/bs_env%n_primitive_cells
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_primitive_cell_for_bandstructure

! **************************************************************************************************
!> \brief ...
!> \param logical_array_3d ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE mpi_AND(logical_array_3d, para_env)
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: logical_array_3d
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mpi_AND'

      INTEGER                                            :: handle, i, j, k, n_1, n_2, n_3
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: integer_array_3d

      CALL timeset(routineN, handle)

      n_1 = SIZE(logical_array_3d, 1)
      n_2 = SIZE(logical_array_3d, 2)
      n_3 = SIZE(logical_array_3d, 3)

      ALLOCATE (integer_array_3d(n_1, n_2, n_3))
      integer_array_3d(:, :, :) = 0

      DO i = 1, n_1
         DO j = 1, n_2
            DO k = 1, n_3
               IF (logical_array_3d(i, j, k)) integer_array_3d(i, j, k) = 1
            END DO
         END DO
      END DO

      CALL para_env%sync()
      CALL para_env%sum(integer_array_3d)
      CALL para_env%sync()

      logical_array_3d(:, :, :) = .FALSE.

      DO i = 1, n_1
         DO j = 1, n_2
            DO k = 1, n_3
               IF (integer_array_3d(i, j, k) == para_env%num_pe) logical_array_3d(i, j, k) = .TRUE.
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE mpi_AND

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param eigenvalues ...
!> \param filename ...
!> \param fm_h_Gamma ...
! **************************************************************************************************
   SUBROUTINE bandstructure_primitive_cell(qs_env, bs_env, eigenvalues, filename, fm_h_Gamma)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvalues
      CHARACTER(LEN=*)                                   :: filename
      TYPE(cp_fm_type)                                   :: fm_h_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bandstructure_primitive_cell'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: h_munu_k, mo_coeff_k, s_munu_k
      INTEGER :: col_global, handle, i, i_atom, i_dim, i_row, ikp, imo, ip, iunit, j, j_atom, &
         j_col, n_ao, n_ao_primitive_cell, n_atom, ncol_local, nrow_local, ref_atom_j, row_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf, first_bf_from_atom, &
                                                            first_bf_of_primit_atom
      INTEGER, DIMENSION(3)                              :: cell_atom_i, cell_atom_j, min_max_cell
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      n_ao_primitive_cell = n_ao/bs_env%n_primitive_cells
      n_atom = bs_env%n_atom

      ALLOCATE (h_munu_k(n_ao_primitive_cell, n_ao_primitive_cell))
      ALLOCATE (s_munu_k(n_ao_primitive_cell, n_ao_primitive_cell))
      ALLOCATE (mo_coeff_k(n_ao_primitive_cell, n_ao_primitive_cell))
      ALLOCATE (eigenvalues(n_ao_primitive_cell, bs_env%nkp_bs))

      ALLOCATE (atom_from_bf(n_ao))
      ALLOCATE (first_bf_from_atom(n_atom))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf, n_ao, "ORB", &
                                                    first_bf_from_atom)

      ALLOCATE (first_bf_of_primit_atom(n_atom))
      CALL get_basis_function_index_of_primitive_atoms(bs_env, first_bf_of_primit_atom, &
                                                       first_bf_from_atom)

      IF (bs_env%para_env%is_source()) THEN
         CALL open_file(filename, unit_number=iunit, file_status="REPLACE", &
                        file_action="WRITE", file_position="APPEND")
      ELSE
         iunit = -1
      END IF

      IF (iunit > 0) THEN

         WRITE (UNIT=iunit, FMT="(2(A,I0),A)") "# ", &
            bs_env%input_kp_bs_n_sp_pts, " special points, ", bs_env%nkp_bs, " k-points"
         DO ip = 1, bs_env%input_kp_bs_n_sp_pts
            WRITE (UNIT=iunit, FMT="(A,I0,T20,T24,3(1X,F14.8),2X,A)") &
               "#  Special point ", ip, bs_env%xkp_special(1:3, ip), &
               ADJUSTL(TRIM(bs_env%kp_special_name(ip)))
         END DO

      END IF

      CALL cp_fm_get_info(matrix=fm_h_Gamma, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO i_dim = 1, 3
         min_max_cell(i_dim) = MIN(MAXVAL(bs_env%cell_of_i_atom(:, i_dim)), &
                                   MAXVAL(-bs_env%cell_of_i_atom(:, i_dim)))
      END DO

      DO ikp = 1, bs_env%nkp_bs

         h_munu_k = z_zero
         s_munu_k = z_zero

         DO i_row = 1, nrow_local
            DO j_col = 1, ncol_local

               row_global = row_indices(i_row)
               col_global = col_indices(j_col)

               i_atom = atom_from_bf(row_global)
               j_atom = atom_from_bf(col_global)

               cell_atom_i = bs_env%cell_of_i_atom(i_atom, 1:3)

               ! atom_i must be in the primitive cell (0,0,0)
               ! (because we calculate h_mu,nu(k) = \sum_R <mu,cell o|h|nu,cell R>
               IF (ANY(cell_atom_i(1:3) .NE. 0)) CYCLE

               cell_atom_j = bs_env%cell_of_i_atom(j_atom, 1:3)

               ! only consider symmetric cell summation, i.e. cell (4,-2,0) needs to have
               ! counterpart (-4,2,0). In case we have 7x7 cell, (-4,2,0) will be absent
               IF (ANY(ABS(cell_atom_j(1:3)) > min_max_cell(1:3))) CYCLE

               arg = (REAL(cell_atom_j(1), dp)*bs_env%kpoints_bandstructure%xkp(1, ikp) + &
                      REAL(cell_atom_j(2), dp)*bs_env%kpoints_bandstructure%xkp(2, ikp) + &
                      REAL(cell_atom_j(3), dp)*bs_env%kpoints_bandstructure%xkp(3, ikp))*twopi

               ref_atom_j = bs_env%ref_atom_primitive_cell(j_atom)

               i = row_global - first_bf_from_atom(i_atom) + first_bf_of_primit_atom(i_atom)
               j = col_global - first_bf_from_atom(j_atom) + first_bf_of_primit_atom(ref_atom_j)

               h_munu_k(i, j) = h_munu_k(i, j) + &
                                COS(arg)*fm_h_Gamma%local_data(i_row, j_col)*z_one + &
                                SIN(arg)*fm_h_Gamma%local_data(i_row, j_col)*gaussi

               s_munu_k(i, j) = s_munu_k(i, j) + &
                                COS(arg)*bs_env%fm_s_Gamma%local_data(i_row, j_col)*z_one + &
                                SIN(arg)*bs_env%fm_s_Gamma%local_data(i_row, j_col)*gaussi
            END DO ! j_col
         END DO ! i_row

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(h_munu_k)
         CALL bs_env%para_env%sum(s_munu_k)
         CALL bs_env%para_env%sync()

         CALL complex_geeig(h_munu_k, s_munu_k, mo_coeff_k, eigenvalues(:, ikp))

         IF (iunit > 0) THEN

            WRITE (UNIT=iunit, FMT="(A,I0,T15,A,T24,3(1X,F14.8))") &
               "#  Point ", ikp, ":", bs_env%kpoints_bandstructure%xkp(1:3, ikp)
            WRITE (UNIT=iunit, FMT="(A)") "#   Band    Energy [eV]"
            DO imo = 1, n_ao_primitive_cell
               WRITE (UNIT=iunit, FMT="(T2,I7,1X,F14.8)") imo, eigenvalues(imo, ikp)*evolt
            END DO

         END IF

      END DO ! ikp

      IF (bs_env%para_env%is_source()) CALL close_file(unit_number=iunit)

      CALL timestop(handle)

   END SUBROUTINE bandstructure_primitive_cell

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param eigenvalues ...
!> \param filename ...
!> \param cfm_h_Gamma_spinor ...
! **************************************************************************************************
   SUBROUTINE bandstructure_primitive_cell_spinor(qs_env, bs_env, eigenvalues, filename, &
                                                  cfm_h_Gamma_spinor)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvalues
      CHARACTER(LEN=*)                                   :: filename
      TYPE(cp_cfm_type)                                  :: cfm_h_Gamma_spinor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bandstructure_primitive_cell_spinor'

      COMPLEX(KIND=dp)                                   :: arg, s_z
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: s_dot_mo_coeff_down, s_dot_mo_coeff_up
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: h_munu_k, mo_coeff_k, s_munu_k, &
                                                            s_munu_k_single
      INTEGER :: col_global, handle, i, i_atom, i_atom_non_spinor, i_dim, i_row, ikp, imo, ip, &
         iunit, j, j_atom, j_atom_non_spinor, j_col, n_ao, n_ao_primitive_cell, n_atom, &
         n_atom_primitive_cell, ncol_local, nrow_local, ref_atom_j, row_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf, first_bf_from_atom, &
                                                            first_bf_of_primit_atom
      INTEGER, DIMENSION(3)                              :: cell_atom_i, cell_atom_j, min_max_cell
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      n_ao_primitive_cell = n_ao/bs_env%n_primitive_cells
      n_atom = bs_env%n_atom
      n_atom_primitive_cell = n_atom/bs_env%n_primitive_cells

      ALLOCATE (h_munu_k(2*n_ao_primitive_cell, 2*n_ao_primitive_cell))
      ALLOCATE (s_munu_k(2*n_ao_primitive_cell, 2*n_ao_primitive_cell))
      ALLOCATE (s_munu_k_single(n_ao_primitive_cell, n_ao_primitive_cell))
      ALLOCATE (mo_coeff_k(2*n_ao_primitive_cell, 2*n_ao_primitive_cell))
      ALLOCATE (eigenvalues(2*n_ao_primitive_cell, bs_env%nkp_bs))
      ALLOCATE (s_dot_mo_coeff_up(n_ao_primitive_cell))
      ALLOCATE (s_dot_mo_coeff_down(n_ao_primitive_cell))

      ALLOCATE (atom_from_bf(2*n_ao))
      ALLOCATE (first_bf_from_atom(2*n_atom))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf, n_ao, "ORB", &
                                                    first_bf_from_atom)
      atom_from_bf(n_ao + 1:2*n_ao) = atom_from_bf(1:n_ao) + n_atom
      first_bf_from_atom(n_atom + 1:2*n_atom) = first_bf_from_atom(1:n_atom) + n_ao

      ALLOCATE (first_bf_of_primit_atom(2*n_atom))
      CALL get_basis_function_index_of_primitive_atoms(bs_env, first_bf_of_primit_atom, &
                                                       first_bf_from_atom)
      first_bf_of_primit_atom(n_atom + 1:2*n_atom) = first_bf_of_primit_atom(1:n_atom) &
                                                     + n_ao_primitive_cell

      IF (bs_env%para_env%is_source()) THEN
         CALL open_file(filename, unit_number=iunit, file_status="REPLACE", &
                        file_action="WRITE", file_position="APPEND")
      ELSE
         iunit = -1
      END IF

      IF (iunit > 0) THEN

         WRITE (UNIT=iunit, FMT="(2(A,I0),A)") "# ", &
            bs_env%input_kp_bs_n_sp_pts, " special points, ", bs_env%nkp_bs, " k-points"
         DO ip = 1, bs_env%input_kp_bs_n_sp_pts
            WRITE (UNIT=iunit, FMT="(A,I0,T20,T24,3(1X,F14.8),2X,A)") &
               "#  Special point ", ip, bs_env%xkp_special(1:3, ip), &
               ADJUSTL(TRIM(bs_env%kp_special_name(ip)))
         END DO

      END IF

      CALL cp_cfm_get_info(matrix=cfm_h_Gamma_spinor, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO i_dim = 1, 3
         min_max_cell(i_dim) = MIN(MAXVAL(bs_env%cell_of_i_atom(:, i_dim)), &
                                   MAXVAL(-bs_env%cell_of_i_atom(:, i_dim)))
      END DO

      DO ikp = 1, bs_env%nkp_bs

         h_munu_k = z_zero
         s_munu_k = z_zero

         DO i_row = 1, nrow_local
            DO j_col = 1, ncol_local

               row_global = row_indices(i_row)
               col_global = col_indices(j_col)

               i_atom = atom_from_bf(row_global)
               j_atom = atom_from_bf(col_global)

               IF (i_atom > n_atom) THEN
                  i_atom_non_spinor = i_atom - n_atom
               ELSE
                  i_atom_non_spinor = i_atom
               END IF

               IF (j_atom > n_atom) THEN
                  j_atom_non_spinor = j_atom - n_atom
               ELSE
                  j_atom_non_spinor = j_atom
               END IF

               cell_atom_i = bs_env%cell_of_i_atom(i_atom_non_spinor, 1:3)

               ! atom_i must be in the primitive cell (0,0,0)
               ! (because we calculate h_mu,nu(k) = \sum_R <mu,cell o|h|nu,cell R>
               IF (ANY(cell_atom_i(1:3) .NE. 0)) CYCLE

               cell_atom_j = bs_env%cell_of_i_atom(j_atom_non_spinor, 1:3)

               ! only consider symmetric cell summation, i.e. cell (4,-2,0) needs to have
               ! counterpart (-4,2,0). In case we have 7x7 cell, (-4,2,0) will be absent
               IF (ANY(ABS(cell_atom_j(1:3)) > min_max_cell(1:3))) CYCLE

               arg = (REAL(cell_atom_j(1), dp)*bs_env%kpoints_bandstructure%xkp(1, ikp) + &
                      REAL(cell_atom_j(2), dp)*bs_env%kpoints_bandstructure%xkp(2, ikp) + &
                      REAL(cell_atom_j(3), dp)*bs_env%kpoints_bandstructure%xkp(3, ikp))*twopi

               IF (j_atom > n_atom) THEN
                  ref_atom_j = bs_env%ref_atom_primitive_cell(j_atom_non_spinor) + n_atom
               ELSE
                  ref_atom_j = bs_env%ref_atom_primitive_cell(j_atom)
               END IF

               i = row_global - first_bf_from_atom(i_atom) + first_bf_of_primit_atom(i_atom)
               j = col_global - first_bf_from_atom(j_atom) + first_bf_of_primit_atom(ref_atom_j)

               h_munu_k(i, j) = h_munu_k(i, j) + &
                                COS(arg)*cfm_h_Gamma_spinor%local_data(i_row, j_col) + &
                                SIN(arg)*cfm_h_Gamma_spinor%local_data(i_row, j_col)*gaussi

               s_munu_k(i, j) = s_munu_k(i, j) + &
                                COS(arg)*bs_env%cfm_s_spinor_Gamma%local_data(i_row, j_col) + &
                                SIN(arg)*bs_env%cfm_s_spinor_Gamma%local_data(i_row, j_col)*gaussi

            END DO ! j_col
         END DO ! i_row

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(h_munu_k)
         CALL bs_env%para_env%sum(s_munu_k)
         CALL bs_env%para_env%sync()

         CALL complex_geeig(h_munu_k, s_munu_k, mo_coeff_k, eigenvalues(:, ikp))

         IF (iunit > 0) THEN

            s_munu_k_single(:, :) = s_munu_k(1:n_ao_primitive_cell, 1:n_ao_primitive_cell)

            WRITE (UNIT=iunit, FMT="(A,I0,T15,A,T24,3(1X,F14.8))") &
               "#  Point ", ikp, ":", bs_env%kpoints_bandstructure%xkp(1:3, ikp)
            WRITE (UNIT=iunit, FMT="(A)") "#   Band    Energy [eV]     <S_z> / (ħ/2) "
            DO imo = 1, 2*n_ao_primitive_cell
               s_dot_mo_coeff_up(:) = MATMUL(s_munu_k_single, &
                                             mo_coeff_k(1:n_ao_primitive_cell, imo))
               s_dot_mo_coeff_down(:) = MATMUL(s_munu_k_single, &
                                               mo_coeff_k(n_ao_primitive_cell + 1:, imo))
               s_z = SUM(CONJG(mo_coeff_k(1:n_ao_primitive_cell, imo))*s_dot_mo_coeff_up) - &
                     SUM(CONJG(mo_coeff_k(n_ao_primitive_cell + 1:, imo))*s_dot_mo_coeff_down)
               WRITE (UNIT=iunit, FMT="(T2,I7,1X,2F14.8)") imo, eigenvalues(imo, ikp)*evolt, &
                  REAL(s_z, KIND=dp)
            END DO

         END IF

      END DO ! ikp

      IF (bs_env%para_env%is_source()) CALL close_file(unit_number=iunit)

      CALL timestop(handle)

   END SUBROUTINE bandstructure_primitive_cell_spinor

! **************************************************************************************************
!> \brief Solves generalized, complex eigenvalue problem, HC = SCε by diagonalizing S^-0.5*H*S^-0.5
!> \param matrix ...
!> \param overlap ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE complex_geeig(matrix, overlap, eigenvectors, eigenvalues)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: matrix, overlap, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: overlap_sqrt_inv, work_1, work_2
      INTEGER                                            :: i, n
      LOGICAL                                            :: check_size

      n = SIZE(matrix, 1)

      check_size = SIZE(matrix, 2) == n .AND. SIZE(overlap, 1) == n .AND. &
                   SIZE(eigenvalues) == n .AND. &
                   SIZE(eigenvectors, 1) == n .AND. SIZE(eigenvectors, 2) == n
      CPASSERT(check_size)

      ALLOCATE (work_1(n, n), work_2(n, n), overlap_sqrt_inv(n, n))

      CALL complex_diag(overlap, eigenvectors, eigenvalues)

      work_1(:, :) = z_zero
      DO i = 1, n
         IF (eigenvalues(i) > 1.0E-5_dp) THEN
            work_1(i, i) = eigenvalues(i)**(-0.5_dp)
         END IF
      END DO
      work_2(:, :) = MATMUL(work_1, TRANSPOSE(CONJG(eigenvectors)))
      overlap_sqrt_inv(:, :) = MATMUL(eigenvectors, work_2)

      work_1(:, :) = MATMUL(matrix, overlap_sqrt_inv)
      work_2(:, :) = MATMUL(overlap_sqrt_inv, work_1)

      CALL complex_diag(work_2, eigenvectors, eigenvalues)

      work_2(:, :) = MATMUL(overlap_sqrt_inv, eigenvectors)
      eigenvectors(:, :) = work_2(:, :)

   END SUBROUTINE complex_geeig

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param first_bf_of_primit_atom ...
!> \param first_bf_from_atom ...
! **************************************************************************************************
   SUBROUTINE get_basis_function_index_of_primitive_atoms(bs_env, first_bf_of_primit_atom, &
                                                          first_bf_from_atom)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_bf_of_primit_atom, &
                                                            first_bf_from_atom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_basis_function_index_of_primitive_atoms'

      INTEGER                                            :: handle, i_atom, n_atom, n_bf_of_atom_i

      CALL timeset(routineN, handle)

      first_bf_of_primit_atom(:) = 1

      n_atom = bs_env%n_atom

      DO i_atom = 2, n_atom
         IF (ANY(bs_env%atoms_i_primitive_cell(:) == i_atom)) THEN
            n_bf_of_atom_i = first_bf_from_atom(i_atom) - first_bf_from_atom(i_atom - 1)
            first_bf_of_primit_atom(i_atom:n_atom) = first_bf_of_primit_atom(i_atom:n_atom) &
                                                     + n_bf_of_atom_i
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_basis_function_index_of_primitive_atoms

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE dos_pdos_ldos(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dos_pdos_ldos'

      INTEGER                                            :: handle, homo, homo_spinor, i_E, i_kind, &
                                                            i_mo, ikp, n_ao, n_E, &
                                                            n_primitive_cells, nkind
      REAL(KIND=dp) :: broadening, E_max, E_max_G0W0, E_min, E_min_G0W0, E_total_window, energy, &
         energy_step_DOS, energy_window_DOS, ev, wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: DOS_G0W0, DOS_G0W0_SOC, DOS_scf, &
                                                            DOS_scf_SOC, eigenval, eigenval_spinor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: PDOS_G0W0, PDOS_G0W0_SOC, PDOS_scf, &
                                                            PDOS_scf_SOC, proj_mo_on_kind, &
                                                            proj_mo_on_kind_spinor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_G0W0_2d, LDOS_scf_2d, &
                                                            LDOS_scf_2d_SOC
      TYPE(band_edges_type)                              :: band_edges_G0W0, band_edges_G0W0_SOC, &
                                                            band_edges_scf, band_edges_scf_guess, &
                                                            band_edges_scf_SOC
      TYPE(cp_cfm_type) :: cfm_ks_ikp, cfm_ks_ikp_spinor, cfm_mos_ikp, cfm_mos_ikp_spinor, &
         cfm_s_ikp, cfm_s_ikp_copy, cfm_s_ikp_spinor, cfm_s_ikp_spinor_copy, cfm_work_ikp, &
         cfm_work_ikp_spinor

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      homo = bs_env%n_occ(1)

      energy_window_DOS = bs_env%energy_window_DOS
      energy_step_DOS = bs_env%energy_step_DOS
      broadening = bs_env%broadening_DOS

      ! if we have done GW, we already have the band edges
      IF (bs_env%do_gw) THEN
         band_edges_scf = bs_env%band_edges_scf
         band_edges_scf_guess = band_edges_scf
      ELSE
         band_edges_scf_guess%VBM = bs_env%eigenval_scf_Gamma(homo)
         band_edges_scf_guess%CBM = bs_env%eigenval_scf_Gamma(homo + 1)

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
      ALLOCATE (proj_mo_on_kind_spinor(2*n_ao, nkind))
      proj_mo_on_kind_spinor(:, :) = 0.0_dp

      ALLOCATE (eigenval(n_ao))
      ALLOCATE (eigenval_spinor(2*n_ao))

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

      CALL cp_cfm_create(cfm_mos_ikp, bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_cfm_create(cfm_work_ikp, bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_cfm_create(cfm_s_ikp_copy, bs_env%fm_ks_Gamma(1)%matrix_struct)

      IF (bs_env%do_soc) THEN

         CALL cp_cfm_create(cfm_mos_ikp_spinor, bs_env%cfm_ks_spinor_ao_Gamma%matrix_struct)
         CALL cp_cfm_create(cfm_work_ikp_spinor, bs_env%cfm_ks_spinor_ao_Gamma%matrix_struct)
         CALL cp_cfm_create(cfm_s_ikp_spinor_copy, bs_env%cfm_ks_spinor_ao_Gamma%matrix_struct)
         CALL cp_cfm_create(cfm_ks_ikp_spinor, bs_env%cfm_ks_spinor_ao_Gamma%matrix_struct)
         CALL cp_cfm_create(cfm_s_ikp_spinor, bs_env%cfm_ks_spinor_ao_Gamma%matrix_struct)

         homo_spinor = 2*bs_env%n_occ(1)

         band_edges_scf_SOC%VBM = -1000.0_dp
         band_edges_scf_SOC%CBM = 1000.0_dp
         band_edges_scf_SOC%DBG = 1000.0_dp

         ! JW for later: G0W0+SOC
         IF (bs_env%do_gw) THEN
            band_edges_G0W0_SOC%VBM = -1000.0_dp
            band_edges_G0W0_SOC%CBM = 1000.0_dp
            band_edges_G0W0_SOC%DBG = 1000.0_dp
         END IF

      END IF

      IF (bs_env%do_local_vbm_cbm_gap) THEN
         CPASSERT(bs_env%int_ldos_xyz == int_ldos_z)
      END IF

      IF (bs_env%do_bs) THEN
         n_primitive_cells = bs_env%n_primitive_cells
      ELSE
         n_primitive_cells = 1
      END IF

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ''
      END IF

      DO ikp = 1, bs_env%kpoints_DOS%nkp

         bs_env%t1 = m_walltime()

         ! 1. get H^KS_µν(k_i) from H^KS_µν(k=0)
         CALL cfm_ikp_from_fm_Gamma(cfm_ks_ikp, bs_env%fm_ks_Gamma(1), &
                                    ikp, qs_env, bs_env%kpoints_DOS, "ORB")

         ! 2. get S_µν(k_i) from S_µν(k=0)
         CALL cfm_ikp_from_fm_Gamma(cfm_s_ikp, bs_env%fm_s_Gamma, &
                                    ikp, qs_env, bs_env%kpoints_DOS, "ORB")
         CALL cp_cfm_to_cfm(cfm_s_ikp, cfm_s_ikp_copy)

         ! 3. Diagonalize (Roothaan-Hall): H_KS(k_i)*C(k_i) = S(k_i)*C(k_i)*ϵ(k_i)
         CALL cp_cfm_geeig(cfm_ks_ikp, cfm_s_ikp_copy, cfm_mos_ikp, eigenval, cfm_work_ikp)

         ! 4. Projection p_nk^A of MO ψ_nk(r) on atom type A (inspired by Mulliken charge)
         !    p_nk^A = sum_µ^A,ν C*_µ^A,n(k) S_µ^A,ν(k) C_ν,n(k)
         CALL compute_proj_mo_on_kind(proj_mo_on_kind, qs_env, cfm_mos_ikp, cfm_s_ikp)

         ! 5. DOS and PDOS

         ! normalize to the primitive cell
         wkp = bs_env%kpoints_DOS%wkp(ikp)/n_primitive_cells
         DO i_E = 1, n_E
            energy = E_min + i_E*energy_step_DOS
            DO i_mo = 1, n_ao
               ! DOS, factor 2 for closed shell
               DOS_scf(i_E) = DOS_scf(i_E) + 2.0_dp*wkp*Gaussian(energy - eigenval(i_mo), &
                                                                 broadening)
               IF (bs_env%do_gw) THEN
                  ev = bs_env%eigenval_G0W0(i_mo, ikp)
                  DOS_G0W0(i_E) = DOS_G0W0(i_E) + 2.0_dp*wkp*Gaussian(energy - ev, broadening)
               END IF
               ! PDOS
               DO i_kind = 1, nkind
                  IF (proj_mo_on_kind(i_mo, i_kind) > 0.0_dp) THEN
                     PDOS_scf(i_E, i_kind) = PDOS_scf(i_E, i_kind) + &
                                             2.0_dp*proj_mo_on_kind(i_mo, i_kind)*wkp* &
                                             Gaussian(energy - eigenval(i_mo), broadening)

                     IF (bs_env%do_gw) THEN
                        PDOS_G0W0(i_E, i_kind) = PDOS_scf(i_E, i_kind) + &
                                                 2.0_dp*proj_mo_on_kind(i_mo, i_kind)*wkp* &
                                                 Gaussian(energy - ev, broadening)
                     END IF
                  END IF
               END DO
            END DO
         END DO

         IF (bs_env%do_local_vbm_cbm_gap) THEN
            CALL add_to_LDOS_2d(LDOS_scf_2d, qs_env, ikp, bs_env, cfm_mos_ikp, &
                                eigenval(:), band_edges_scf_guess)
            IF (bs_env%do_gw) THEN
               CALL add_to_LDOS_2d(LDOS_G0W0_2d, qs_env, ikp, bs_env, cfm_mos_ikp, &
                                   bs_env%eigenval_G0W0(:, ikp), band_edges_G0W0)
            END IF
         END IF

         band_edges_scf%VBM = MAX(band_edges_scf%VBM, eigenval(homo))
         band_edges_scf%CBM = MIN(band_edges_scf%CBM, eigenval(homo + 1))
         band_edges_scf%DBG = MIN(band_edges_scf%DBG, eigenval(homo + 1) - eigenval(homo))

         ! now the same with spin-orbit coupling
         IF (bs_env%do_soc) THEN

            ! 1. get H^KS_µν,σσ'(k_i) from H^KS_µν,σσ'(k=0)
            CALL cfm_ikp_from_cfm_spinor_Gamma(cfm_ks_ikp_spinor, bs_env%cfm_ks_spinor_ao_Gamma, &
                                               bs_env%fm_s_Gamma%matrix_struct, &
                                               ikp, qs_env, bs_env%kpoints_DOS, "ORB")

            ! 2. get S_µν,σσ'(k_i) from S_µν,σσ'(k=0)
            CALL cfm_ikp_from_cfm_spinor_Gamma(cfm_s_ikp_spinor, bs_env%cfm_s_spinor_Gamma, &
                                               bs_env%fm_s_Gamma%matrix_struct, &
                                               ikp, qs_env, bs_env%kpoints_DOS, "ORB")
            CALL cp_cfm_to_cfm(cfm_s_ikp_spinor, cfm_s_ikp_spinor_copy)

            ! 3. Diagonalize (Roothaan-Hall): H_KS(k_i)*C(k_i) = S(k_i)*C(k_i)*ϵ(k_i)
            CALL cp_cfm_geeig(cfm_ks_ikp_spinor, cfm_s_ikp_spinor_copy, cfm_mos_ikp_spinor, &
                              eigenval_spinor, cfm_work_ikp_spinor)

            ! 4. Projection p_nk^A of spinor ψ_nk(r) on atom type A (inspired by Mulliken charge)
            !    p_nk^A = sum_µ^A,ν,σ,σ' C*_µ^Aσ,n(k) S_µ^A,ν,σσ'(k) C_νσ',nσ(k)
            CALL compute_proj_mo_on_kind(proj_mo_on_kind_spinor, qs_env, &
                                         cfm_mos_ikp_spinor, cfm_s_ikp_spinor)

            ! 5. DOS and PDOS from spinors
            DO i_E = 1, n_E
               energy = E_min + i_E*energy_step_DOS
               DO i_mo = 1, 2*n_ao
                  ! DOS from spinors
                  DOS_scf_SOC(i_E) = DOS_scf_SOC(i_E) + wkp*Gaussian(energy - eigenval_spinor(i_mo), &
                                                                     broadening)
                  ! PDOS from spinors
                  DO i_kind = 1, nkind
                     IF (proj_mo_on_kind_spinor(i_mo, i_kind) > 0.0_dp) THEN
                        PDOS_scf_SOC(i_E, i_kind) = PDOS_scf_SOC(i_E, i_kind) + &
                                                    proj_mo_on_kind_spinor(i_mo, i_kind)*wkp* &
                                                    Gaussian(energy - eigenval_spinor(i_mo), &
                                                             broadening)
                     END IF
                  END DO
               END DO
            END DO

            ! 6. VBM and CBM
            band_edges_scf_SOC%VBM = MAX(band_edges_scf_SOC%VBM, eigenval_spinor(homo_spinor))
            band_edges_scf_SOC%CBM = MIN(band_edges_scf_SOC%CBM, eigenval_spinor(homo_spinor + 1))
            band_edges_scf_SOC%DBG = MIN(band_edges_scf_SOC%DBG, eigenval_spinor(homo_spinor + 1) &
                                         - eigenval_spinor(homo_spinor))

            CALL write_SOC_eigenvalues(eigenval_spinor, "SCF", ikp, bs_env)

            IF (bs_env%do_local_vbm_cbm_gap) THEN
               CALL add_to_LDOS_2d(LDOS_scf_2d_SOC, qs_env, ikp, bs_env, cfm_mos_ikp_spinor, &
                                   eigenval_spinor, band_edges_scf_guess, .TRUE., cfm_work_ikp)
            END IF

         END IF ! do_soc

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T43,I5,A,I3,A,F7.1,A)') &
               'Compute DOS, LDOS for k-point ', ikp, ' /', bs_env%kpoints_DOS%nkp, &
               ',    Execution time', m_walltime() - bs_env%t1, ' s'
         END IF

      END DO ! ikp_DOS

      band_edges_scf%IDBG = band_edges_scf%CBM - band_edges_scf%VBM
      IF (bs_env%do_soc) THEN
         band_edges_scf_SOC%IDBG = band_edges_scf_SOC%CBM - band_edges_scf_SOC%VBM
      END IF

      CALL write_band_edges(band_edges_scf, "SCF", bs_env)
      CALL write_dos_pdos(DOS_scf, PDOS_scf, bs_env, qs_env, "SCF", E_min, band_edges_scf%VBM)
      IF (bs_env%do_local_vbm_cbm_gap) THEN
         CALL compute_print_local_VBM_CBM_gap(LDOS_scf_2d, bs_env, band_edges_scf, "SCF")
      END IF

      IF (bs_env%do_soc) THEN
         CALL write_band_edges(band_edges_scf_SOC, "SCF+SOC", bs_env)
         CALL write_dos_pdos(DOS_scf_SOC, PDOS_scf_SOC, bs_env, qs_env, "SCF_SOC", &
                             E_min, band_edges_scf_SOC%VBM)
         IF (bs_env%do_local_vbm_cbm_gap) THEN
            ! argument band_edges_scf is actually correct because the non-SOC band edges
            ! have been used as reference in add_to_LDOS_2d
            CALL compute_print_local_VBM_CBM_gap(LDOS_scf_2d_SOC, bs_env, band_edges_scf, &
                                                 "SCF_SOC")
         END IF
      END IF

      IF (bs_env%do_gw) THEN
         CALL write_band_edges(band_edges_G0W0, "G0W0", bs_env)
         CALL write_dos_pdos(DOS_G0W0, PDOS_G0W0, bs_env, qs_env, "G0W0", E_min, &
                             band_edges_G0W0%VBM)
         IF (bs_env%do_local_vbm_cbm_gap) THEN
            CALL compute_print_local_VBM_CBM_gap(LDOS_G0W0_2d, bs_env, band_edges_G0W0, "G0W0")
         END IF
      END IF

      CALL cp_cfm_release(cfm_s_ikp)
      CALL cp_cfm_release(cfm_ks_ikp)
      CALL cp_cfm_release(cfm_mos_ikp)
      CALL cp_cfm_release(cfm_work_ikp)
      CALL cp_cfm_release(cfm_s_ikp_copy)

      CALL cp_cfm_release(cfm_s_ikp_spinor)
      CALL cp_cfm_release(cfm_ks_ikp_spinor)
      CALL cp_cfm_release(cfm_mos_ikp_spinor)
      CALL cp_cfm_release(cfm_work_ikp_spinor)
      CALL cp_cfm_release(cfm_s_ikp_spinor_copy)

      CALL timestop(handle)

   END SUBROUTINE dos_pdos_ldos

! **************************************************************************************************
!> \brief ...
!> \param LDOS_2d ...
!> \param bs_env ...
!> \param band_edges ...
!> \param scf_gw_soc ...
! **************************************************************************************************
   SUBROUTINE compute_print_local_VBM_CBM_gap(LDOS_2d, bs_env, band_edges, scf_gw_soc)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_2d
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(band_edges_type)                              :: band_edges
      CHARACTER(LEN=*)                                   :: scf_gw_soc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_print_local_VBM_CBM_gap'

      INTEGER :: handle, i_E, i_E_occ_end, i_E_occ_start, i_E_vir_end, i_E_vir_start, i_eps, &
         i_window, i_x, i_x_bin, i_x_end, i_x_end_bin, i_x_end_glob, i_x_start, i_x_start_bin, &
         i_x_start_glob, i_y, i_y_bin, i_y_end, i_y_end_bin, i_y_end_glob, i_y_start, &
         i_y_start_bin, i_y_start_glob, n_E, n_E_occ, n_E_vir, n_eps, n_windows
      INTEGER, DIMENSION(2)                              :: bin_mesh
      LOGICAL                                            :: do_xy_bins
      REAL(KIND=dp)                                      :: E_max_i_window, E_min, E_min_i_window, &
                                                            energy, energy_step_ldos, eps, &
                                                            max_energy_window_ldos, win
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: avg_LDOS_occ, avg_LDOS_vir, local_CBM, &
                                                            local_gap, local_VBM
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: LDOS_2d_bins

      CALL timeset(routineN, handle)

      n_eps = SIZE(bs_env%eps_ldos)
      n_windows = SIZE(bs_env%energy_window_ldos)
      n_E = SIZE(LDOS_2d, 3)

      max_energy_window_ldos = MAXVAL(bs_env%energy_window_ldos)
      energy_step_ldos = bs_env%energy_step_ldos
      E_min = band_edges%VBM - 0.5_dp*max_energy_window_ldos

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

      ALLOCATE (avg_LDOS_occ(i_x_start_bin:i_x_end_bin, i_y_start_bin:i_y_end_bin))
      ALLOCATE (avg_LDOS_vir(i_x_start_bin:i_x_end_bin, i_y_start_bin:i_y_end_bin))
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

         ! transform interval [i_x_start, i_x_end] to [1, bin_mesh(1)] (and same for y)
         DO i_x = i_x_start, i_x_end
            DO i_y = i_y_start, i_y_end
               i_x_bin = bin_mesh(1)*(i_x - i_x_start_glob)/(i_x_end_glob - i_x_start_glob + 1) + 1
               i_y_bin = bin_mesh(2)*(i_y - i_y_start_glob)/(i_y_end_glob - i_y_start_glob + 1) + 1
               LDOS_2d_bins(i_x_bin, i_y_bin, :) = LDOS_2d_bins(i_x_bin, i_y_bin, :) + &
                                                   LDOS_2d(i_x, i_y, :)
            END DO
         END DO

         CALL bs_env%para_env%sum(LDOS_2d_bins)

      ELSE

         LDOS_2d_bins(:, :, :) = LDOS_2d(:, :, :)

      END IF

      IF (bin_mesh(1)*bin_mesh(2) < bs_env%n_bins_max_for_printing) THEN
         CALL print_LDOS_2d_bins(LDOS_2d_bins, bs_env, E_min, scf_gw_soc)
      END IF

      ALLOCATE (local_VBM(i_x_start_bin:i_x_end_bin, i_y_start_bin:i_y_end_bin))
      ALLOCATE (local_CBM(i_x_start_bin:i_x_end_bin, i_y_start_bin:i_y_end_bin))
      ALLOCATE (local_gap(i_x_start_bin:i_x_end_bin, i_y_start_bin:i_y_end_bin))

      DO i_window = 1, n_windows

         win = bs_env%energy_window_ldos(i_window)

         i_E_occ_start = 1E7
         i_E_occ_end = -1E7
         i_E_vir_start = 1E7
         i_E_vir_end = -1E7

         DO i_E = 1, n_E

            energy = E_min + i_E*energy_step_ldos

            E_min_i_window = band_edges%VBM - win/2.0_dp
            E_max_i_window = band_edges%CBM + win/2.0_dp

            IF (energy >= E_min_i_window) i_E_occ_start = MIN(i_E_occ_start, i_E)
            IF (energy <= band_edges%VBM) i_E_occ_end = MAX(i_E_occ_end, i_E)
            IF (energy >= band_edges%CBM) i_E_vir_start = MIN(i_E_vir_start, i_E)
            IF (energy <= E_max_i_window) i_E_vir_end = MAX(i_E_vir_end, i_E)

         END DO

         n_E_occ = i_E_occ_end - i_E_occ_start + 1
         n_E_vir = i_E_vir_end - i_E_vir_start + 1

         DO i_x = i_x_start_bin, i_x_end_bin
            DO i_y = i_y_start_bin, i_y_end_bin
               avg_LDOS_occ(i_x, i_y) = SUM(LDOS_2d_bins(i_x, i_y, i_E_occ_start:i_E_occ_end)) &
                                        /n_E_occ
               avg_LDOS_vir(i_x, i_y) = SUM(LDOS_2d_bins(i_x, i_y, i_E_vir_start:i_E_vir_end)) &
                                        /n_E_vir
            END DO
         END DO

         DO i_eps = 1, n_eps

            local_VBM(:, :) = 0.0_dp
            local_CBM(:, :) = 0.0_dp
            local_gap(:, :) = 0.0_dp

            DO i_x = i_x_start_bin, i_x_end_bin
               DO i_y = i_y_start_bin, i_y_end_bin

                  DO i_E = n_E/2, 1, -1

                     eps = bs_env%eps_ldos(i_eps)*avg_LDOS_occ(i_x, i_y)

                     IF (LDOS_2d_bins(i_x, i_y, i_E) > eps) THEN
                        energy = E_min + i_E*energy_step_ldos
                        local_VBM(i_x, i_y) = energy
                        ! do not look for lower energies as we have found the local VBM
                        EXIT
                     END IF

                  END DO ! i_E

                  DO i_E = n_E/2 + 1, n_E

                     eps = bs_env%eps_ldos(i_eps)*avg_LDOS_vir(i_x, i_y)

                     IF (LDOS_2d_bins(i_x, i_y, i_E) > eps) THEN
                        energy = E_min + i_E*energy_step_ldos
                        local_CBM(i_x, i_y) = energy
                        ! do not look for higher energies as we have found the local CBM
                        EXIT
                     END IF

                  END DO ! i_E

                  local_gap(i_x, i_y) = local_CBM(i_x, i_y) - local_VBM(i_x, i_y)

               END DO ! i_y
            END DO ! i_x

            eps = bs_env%eps_ldos(i_eps)

            CALL sum_and_print_2d_map("VBM", scf_gw_soc, local_VBM, bs_env, eps, win, do_xy_bins)
            CALL sum_and_print_2d_map("CBM", scf_gw_soc, local_CBM, bs_env, eps, win, do_xy_bins)
            CALL sum_and_print_2d_map("gap", scf_gw_soc, local_gap, bs_env, eps, win, do_xy_bins)

         END DO ! i_eps

      END DO ! i_window

      CALL timestop(handle)

   END SUBROUTINE compute_print_local_VBM_CBM_gap

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
      CHARACTER(len=default_path_length)                 :: fname
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

!               IF (coord(1) < -1000/angstrom .AND. coord(2) < -1000/angstrom) THEN
!                  print_format = "(3A,F7.2,A,F7.2,A)"
!               ELSE IF (coord(1) < -100/angstrom .AND. coord(2) < -1000/angstrom) THEN
!                  print_format = "(3A,F6.2,A,F7.2,A)"
!               ELSE IF (coord(1) < -1000/angstrom .AND. coord(2) < -100/angstrom) THEN
!                  print_format = "(3A,F7.2,A,F6.2,A)"
!               ELSE IF (coord(1) < -100/angstrom .AND. coord(2) < -100/angstrom) THEN
!                  print_format = "(3A,F6.2,A,F6.2,A)"
!               ELSE IF (coord(1) < -100/angstrom .AND. coord(2) < -10/angstrom) THEN
!                  print_format = "(3A,F6.2,A,F5.2,A)"
!               ELSE IF (coord(1) < -10/angstrom .AND. coord(2) < -100/angstrom) THEN
!                  print_format = "(3A,F5.2,A,F6.2,A)"
!               ELSE IF (coord(1) < -10/angstrom .AND. coord(2) < -10/angstrom) THEN
!                  print_format = "(3A,F5.2,A,F5.2,A)"
!               ELSE IF (coord(1) < -1/angstrom .AND. coord(2) < -10/angstrom) THEN
!                  print_format = "(3A,F4.2,A,F5.2,A)"
!               ELSE IF (coord(1) < -10/angstrom .AND. coord(2) < -1/angstrom) THEN
!                  print_format = "(3A,F5.2,A,F4.2,A)"
!               ELSE IF (coord(1) < 10/angstrom .AND. coord(2) < 10/angstrom) THEN
!                  print_format = "(3A,F4.2,A,F4.2,A)"
!               ELSE IF (coord(1) < 100/angstrom .AND. coord(2) < 10/angstrom) THEN
!                  print_format = "(3A,F5.2,A,F4.2,A)"
!               ELSE IF (coord(1) < 10/angstrom .AND. coord(2) < 100/angstrom) THEN
!                  print_format = "(3A,F4.2,A,F5.2,A)"
!               ELSE IF (coord(1) < 100/angstrom .AND. coord(2) < 100/angstrom) THEN
!                  print_format = "(3A,F5.2,A,F5.2,A)"
!               ELSE IF (coord(1) < 100/angstrom .AND. coord(2) < 1000/angstrom) THEN
!                  print_format = "(3A,F5.2,A,F6.2,A)"
!               ELSE IF (coord(1) < 1000/angstrom .AND. coord(2) < 100/angstrom) THEN
!                  print_format = "(3A,F6.2,A,F5.2,A)"
!               ELSE
!                  print_format = "(3A,F6.2,A,F6.2,A)"
!               END IF

               WRITE (fname, print_format) "LDOS_", scf_gw_soc, &
                  "_at_x_", coord(1)*angstrom, '_A_and_y_', coord(2)*angstrom, '_A'

               CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", &
                              file_action="WRITE")

               DO i_E = 1, n_E
                  energy = E_min + i_E*bs_env%energy_step_ldos
                  WRITE (iunit, "(2F17.3)") energy*evolt, LDOS_2d_bins(i_x, i_y, i_E)*evolt
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
!> \param vbm_cbm_gap ...
!> \param scf_gw_soc ...
!> \param local_VBM_CBM_gap ...
!> \param bs_env ...
!> \param eps ...
!> \param win ...
!> \param do_xy_bins ...
! **************************************************************************************************
   SUBROUTINE sum_and_print_2d_map(vbm_cbm_gap, scf_gw_soc, local_VBM_CBM_gap, &
                                   bs_env, eps, win, do_xy_bins)
      CHARACTER(LEN=*)                                   :: vbm_cbm_gap, scf_gw_soc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: local_VBM_CBM_gap
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: eps, win
      LOGICAL                                            :: do_xy_bins

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sum_and_print_2d_map'

      CHARACTER(len=default_path_length)                 :: fname
      INTEGER :: handle, i_x, i_x_end, i_x_end_local, i_x_start, i_x_start_local, i_y, i_y_end, &
         i_y_end_local, i_y_start, i_y_start_local, iunit, n_x, n_y
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: local_VBM_CBM_gap_all
      REAL(KIND=dp), DIMENSION(3)                        :: coord, idx

      CALL timeset(routineN, handle)

      i_x_start_local = LBOUND(local_VBM_CBM_gap, 1)
      i_x_end_local = UBOUND(local_VBM_CBM_gap, 1)
      i_y_start_local = LBOUND(local_VBM_CBM_gap, 2)
      i_y_end_local = UBOUND(local_VBM_CBM_gap, 2)

      i_x_start = i_x_start_local
      i_x_end = i_x_end_local
      i_y_start = i_y_start_local
      i_y_end = i_y_end_local

      CALL bs_env%para_env%min(i_x_start)
      CALL bs_env%para_env%max(i_x_end)
      CALL bs_env%para_env%min(i_y_start)
      CALL bs_env%para_env%max(i_y_end)

      n_x = i_x_end - i_x_start + 1
      n_y = i_y_end - i_y_start + 1

      ALLOCATE (local_VBM_CBM_gap_all(i_x_start:i_x_end, i_y_start:i_y_end))
      local_VBM_CBM_gap_all(:, :) = 0.0_dp

      local_VBM_CBM_gap_all(i_x_start_local:i_x_end_local, i_y_start_local:i_y_end_local) = &
         local_VBM_CBM_gap(:, :)

      IF (.NOT. do_xy_bins) THEN

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(local_VBM_CBM_gap_all)
         CALL bs_env%para_env%sync()

      END IF

      IF (bs_env%para_env%is_source()) THEN

         WRITE (fname, "(4A,F5.3,A,F5.3)") vbm_cbm_gap, "_", scf_gw_soc, "_EPS_LDOS_", eps, &
            '_ENERGY_WINDOW_', win*evolt

         CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", file_action="WRITE")

         DO i_x = i_x_start, i_x_end
            DO i_y = i_y_start, i_y_end
               idx(1) = (REAL(i_x, KIND=dp) - 0.5_dp)/REAL(n_x, KIND=dp)
               idx(2) = (REAL(i_y, KIND=dp) - 0.5_dp)/REAL(n_y, KIND=dp)
               idx(3) = 0.0_dp
               coord(1:3) = MATMUL(bs_env%hmat, idx)
               WRITE (iunit, "(3F17.3)") coord(1)*angstrom, coord(2)*angstrom, &
                  local_VBM_CBM_gap_all(i_x, i_y)*evolt
            END DO
         END DO

         CALL close_file(iunit)

      END IF

      CALL timestop(handle)

   END SUBROUTINE sum_and_print_2d_map

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

      INTEGER :: handle, i_E, i_x, i_x_end, i_x_start, i_y, i_y_end, i_y_start, i_z, i_z_end, &
         i_z_start, j_col, j_mo, n_E, n_mo, ncol_local, nimages
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      LOGICAL                                            :: is_any_weight_non_zero, my_do_spinor
      REAL(KIND=dp)                                      :: arg, E_max, E_min, E_total_window, &
                                                            energy, energy_step_ldos, &
                                                            max_energy_window_ldos, weight, wkp
      TYPE(cp_cfm_type)                                  :: cfm_weighted_dm_ikp, cfm_work
      TYPE(cp_fm_type)                                   :: fm_non_spinor, fm_non_spinor_2, &
                                                            fm_weighted_dm_ikp_re
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: weighted_dm_ikp_re
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: LDOS_3d, rho_g
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      my_do_spinor = .FALSE.
      IF (PRESENT(do_spinor)) my_do_spinor = do_spinor

      CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env, dft_control=dft_control)

      ! previously, dft_control%nimages set to # neighbor cells, revert for Γ-only KS matrix
      nimages = dft_control%nimages
      dft_control%nimages = 1

      max_energy_window_ldos = MAXVAL(bs_env%energy_window_ldos)
      energy_step_ldos = bs_env%energy_step_ldos

      E_min = band_edges%VBM - 0.5_dp*max_energy_window_ldos
      E_max = band_edges%CBM + 0.5_dp*max_energy_window_ldos
      E_total_window = E_max - E_min

      n_E = INT(E_total_window/energy_step_ldos)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL auxbas_pw_pool%create_pw(LDOS_3d, use_data=REALDATA3D, in_space=REALSPACE)
      CALL auxbas_pw_pool%create_pw(rho_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      i_x_start = LBOUND(LDOS_3d%cr3d, 1)
      i_x_end = UBOUND(LDOS_3d%cr3d, 1)
      i_y_start = LBOUND(LDOS_3d%cr3d, 2)
      i_y_end = UBOUND(LDOS_3d%cr3d, 2)
      i_z_start = LBOUND(LDOS_3d%cr3d, 3)
      i_z_end = UBOUND(LDOS_3d%cr3d, 3)

      IF (ikp == 1) THEN
         ALLOCATE (LDOS_2d(i_x_start:i_x_end, i_y_start:i_y_end, n_E))
         LDOS_2d(:, :, :) = 0.0_dp
      END IF

      CALL cp_cfm_create(cfm_work, cfm_mos_ikp%matrix_struct)
      CALL cp_cfm_create(cfm_weighted_dm_ikp, cfm_mos_ikp%matrix_struct)
      CALL cp_fm_create(fm_weighted_dm_ikp_re, cfm_mos_ikp%matrix_struct)
      IF (my_do_spinor) THEN
         CALL cp_fm_create(fm_non_spinor, cfm_non_spinor%matrix_struct)
         CALL cp_fm_create(fm_non_spinor_2, cfm_non_spinor%matrix_struct)
      END IF

      CALL cp_cfm_get_info(matrix=cfm_mos_ikp, &
                           ncol_global=n_mo, &
                           ncol_local=ncol_local, &
                           col_indices=col_indices)

      NULLIFY (weighted_dm_ikp_re)
      CALL dbcsr_allocate_matrix_set(weighted_dm_ikp_re, 1)
      ALLOCATE (weighted_dm_ikp_re(1)%matrix)
      CALL dbcsr_create(weighted_dm_ikp_re(1)%matrix, template=bs_env%mat_ao_ao%matrix, &
                        matrix_type=dbcsr_type_symmetric)

      wkp = 0
      DO i_x = -bs_env%periodic(1), bs_env%periodic(1)
         DO i_y = -bs_env%periodic(2), bs_env%periodic(2)
            DO i_z = -bs_env%periodic(3), bs_env%periodic(3)
               arg = i_x*bs_env%kpoints_DOS%xkp(1, ikp) + &
                     i_y*bs_env%kpoints_DOS%xkp(2, ikp) + &
                     i_z*bs_env%kpoints_DOS%xkp(3, ikp)
               ! only cosine weight necessary; sine weight will cancel because of R <-> -R
               wkp = wkp + bs_env%kpoints_DOS%wkp(ikp)*COS(twopi*arg)
            END DO
         END DO
      END DO

      DO i_E = 1, n_E

         energy = E_min + i_E*energy_step_ldos

         is_any_weight_non_zero = .FALSE.

         DO j_col = 1, ncol_local

            j_mo = col_indices(j_col)

            weight = wkp*Gaussian(energy - eigenval(j_mo), energy_step_ldos)

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
               CALL get_cfm_submat(cfm_non_spinor, cfm_weighted_dm_ikp, 1, 1)
               CALL cp_cfm_to_fm(cfm_non_spinor, fm_non_spinor)
               CALL get_cfm_submat(cfm_non_spinor, cfm_weighted_dm_ikp, n_mo/2, n_mo/2)
               CALL cp_cfm_to_fm(cfm_non_spinor, fm_non_spinor_2)
               CALL cp_fm_scale_and_add(1.0_dp, fm_non_spinor, 1.0_dp, fm_non_spinor_2)
               CALL copy_fm_to_dbcsr(fm_non_spinor, weighted_dm_ikp_re(1)%matrix, &
                                     keep_sparsity=.FALSE.)
            ELSE
               CALL cp_cfm_to_fm(cfm_weighted_dm_ikp, fm_weighted_dm_ikp_re)
               CALL copy_fm_to_dbcsr(fm_weighted_dm_ikp_re, weighted_dm_ikp_re(1)%matrix, &
                                     keep_sparsity=.FALSE.)
            END IF

            LDOS_3d%cr3d(:, :, :) = 0.0_dp

            CALL calculate_rho_elec(matrix_p_kp=weighted_dm_ikp_re, &
                                    rho=LDOS_3d, &
                                    rho_gspace=rho_g, &
                                    ks_env=ks_env)

            DO i_z = i_z_start, i_z_end
               LDOS_2d(:, :, i_E) = LDOS_2d(:, :, i_E) + LDOS_3d%cr3d(:, :, i_z)
            END DO

         END IF

      END DO

      ! set back nimages
      dft_control%nimages = nimages

      CALL auxbas_pw_pool%give_back_pw(LDOS_3d)
      CALL auxbas_pw_pool%give_back_pw(rho_g)

      CALL cp_cfm_release(cfm_work)
      CALL cp_cfm_release(cfm_weighted_dm_ikp)

      CALL cp_fm_release(fm_weighted_dm_ikp_re)

      CALL dbcsr_deallocate_matrix_set(weighted_dm_ikp_re)

      IF (my_do_spinor) THEN
         CALL cp_fm_release(fm_non_spinor)
         CALL cp_fm_release(fm_non_spinor_2)
      END IF

      CALL timestop(handle)

   END SUBROUTINE add_to_LDOS_2d

! **************************************************************************************************
!> \brief ...
!> \param eigenval_spinor ...
!> \param scf_gw ...
!> \param ikp ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE write_SOC_eigenvalues(eigenval_spinor, scf_gw, ikp, bs_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_spinor
      CHARACTER(LEN=*)                                   :: scf_gw
      INTEGER                                            :: ikp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_SOC_eigenvalues'

      CHARACTER(len=3)                                   :: occ_vir
      CHARACTER(LEN=default_string_length)               :: fname
      INTEGER                                            :: handle, i_mo, iunit

      CALL timeset(routineN, handle)

      IF (ikp < 10) THEN
         WRITE (fname, "(2A,I1)") TRIM(scf_gw), "+SOC_band_structure_for_kpoint_000", ikp
      ELSE IF (ikp < 100) THEN
         WRITE (fname, "(2A,I2)") TRIM(scf_gw), "+SOC_band_structure_for_kpoint_00", ikp
      ELSE IF (ikp < 1000) THEN
         WRITE (fname, "(2A,I3)") TRIM(scf_gw), "+SOC_band_structure_for_kpoint_0", ikp
      ELSE
         WRITE (fname, "(2A,I4)") TRIM(scf_gw), "+SOC_band_structure_for_kpoint_", ikp
      END IF

      IF (bs_env%para_env%is_source()) THEN

         CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", file_action="WRITE")

         WRITE (iunit, "(A)") " "
         WRITE (iunit, "(A10,3F10.4)") "kpoint: ", bs_env%kpoints_DOS%xkp(:, ikp)
         WRITE (iunit, "(A)") " "

         DO i_mo = 1, SIZE(eigenval_spinor)
            IF (i_mo .LE. 2*bs_env%n_occ(1)) occ_vir = 'occ'
            IF (i_mo > 2*bs_env%n_occ(1)) occ_vir = 'vir'
            WRITE (iunit, "(I5,3A,4F16.3,F17.3)") i_mo, ' (', occ_vir, ') ', &
               eigenval_spinor(i_mo)*evolt
         END DO

         CALL close_file(iunit)

      END IF

      CALL timestop(handle)

   END SUBROUTINE write_SOC_eigenvalues

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

         WRITE (iunit, TRIM(output_string)) "Energy (eV)    DOS (1/eV)    PDOS (1/eV) ", &
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

      IF (energy < 5*broadening) THEN
         Gaussian = 1/broadening/SQRT(twopi)*EXP(-energy**2/broadening**2)
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

END MODULE post_scf_bandstructure_utils

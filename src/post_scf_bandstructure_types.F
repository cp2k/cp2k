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
MODULE post_scf_bandstructure_types
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cp_cfm_types,                    ONLY: cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_release
   USE dbt_api,                         ONLY: dbt_destroy,&
                                              dbt_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_release,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE qs_tensors_types,                ONLY: neighbor_list_3c_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'post_scf_bandstructure_types'

   PUBLIC :: post_scf_bandstructure_type, band_edges_type, bs_env_release

   ! valence band maximum (VBM), conduction band minimum (CBM), direct band gap (DBG),
   ! indirect band gap (IDBG)
   TYPE band_edges_type
      REAL(KIND=dp)                                   :: VBM = -1.0_dp, &
                                                         CBM = -1.0_dp, &
                                                         DBG = -1.0_dp, &
                                                         IDBG = -1.0_dp
   END TYPE band_edges_type

   TYPE post_scf_bandstructure_type

      ! decide which calculations will be done
      LOGICAL                                         :: do_gw = .FALSE., &
                                                         do_soc = .FALSE., &
                                                         do_bs = .FALSE., &
                                                         do_bs_primitive_cell = .FALSE., &
                                                         do_local_vbm_cbm_gap = .FALSE.

      ! various eigenvalues computed in GW code, some depend on k-points and
      ! have therefore two dimensions (spin dependence yet to do)
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: eigenval_scf_Gamma
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: eigenval_scf, &
                                                         eigenval_G0W0
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE  :: eigenval_scGW0
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: eigenval_prim_cell_scf, &
                                                         eigenval_prim_cell_scf_soc, &
                                                         eigenval_prim_cell_G0W0, &
                                                         eigenval_prim_cell_G0W0_soc
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: eigenval_scf_soc, &
                                                         eigenval_G0W0_soc
      TYPE(band_edges_type)                           :: band_edges_scf, &
                                                         band_edges_G0W0

      ! general parameters on molecular orbitals and basis sets
      INTEGER                                         :: n_ao = -1, &
                                                         n_RI = -1, &
                                                         n_spin = -1, &
                                                         n_atom = -1, &
                                                         max_AO_bf_per_atom = -1
      INTEGER, DIMENSION(:), ALLOCATABLE              :: i_ao_start_from_atom, &
                                                         i_ao_end_from_atom, &
                                                         i_RI_start_from_atom, &
                                                         i_RI_end_from_atom
      INTEGER, DIMENSION(2)                           :: n_occ = -1, &
                                                         n_vir = -1
      REAL(KIND=dp)                                   :: spin_degeneracy = -1.0_dp
      REAL(KIND=dp), DIMENSION(2)                     :: e_fermi = -1.0_dp

      ! kpoint mesh for chi, eps, W
      INTEGER, DIMENSION(:), POINTER                  :: nkp_grid_DOS_input => NULL()
      INTEGER, DIMENSION(3)                           :: nkp_grid_chi_eps_W_orig = -1, &
                                                         nkp_grid_chi_eps_W_extra = -1
      INTEGER                                         :: nkp_chi_eps_W_orig = -1, &
                                                         nkp_chi_eps_W_extra = -1, &
                                                         nkp_chi_eps_W_orig_plus_extra = -1, &
                                                         nkp_chi_eps_W_batch = -1, &
                                                         num_chi_eps_W_batches = -1, &
                                                         size_lattice_sum_V = -1
      TYPE(kpoint_type), POINTER                      :: kpoints_chi_eps_W => NULL(), &
                                                         kpoints_DOS => NULL(), &
                                                         kpoints_bandstructure => NULL()
      LOGICAL                                         :: approx_kp_extrapol = .FALSE.
      REAL(KIND=dp)                                   :: wkp_no_extra = -1.0_dp
      INTEGER                                         :: input_kp_bs_npoints = -1, &
                                                         input_kp_bs_n_sp_pts = -1, &
                                                         nkp_bs = -1, &
                                                         nkp_DOS = -1
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: xkp_special
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), ALLOCATABLE                    :: kp_special_name

      ! parameters for GW band structure calculation of small unit cell (with multiple unit cell)
      INTEGER, DIMENSION(3)                           :: periodic = -1, &
                                                         multiple_unit_cell = -1
      LOGICAL                                         :: calculate_bandstructure_of_primitive_cell &
                                                         = .FALSE.
      INTEGER                                         :: n_atom_in_primitive_cell = -1, &
                                                         n_primitive_cells = -1
      INTEGER, DIMENSION(:), ALLOCATABLE              :: atoms_i_primitive_cell, &
                                                         ref_atom_primitive_cell
      INTEGER, DIMENSION(:, :), ALLOCATABLE           :: cell_of_i_atom
      REAL(KIND=dp), DIMENSION(3, 3)                  :: hmat_primitive_cell = -1.0_dp, &
                                                         hinv_primitive_cell = -1.0_dp, &
                                                         hmat = -1.0_dp

      ! imaginary time and frequency grids
      INTEGER                                         :: num_time_freq_points = -1, &
                                                         num_freq_points_fit = -1
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: imag_time_points, &
                                                         imag_freq_points, &
                                                         imag_freq_points_fit
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: weights_cos_t_to_w, &
                                                         weights_cos_w_to_t, &
                                                         weights_sin_t_to_w
      INTEGER                                         :: nparam_pade = -1, &
                                                         num_points_per_magnitude = -1
      REAL(KIND=dp)                                   :: freq_max_fit = -1.0_dp, &
                                                         regularization_minimax = -1.0_dp, &
                                                         stabilize_exp = -1.0_dp

      ! filter threshold for matrix-tensor operations
      REAL(KIND=dp)                                   :: eps_filter = -1.0_dp, &
                                                         eps_3c_int = -1.0_dp, &
                                                         eps_atom_grid_2d_mat = -1.0_dp

      ! threshold for inverting ao overlap matrix, RI matrices
      REAL(KIND=dp)                                   :: eps_eigval_mat_s = -1.0_dp, &
                                                         eps_eigval_mat_RI = -1.0_dp, &
                                                         regularization_RI = -1.0_dp

      ! global full matrices used in GW
      TYPE(cp_fm_type)                                :: fm_s_Gamma
      TYPE(cp_fm_type), DIMENSION(2)                  :: fm_ks_Gamma, &
                                                         fm_V_xc_Gamma, &
                                                         fm_mo_coeff_Gamma
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE     :: fm_Gocc, &
                                                         fm_Gvir
      TYPE(cp_fm_type), DIMENSION(4)                  :: fm_work_mo
      TYPE(cp_fm_type)                                :: fm_RI_RI, &
                                                         fm_chi_Gamma_freq, &
                                                         fm_W_MIC_freq, &
                                                         fm_W_MIC_freq_1_extra, &
                                                         fm_W_MIC_freq_1_no_extra, &
                                                         fm_h_G0W0_Gamma

      ! global dbcsr matrices used in GW
      TYPE(dbcsr_p_type)                              :: mat_ao_ao, &
                                                         mat_RI_RI
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: mat_chi_Gamma_tau => NULL()

      ! local dbcsr matrices used in GW (local in tensor group)
      TYPE(dbcsr_p_type)                              :: mat_ao_ao_tensor, &
                                                         mat_RI_RI_tensor, &
                                                         mat_Sigma_from_Gocc_tensor, &
                                                         mat_Sigma_from_Gvir_tensor, &
                                                         mat_W_MIC_time_tensor

      TYPE(dbt_type)                                  :: t_G, &
                                                         t_chi, &
                                                         t_W, &
                                                         t_RI_AO__AO, &
                                                         t_RI__AO_AO

      ! parameters and data for parallelization
      INTEGER                                         :: group_size_tensor = -1, &
                                                         tensor_group_color = -1, &
                                                         num_tensor_groups = -1, &
                                                         diag_group_color = -1, &
                                                         num_diag_groups = -1, &
                                                         min_block_size = -1, &
                                                         input_memory_per_proc_GB = -1
      INTEGER(KIND=int_8)                             :: input_memory_per_proc = -1
      TYPE(mp_para_env_type), POINTER                 :: para_env => NULL(), &
                                                         para_env_tensor => NULL()
      REAL(KIND=dp)                                   :: occupation_3c_int = -1.0_dp, &
                                                         max_dist_AO_atoms = -1.0_dp, &
                                                         safety_factor_memory = -1.0_dp
      ! parallelization: atom range i and atom range j for tensor group
      INTEGER, DIMENSION(2)                           :: atoms_i = -1, &
                                                         atoms_j = -1
      INTEGER                                         :: n_atom_i = -1, &
                                                         n_intervals_i = -1, &
                                                         n_atom_j = -1, &
                                                         n_intervals_j = -1, &
                                                         n_atom_per_interval_ij = -1, &
                                                         n_intervals_inner_loop_atoms = -1, &
                                                         n_atom_per_IL_interval = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE           :: i_atom_intervals, &
                                                         j_atom_intervals, &
                                                         inner_loop_atom_intervals
      INTEGER, DIMENSION(:, :), ALLOCATABLE           :: atoms_i_t_group, &
                                                         atoms_j_t_group

      ! check-arrays and names for restarting
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: read_chi_and_G_occ_vir, &
                                                         calc_chi_and_G_occ_vir, &
                                                         Sigma_c_exists
      LOGICAL                                         :: all_W_exist = .FALSE., &
                                                         Sigma_x_exists = .FALSE.
      CHARACTER(LEN=3)                                :: chi_name = "chi"
      CHARACTER(LEN=5)                                :: Gocc_name = "G_occ", &
                                                         Gvir_name = "G_vir"
      CHARACTER(LEN=6)                                :: W_time_name = "W_time"
      CHARACTER(LEN=7)                                :: Sigma_x_name = "Sigma_x"
      CHARACTER(LEN=13)                               :: Sigma_p_name = "Sigma_pos_tau", &
                                                         Sigma_n_name = "Sigma_neg_tau"
      CHARACTER(LEN=default_path_length)              :: prefix = ""

      REAL(KIND=dp)                                   :: t1 = -1.0_dp, &
                                                         t2 = -1.0_dp
      INTEGER                                         :: unit_nr = -1

      ! parameters and data for basis sets
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                 :: basis_set_AO, &
                                                         basis_set_RI
      TYPE(neighbor_list_3c_type)                     :: nl_3c
      TYPE(libint_potential_type)                     :: ri_metric, &
                                                         trunc_coulomb

      ! parameters for SOC calculation (SOC is always computed by default because it is cheap)
      REAL(KIND=dp)                                   :: energy_window_soc = -1.0_dp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER    :: mat_V_SOC_xyz => NULL()
      TYPE(cp_fm_type), DIMENSION(3)                  :: fm_V_SOC_xyz_mo
      TYPE(cp_cfm_type)                               :: cfm_ks_spinor_ao_Gamma, &
                                                         cfm_h_G0W0_spinor_ao_Gamma, &
                                                         cfm_s_spinor_Gamma
      TYPE(band_edges_type)                           :: band_edges_scf_SOC, &
                                                         band_edges_G0W0_SOC
      ! parameters for DOS and PDOS calculation
      REAL(KIND=dp)                                   :: energy_window_DOS = -1.0_dp, &
                                                         energy_step_DOS = -1.0_dp, &
                                                         broadening_DOS = -1.0_dp

      ! parameters for local VBM, CBM, gap calculation (from local density of states, ldos)
      INTEGER                                         :: int_ldos_xyz = -1
      REAL(KIND=dp), DIMENSION(:), POINTER            :: energy_window_ldos => NULL(), &
                                                         eps_ldos => NULL()
      REAL(KIND=dp)                                   :: energy_step_ldos = -1.0_dp
      INTEGER, DIMENSION(:), POINTER                  :: bin_mesh => NULL()
      INTEGER                                         :: n_bins_max_for_printing = -1

   END TYPE post_scf_bandstructure_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE bs_env_release(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bs_env_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(bs_env))

      CALL safe_kpoints_release(bs_env%kpoints_chi_eps_W)
      CALL safe_kpoints_release(bs_env%kpoints_DOS)
      CALL safe_kpoints_release(bs_env%kpoints_bandstructure)

      IF (ALLOCATED(bs_env%xkp_special)) DEALLOCATE (bs_env%xkp_special)
      IF (ALLOCATED(bs_env%kp_special_name)) DEALLOCATE (bs_env%kp_special_name)
      IF (ALLOCATED(bs_env%atoms_i_primitive_cell)) DEALLOCATE (bs_env%atoms_i_primitive_cell)
      IF (ALLOCATED(bs_env%ref_atom_primitive_cell)) DEALLOCATE (bs_env%ref_atom_primitive_cell)
      IF (ALLOCATED(bs_env%cell_of_i_atom)) DEALLOCATE (bs_env%cell_of_i_atom)
      IF (ALLOCATED(bs_env%imag_time_points)) DEALLOCATE (bs_env%imag_time_points)
      IF (ALLOCATED(bs_env%imag_freq_points)) DEALLOCATE (bs_env%imag_freq_points)
      IF (ALLOCATED(bs_env%eigenval_scf_Gamma)) DEALLOCATE (bs_env%eigenval_scf_Gamma)
      IF (ALLOCATED(bs_env%eigenval_scf)) DEALLOCATE (bs_env%eigenval_scf)
      IF (ALLOCATED(bs_env%eigenval_G0W0)) DEALLOCATE (bs_env%eigenval_G0W0)
      IF (ALLOCATED(bs_env%eigenval_scGW0)) DEALLOCATE (bs_env%eigenval_scGW0)
      IF (ALLOCATED(bs_env%eigenval_prim_cell_scf)) DEALLOCATE (bs_env%eigenval_prim_cell_scf)
      IF (ALLOCATED(bs_env%eigenval_prim_cell_scf_soc)) DEALLOCATE (bs_env%eigenval_prim_cell_scf_soc)
      IF (ALLOCATED(bs_env%eigenval_prim_cell_G0W0)) DEALLOCATE (bs_env%eigenval_prim_cell_G0W0)
      IF (ALLOCATED(bs_env%eigenval_prim_cell_G0W0_soc)) DEALLOCATE (bs_env%eigenval_prim_cell_G0W0_soc)
      IF (ALLOCATED(bs_env%eigenval_scf_soc)) DEALLOCATE (bs_env%eigenval_scf_soc)
      IF (ALLOCATED(bs_env%eigenval_G0W0_soc)) DEALLOCATE (bs_env%eigenval_G0W0_soc)
      IF (ALLOCATED(bs_env%i_ao_start_from_atom)) DEALLOCATE (bs_env%i_ao_start_from_atom)
      IF (ALLOCATED(bs_env%i_ao_end_from_atom)) DEALLOCATE (bs_env%i_ao_end_from_atom)
      IF (ALLOCATED(bs_env%i_RI_start_from_atom)) DEALLOCATE (bs_env%i_RI_start_from_atom)
      IF (ALLOCATED(bs_env%i_RI_end_from_atom)) DEALLOCATE (bs_env%i_RI_end_from_atom)
      IF (ALLOCATED(bs_env%i_atom_intervals)) DEALLOCATE (bs_env%i_atom_intervals)
      IF (ALLOCATED(bs_env%j_atom_intervals)) DEALLOCATE (bs_env%j_atom_intervals)
      IF (ALLOCATED(bs_env%atoms_i_t_group)) DEALLOCATE (bs_env%atoms_i_t_group)
      IF (ALLOCATED(bs_env%atoms_j_t_group)) DEALLOCATE (bs_env%atoms_j_t_group)
      IF (ALLOCATED(bs_env%read_chi_and_G_occ_vir)) DEALLOCATE (bs_env%read_chi_and_G_occ_vir)
      IF (ALLOCATED(bs_env%calc_chi_and_G_occ_vir)) DEALLOCATE (bs_env%calc_chi_and_G_occ_vir)
      IF (ALLOCATED(bs_env%Sigma_c_exists)) DEALLOCATE (bs_env%Sigma_c_exists)

      CALL cp_fm_release(bs_env%fm_s_Gamma)
      CALL cp_fm_release(bs_env%fm_ks_Gamma(1))
      CALL cp_fm_release(bs_env%fm_ks_Gamma(2))
      CALL cp_fm_release(bs_env%fm_V_xc_Gamma(1))
      CALL cp_fm_release(bs_env%fm_V_xc_Gamma(2))
      CALL cp_fm_release(bs_env%fm_mo_coeff_Gamma(1))
      CALL cp_fm_release(bs_env%fm_mo_coeff_Gamma(2))
      CALL cp_fm_release(bs_env%fm_Gocc)
      CALL cp_fm_release(bs_env%fm_Gvir)
      CALL cp_fm_release(bs_env%fm_work_mo(1))
      CALL cp_fm_release(bs_env%fm_work_mo(2))
      CALL cp_fm_release(bs_env%fm_work_mo(3))
      CALL cp_fm_release(bs_env%fm_work_mo(4))
      CALL cp_fm_release(bs_env%fm_RI_RI)
      CALL cp_fm_release(bs_env%fm_chi_Gamma_freq)
      CALL cp_fm_release(bs_env%fm_W_MIC_freq)
      CALL cp_fm_release(bs_env%fm_W_MIC_freq_1_extra)
      CALL cp_fm_release(bs_env%fm_W_MIC_freq_1_no_extra)
      CALL cp_fm_release(bs_env%fm_h_G0W0_Gamma)

      CALL release_dbcsr_p_type(bs_env%mat_ao_ao)
      CALL release_dbcsr_p_type(bs_env%mat_RI_RI)
      CALL safe_dbcsr_deallocate_matrix_set_1d(bs_env%mat_chi_Gamma_tau)

      CALL release_dbcsr_p_type(bs_env%mat_ao_ao_tensor)
      CALL release_dbcsr_p_type(bs_env%mat_RI_RI_tensor)
      CALL release_dbcsr_p_type(bs_env%mat_Sigma_from_Gocc_tensor)
      CALL release_dbcsr_p_type(bs_env%mat_Sigma_from_Gvir_tensor)
      CALL release_dbcsr_p_type(bs_env%mat_W_MIC_time_tensor)

      CALL mp_para_env_release(bs_env%para_env)
      IF (ASSOCIATED(bs_env%para_env_tensor)) CALL mp_para_env_release(bs_env%para_env_tensor)

      CALL safe_dbt_destroy(bs_env%t_G)
      CALL safe_dbt_destroy(bs_env%t_chi)
      CALL safe_dbt_destroy(bs_env%t_W)
      CALL safe_dbt_destroy(bs_env%t_RI_AO__AO)
      CALL safe_dbt_destroy(bs_env%t_RI__AO_AO)

      IF (ALLOCATED(bs_env%basis_set_AO)) DEALLOCATE (bs_env%basis_set_AO)
      IF (ALLOCATED(bs_env%basis_set_RI)) DEALLOCATE (bs_env%basis_set_RI)

      ! SOC matrices and arrays
      CALL safe_dbcsr_deallocate_matrix_set_2d(bs_env%mat_V_SOC_xyz)
      CALL cp_fm_release(bs_env%fm_V_SOC_xyz_mo(1))
      CALL cp_fm_release(bs_env%fm_V_SOC_xyz_mo(2))
      CALL cp_fm_release(bs_env%fm_V_SOC_xyz_mo(3))
      CALL cp_cfm_release(bs_env%cfm_ks_spinor_ao_Gamma)
      CALL cp_cfm_release(bs_env%cfm_h_G0W0_spinor_ao_Gamma)
      CALL cp_cfm_release(bs_env%cfm_s_spinor_Gamma)

      DEALLOCATE (bs_env)

      CALL timestop(handle)

   END SUBROUTINE bs_env_release

! **************************************************************************************************
!> \brief ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE safe_kpoints_release(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      IF (ASSOCIATED(kpoints)) CALL kpoint_release(kpoints)

   END SUBROUTINE safe_kpoints_release

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_p_type_matrix ...
! **************************************************************************************************
   SUBROUTINE release_dbcsr_p_type(dbcsr_p_type_matrix)
      TYPE(dbcsr_p_type)                                 :: dbcsr_p_type_matrix

      IF (ASSOCIATED(dbcsr_p_type_matrix%matrix)) THEN
         CALL dbcsr_release(dbcsr_p_type_matrix%matrix)
         DEALLOCATE (dbcsr_p_type_matrix%matrix)
      END IF

   END SUBROUTINE release_dbcsr_p_type

! **************************************************************************************************
!> \brief ...
!> \param t ...
! **************************************************************************************************
   SUBROUTINE safe_dbt_destroy(t)
      TYPE(dbt_type)                                     :: t

      IF (ASSOCIATED(t%matrix_rep)) CALL dbt_destroy(t)

   END SUBROUTINE safe_dbt_destroy

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_p_type_matrix_array ...
! **************************************************************************************************
   SUBROUTINE safe_dbcsr_deallocate_matrix_set_1d(dbcsr_p_type_matrix_array)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_p_type_matrix_array

      IF (ASSOCIATED(dbcsr_p_type_matrix_array)) CALL dbcsr_deallocate_matrix_set(dbcsr_p_type_matrix_array)

   END SUBROUTINE safe_dbcsr_deallocate_matrix_set_1d
! **************************************************************************************************
!> \brief ...
!> \param dbcsr_p_type_matrix_array ...
! **************************************************************************************************
   SUBROUTINE safe_dbcsr_deallocate_matrix_set_2d(dbcsr_p_type_matrix_array)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: dbcsr_p_type_matrix_array

      IF (ASSOCIATED(dbcsr_p_type_matrix_array)) CALL dbcsr_deallocate_matrix_set(dbcsr_p_type_matrix_array)

   END SUBROUTINE safe_dbcsr_deallocate_matrix_set_2d

END MODULE post_scf_bandstructure_types

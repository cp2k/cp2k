!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for GW calculation of the electronic band structure
!> \par History
!>       2023.07 created [Jan Wilhelm]
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE gw_types
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dbt_api,                         ONLY: dbt_destroy,&
                                              dbt_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_release,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_tensors_types,                ONLY: neighbor_list_3c_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_types'

   PUBLIC :: gw_type, gw_env_release

   TYPE gw_type

      ! various eigenvalues computed in GW code, some depend on k-points and
      ! have therefore two dimensions (spin dependence yet to do)
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: eigenvalues_scf_Gamma
      ! JW TODO: fully distribute these at given point to save memory
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: eigenvalues_scf, &
                                                         eigenvalues_G0W0
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE  :: eigenvalues_scGW0
      REAL(KIND=dp)                                   :: E_VBM_scf, &
                                                         E_CBM_scf, &
                                                         E_DBG_scf, &
                                                         E_IDBG_scf, &
                                                         E_VBM_G0W0, &
                                                         E_CBM_G0W0, &
                                                         E_DBG_G0W0, &
                                                         E_IDBG_G0W0

      ! general parameters on molecular orbitals and basis sets
      INTEGER                                         :: n_ao, &
                                                         n_RI, &
                                                         n_spin, &
                                                         n_atom
      INTEGER, DIMENSION(2)                           :: n_occ, &
                                                         n_vir
      REAL(KIND=dp)                                   :: spin_degeneracy
      REAL(KIND=dp), DIMENSION(2)                     :: e_fermi

      ! kpoint mesh for chi, eps, W
      INTEGER, DIMENSION(:), POINTER                  :: nkp_grid_chi_eps_W_input, &
                                                         nkp_grid_Sigma_input
      INTEGER, DIMENSION(3)                           :: nkp_grid_chi_eps_W_orig, &
                                                         nkp_grid_chi_eps_W_extra
      INTEGER                                         :: nkp_chi_eps_W_orig, &
                                                         nkp_chi_eps_W_extra, &
                                                         ikp_local_start, &
                                                         ikp_local_end
      TYPE(kpoint_type), POINTER                      :: kpoints_chi_eps_W, &
                                                         kpoints_Sigma
      INTEGER                                         :: size_lattice_sum_V

      ! imaginary time and frequency grids
      INTEGER                                         :: num_time_freq_points, &
                                                         num_freq_points_fit
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: imag_time_points, &
                                                         imag_freq_points, &
                                                         imag_freq_points_fit
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE     :: weights_cos_t_to_w, &
                                                         weights_cos_w_to_t, &
                                                         weights_sin_t_to_w
      INTEGER                                         :: nparam_pade, &
                                                         num_points_per_magnitude
      REAL(KIND=dp)                                   :: freq_max_fit, &
                                                         regularization_minimax, &
                                                         stabilize_exp

      ! filter threshold for matrix-tensor operations
      REAL(KIND=dp)                                   :: eps_filter, &
                                                         eps_3c_int, &
                                                         eps_atom_grid_2d_mat

      ! threshold for inverting ao overlap matrix, RI matrices
      REAL(KIND=dp)                                   :: eps_eigval_mat_s, &
                                                         eps_eigval_mat_RI, &
                                                         regularization_RI

      ! global full matrices used in GW
      TYPE(cp_fm_type)                                :: fm_s_Gamma
      TYPE(cp_fm_type), DIMENSION(2)                  :: fm_ks_Gamma, &
                                                         fm_mo_coeff_Gamma
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE     :: fm_Gocc, &
                                                         fm_Gvir, &
                                                         fm_chi_Gamma_time
      TYPE(cp_fm_type), DIMENSION(3)                  :: fm_work_mo
! JW 2del
!      TYPE(cp_cfm_type), DIMENSION(:), POINTER        :: cfm_V_kp
!      TYPE(cp_fm_type), DIMENSION(:,:), ALLOCATABLE   :: fm_V_kp
!      TYPE(cp_fm_type)                                :: fm_RI_RI
!      TYPE(cp_cfm_type)                                :: cfm_RI_RI

      ! global dbcsr matrices used in GW
      TYPE(dbcsr_p_type)                              :: mat_ao_ao, &
                                                         mat_RI_RI
! JW 2del creates corrupted list
!      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: mat_chi_Gamma_time

      ! local dbcsr matrices used in GW (local in tensor group)
      TYPE(dbcsr_p_type)                              :: mat_ao_ao_tensor, &
                                                         mat_RI_RI_tensor, &
                                                         mat_Sigma_from_Gocc_tensor, &
                                                         mat_Sigma_from_Gvir_tensor, &
                                                         mat_W_MIC_time_tensor

      TYPE(dbt_type)                                  :: t_G, &
                                                         t_chi, &
                                                         t_W, &
                                                         t_Sigma, &
                                                         t_RI_AO__AO, &
                                                         t_RI__AO_AO, &
                                                         ! maybe for later: optimize for 3c operations
                                                         t_RI_AO__AO_opt

      ! local dbcsr matrices used in GW (local in kp diag group)
!      TYPE(dbcsr_type), POINTER                       :: mat_V_kp_loc, &
!                                                         mat_M_kp_loc, &
!                                                         mat_M_inv_kp_loc

      ! parameters and data for parallelization
      INTEGER                                         :: tensor_group_size, &
                                                         tensor_group_color, &
                                                         num_tensor_groups, &
                                                         kp_diag_group_size, &
                                                         kp_diag_group_color, &
                                                         num_kp_diag_groups, &
                                                         min_block_size
      TYPE(mp_para_env_type), POINTER                 :: para_env, &
                                                         para_env_tensor, &
                                                         para_env_kp_diag
      ! parallelization: two different atom ranges for tensor group
! JW 2del
!      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_subset_i, &
!                                                         particle_subset_j
      INTEGER                                         :: atom_i_start, &
                                                         atom_i_end, &
                                                         atom_j_start, &
                                                         atom_j_end, &
                                                         n_atom_ij

      REAL(KIND=dp)                                   :: t1, t2
      INTEGER                                         :: unit_nr

      ! parameters and data for basis sets
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                 :: basis_set_AO, &
                                                         basis_set_RI
      TYPE(neighbor_list_3c_type)                     :: nl_3c
      TYPE(libint_potential_type)                     :: ri_metric, &
                                                         trunc_coulomb

   END TYPE gw_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE gw_env_release(gw_env)
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gw_env_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(gw_env))

      IF (ASSOCIATED(gw_env%nkp_grid_chi_eps_W_input)) DEALLOCATE (gw_env%nkp_grid_chi_eps_W_input)
      IF (ASSOCIATED(gw_env%kpoints_chi_eps_W)) CALL kpoint_release(gw_env%kpoints_chi_eps_W)
      IF (ASSOCIATED(gw_env%kpoints_Sigma)) CALL kpoint_release(gw_env%kpoints_Sigma)

      IF (ALLOCATED(gw_env%imag_time_points)) DEALLOCATE (gw_env%imag_time_points)
      IF (ALLOCATED(gw_env%imag_freq_points)) DEALLOCATE (gw_env%imag_freq_points)
      IF (ALLOCATED(gw_env%eigenvalues_scf_Gamma)) DEALLOCATE (gw_env%eigenvalues_scf_Gamma)
      IF (ALLOCATED(gw_env%eigenvalues_scf)) DEALLOCATE (gw_env%eigenvalues_scf)
      IF (ALLOCATED(gw_env%eigenvalues_G0W0)) DEALLOCATE (gw_env%eigenvalues_G0W0)
      IF (ALLOCATED(gw_env%eigenvalues_scGW0)) DEALLOCATE (gw_env%eigenvalues_scGW0)

      CALL cp_fm_release(gw_env%fm_s_Gamma)
      CALL cp_fm_release(gw_env%fm_ks_Gamma(1))
      CALL cp_fm_release(gw_env%fm_ks_Gamma(2))
      CALL cp_fm_release(gw_env%fm_mo_coeff_Gamma(1))
      CALL cp_fm_release(gw_env%fm_mo_coeff_Gamma(2))
      CALL cp_fm_release(gw_env%fm_Gocc)
      CALL cp_fm_release(gw_env%fm_Gvir)
      CALL cp_fm_release(gw_env%fm_chi_Gamma_time)
      CALL cp_fm_release(gw_env%fm_work_mo(1))
      CALL cp_fm_release(gw_env%fm_work_mo(2))
      CALL cp_fm_release(gw_env%fm_work_mo(3))

      CALL release_dbcsr_p_type(gw_env%mat_ao_ao)
      CALL release_dbcsr_p_type(gw_env%mat_RI_RI)

      CALL release_dbcsr_p_type(gw_env%mat_ao_ao_tensor)
      CALL release_dbcsr_p_type(gw_env%mat_RI_RI_tensor)
      CALL release_dbcsr_p_type(gw_env%mat_Sigma_from_Gocc_tensor)
      CALL release_dbcsr_p_type(gw_env%mat_Sigma_from_Gvir_tensor)
      CALL release_dbcsr_p_type(gw_env%mat_W_MIC_time_tensor)

      CALL mp_para_env_release(gw_env%para_env)
      CALL mp_para_env_release(gw_env%para_env_tensor)
      CALL mp_para_env_release(gw_env%para_env_kp_diag)

      CALL safe_dbt_destroy(gw_env%t_G)
      CALL safe_dbt_destroy(gw_env%t_chi)
      CALL safe_dbt_destroy(gw_env%t_W)
      CALL safe_dbt_destroy(gw_env%t_Sigma)
      CALL safe_dbt_destroy(gw_env%t_RI_AO__AO)
      CALL safe_dbt_destroy(gw_env%t_RI__AO_AO)

      IF (ALLOCATED(gw_env%basis_set_AO)) DEALLOCATE (gw_env%basis_set_AO)
      IF (ALLOCATED(gw_env%basis_set_RI)) DEALLOCATE (gw_env%basis_set_RI)

      DEALLOCATE (gw_env)

      CALL timestop(handle)

   END SUBROUTINE gw_env_release

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

END MODULE gw_types

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX and K-points.
!> \auhtor Augustin Bussy (01.2023)
! **************************************************************************************************

MODULE hfx_ri_kp
   USE admm_types,                      ONLY: get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              real_to_scaled,&
                                              scaled_to_real
   USE cp_array_utils,                  ONLY: cp_1d_logical_p_type,&
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_clear, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_dot, &
        dbcsr_filter, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_put_block, dbcsr_release, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE dbt_api,                         ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_clear, dbt_contract, &
        dbt_copy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_create, dbt_destroy, &
        dbt_distribution_destroy, dbt_distribution_new, dbt_distribution_type, dbt_filter, &
        dbt_finalize, dbt_get_block, dbt_get_info, dbt_get_stored_coordinates, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_mp_environ_pgrid, dbt_pgrid_create, dbt_pgrid_destroy, &
        dbt_pgrid_type, dbt_put_block, dbt_scale, dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE hfx_ri,                          ONLY: get_idx_to_atom,&
                                              hfx_ri_pre_scf_calc_tensors
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: do_potential_short,&
                                              hfx_ri_do_2c_cholesky,&
                                              hfx_ri_do_2c_diag,&
                                              hfx_ri_do_2c_iter
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE iterate_matrix,                  ONLY: invert_hotelling
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: erfc_cutoff
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_para_env_type,&
                                              mp_request_type,&
                                              mp_waitall
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tensors,                      ONLY: &
        build_2c_derivatives, build_2c_neighbor_lists, build_3c_derivatives, &
        build_3c_neighbor_lists, get_3c_iterator_info, get_tensor_occupancy, &
        neighbor_list_3c_destroy, neighbor_list_3c_iterate, neighbor_list_3c_iterator_create, &
        neighbor_list_3c_iterator_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              create_tensor_batches,&
                                              distribution_2d_create,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_iterator_type,&
                                              neighbor_list_3c_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_num_threads

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_ri_update_ks_kp, hfx_ri_update_forces_kp

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_kp'
CONTAINS

! NOTES: for a start, we do not seek performance, but accuracy. So in this first implementation,
!        we give little consideration to batching, load balance and such.
!        We also put everything here, even if there is some code replication with the original RI_HFX
!        We will only work in the RHO flavor
!        For now, we will also always assume that there is a single para_env, and that there is no
!        K-point subgroup. This might change in the future

! **************************************************************************************************
!> \brief I_1nitialize the ri_data for K-point. For now, we take the normal, usual existing ri_data
!>        and we adapt it to our needs
!> \param dbcsr_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE adapt_ri_data_to_kp(dbcsr_template, ri_data, qs_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i_img, i_RI, i_spin, iatom, natom, &
                                                            nblks_RI, nimg, nkind, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, dist1, dist2, dist3
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (dft_control, para_env)

      !The main thing that we need to do is to allocate more space for the integrals, such that there
      !is room for each periodic image. Note that we only go in 1D, i.e. we store (mu^0 sigma^a|P^0),
      !and (P^0|Q^a) => the RI basis is always in the main cell.

      !Get kpoint info
      CALL get_qs_env(qs_env, dft_control=dft_control, natom=natom, para_env=para_env, nkind=nkind)
      nimg = ri_data%nimg

      !Along the RI direction we have basis elements spread accross ncell_RI images.
      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(nblks_RI*ri_data%ncell_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI_split(:)
      END DO

      ALLOCATE (ri_data%t_3c_int_ctr_1(1, nimg))
      CALL create_3c_tensor(ri_data%t_3c_int_ctr_1(1, 1), dist1, dist2, dist3, &
                            ri_data%pgrid_1, ri_data%bsizes_AO_split, bsizes_RI_ext, &
                            ri_data%bsizes_AO_split, [1, 2], [3], name="(AO RI | AO)")

      DO i_img = 2, nimg
         CALL dbt_create(ri_data%t_3c_int_ctr_1(1, 1), ri_data%t_3c_int_ctr_1(1, i_img))
      END DO
      DEALLOCATE (dist1, dist2, dist3)

      ALLOCATE (ri_data%t_3c_int_ctr_2(1, 1))
      CALL create_3c_tensor(ri_data%t_3c_int_ctr_2(1, 1), dist1, dist2, dist3, &
                            ri_data%pgrid_1, ri_data%bsizes_AO_split, bsizes_RI_ext, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      !We use full block sizes for the 2c quantities
      DEALLOCATE (bsizes_RI_ext)
      nblks_RI = SIZE(ri_data%bsizes_RI)
      ALLOCATE (bsizes_RI_ext(nblks_RI*ri_data%ncell_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI(:)
      END DO

      ALLOCATE (ri_data%t_2c_inv(1, natom), ri_data%t_2c_int(1, natom), ri_data%t_2c_pot(1, natom))
      CALL create_2c_tensor(ri_data%t_2c_inv(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                            bsizes_RI_ext, bsizes_RI_ext, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_create(ri_data%t_2c_inv(1, 1), ri_data%t_2c_int(1, 1))
      CALL dbt_create(ri_data%t_2c_inv(1, 1), ri_data%t_2c_pot(1, 1))
      DO iatom = 2, natom
         CALL dbt_create(ri_data%t_2c_inv(1, 1), ri_data%t_2c_inv(1, iatom))
         CALL dbt_create(ri_data%t_2c_inv(1, 1), ri_data%t_2c_int(1, iatom))
         CALL dbt_create(ri_data%t_2c_inv(1, 1), ri_data%t_2c_pot(1, iatom))
      END DO

      ALLOCATE (ri_data%kp_cost(natom, natom, nimg))
      ri_data%kp_cost = 0.0_dp

      !We store the density and KS matrix in tensor format
      nspins = dft_control%nspins
      ALLOCATE (ri_data%rho_ao_t(nspins, nimg), ri_data%ks_t(nspins, nimg))
      CALL create_2c_tensor(ri_data%rho_ao_t(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      CALL dbt_create(dbcsr_template, ri_data%ks_t(1, 1))

      IF (nspins == 2) THEN
         CALL dbt_create(ri_data%rho_ao_t(1, 1), ri_data%rho_ao_t(2, 1))
         CALL dbt_create(ri_data%ks_t(1, 1), ri_data%ks_t(2, 1))
      END IF
      DO i_img = 2, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ri_data%rho_ao_t(1, 1), ri_data%rho_ao_t(i_spin, i_img))
            CALL dbt_create(ri_data%ks_t(1, 1), ri_data%ks_t(i_spin, i_img))
         END DO
      END DO

   END SUBROUTINE adapt_ri_data_to_kp

! **************************************************************************************************
!> \brief The pre-scf steps for RI-HFX k-points calculation. Namely the calculation of the integrals
!> \param dbcsr_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_kp(dbcsr_template, ri_data, qs_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ri_pre_scf_kp'

      INTEGER                                            :: handle, i_img, iatom, natom, nimg, nkind
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: t_2c_op_pot, t_2c_op_RI
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_int
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, natom=natom, nkind=nkind)

      CALL cleanup_kp(ri_data)

      !We do all the checks on what we allow in this initial implementation
      IF (ri_data%flavor .NE. ri_pmat) CPABORT("K-points RI-HFX only with RHO flavor")
      IF (ri_data%same_op) ri_data%same_op = .FALSE. !force the full calculation with RI metric
      IF (ABS(ri_data%eps_pgf_orb - dft_control%qs_control%eps_pgf_orb) > 1.0E-16_dp) &
         CPABORT("RI%EPS_PGF_ORB and QS%EPS_PGF_ORB must be identical for RI-HFX k-points")

      CALL get_kp_and_ri_images(ri_data, qs_env)
      nimg = ri_data%nimg

      !Calculate the integrals
      ALLOCATE (t_2c_op_pot(nimg), t_2c_op_RI(nimg))
      ALLOCATE (t_3c_int(1, nimg))
      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, t_3c_int, do_kpoints=.TRUE.)

      !Make sure the internals have the k-point format
      CALL adapt_ri_data_to_kp(dbcsr_template, ri_data, qs_env)

      !For each atom i, we calculate the inverse RI metric (P^0 | Q^0)^-1 without external bumping yet
      !Also store the off-diagonal integrals of the RI metric in case of forces, bumped from the left
      DO iatom = 1, natom
         CALL get_ext_2c_int(ri_data%t_2c_inv(1, iatom), t_2c_op_RI, iatom, iatom, 1, ri_data, qs_env, &
                             do_inverse=.TRUE.)
         !for the forces:
         !off-diagonl RI metric bumped from the left
         CALL get_ext_2c_int(ri_data%t_2c_int(1, iatom), t_2c_op_RI, iatom, iatom, 1, ri_data, &
                             qs_env, off_diagonal=.TRUE.)
         CALL apply_bump(ri_data%t_2c_int(1, iatom), iatom, ri_data, qs_env, from_left=.TRUE., from_right=.FALSE.)

         !RI metric with bumped off-diagonal blocks (but not inverted), depumed from left and right
         CALL get_ext_2c_int(ri_data%t_2c_pot(1, iatom), t_2c_op_RI, iatom, iatom, 1, ri_data, qs_env, &
                             do_inverse=.TRUE., skip_inverse=.TRUE.)
         CALL apply_bump(ri_data%t_2c_pot(1, iatom), iatom, ri_data, qs_env, from_left=.TRUE., &
                         from_right=.TRUE., debump=.TRUE.)

      END DO

      DO i_img = 1, nimg
         CALL dbcsr_release(t_2c_op_RI(i_img))
      END DO

      ALLOCATE (ri_data%kp_mat_2c_pot(1, nimg))
      DO i_img = 1, nimg
         CALL dbcsr_create(ri_data%kp_mat_2c_pot(1, i_img), template=t_2c_op_pot(i_img))
         CALL dbcsr_copy(ri_data%kp_mat_2c_pot(1, i_img), t_2c_op_pot(i_img))
         CALL dbcsr_release(t_2c_op_pot(i_img))
      END DO

      !Pre-contract all 3c integrals with the bumped inverse RI metric (P^0|Q^0)^-1,
      !and store in ri_data%t_3c_int_ctr_1
      CALL precontract_3c_ints(t_3c_int, ri_data, qs_env)

      !reorder the 3c integrals such that empty images are bunched up together
      CALL reorder_3c_ints(ri_data%t_3c_int_ctr_1(1, :), ri_data)

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_pre_scf_kp

! **************************************************************************************************
!> \brief clean-up the KP specific data from ri_data
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE cleanup_kp(ri_data)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: i, j

      IF (ALLOCATED(ri_data%kp_cost)) DEALLOCATE (ri_data%kp_cost)
      IF (ALLOCATED(ri_data%idx_to_img)) DEALLOCATE (ri_data%idx_to_img)
      IF (ALLOCATED(ri_data%img_to_idx)) DEALLOCATE (ri_data%img_to_idx)
      IF (ALLOCATED(ri_data%present_images)) DEALLOCATE (ri_data%present_images)
      IF (ALLOCATED(ri_data%img_to_RI_cell)) DEALLOCATE (ri_data%img_to_RI_cell)
      IF (ALLOCATED(ri_data%RI_cell_to_img)) DEALLOCATE (ri_data%RI_cell_to_img)

      IF (ALLOCATED(ri_data%kp_mat_2c_pot)) THEN
         DO j = 1, SIZE(ri_data%kp_mat_2c_pot, 2)
            DO i = 1, SIZE(ri_data%kp_mat_2c_pot, 1)
               CALL dbcsr_release(ri_data%kp_mat_2c_pot(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%kp_mat_2c_pot)
      END IF

      IF (ALLOCATED(ri_data%kp_t_3c_int)) THEN
         DO i = 1, SIZE(ri_data%kp_t_3c_int)
            CALL dbt_destroy(ri_data%kp_t_3c_int(i))
         END DO
         DEALLOCATE (ri_data%kp_t_3c_int)
      END IF

      IF (ALLOCATED(ri_data%t_2c_inv)) THEN
         DO j = 1, SIZE(ri_data%t_2c_inv, 2)
            DO i = 1, SIZE(ri_data%t_2c_inv, 1)
               CALL dbt_destroy(ri_data%t_2c_inv(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%t_2c_inv)
      END IF

      IF (ALLOCATED(ri_data%t_2c_int)) THEN
         DO j = 1, SIZE(ri_data%t_2c_int, 2)
            DO i = 1, SIZE(ri_data%t_2c_int, 1)
               CALL dbt_destroy(ri_data%t_2c_int(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%t_2c_int)
      END IF

      IF (ALLOCATED(ri_data%t_2c_pot)) THEN
         DO j = 1, SIZE(ri_data%t_2c_pot, 2)
            DO i = 1, SIZE(ri_data%t_2c_pot, 1)
               CALL dbt_destroy(ri_data%t_2c_pot(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%t_2c_pot)
      END IF

      IF (ALLOCATED(ri_data%t_3c_int_ctr_1)) THEN
         DO j = 1, SIZE(ri_data%t_3c_int_ctr_1, 2)
            DO i = 1, SIZE(ri_data%t_3c_int_ctr_1, 1)
               CALL dbt_destroy(ri_data%t_3c_int_ctr_1(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%t_3c_int_ctr_1)
      END IF

      IF (ALLOCATED(ri_data%t_3c_int_ctr_2)) THEN
         DO j = 1, SIZE(ri_data%t_3c_int_ctr_2, 2)
            DO i = 1, SIZE(ri_data%t_3c_int_ctr_2, 1)
               CALL dbt_destroy(ri_data%t_3c_int_ctr_2(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%t_3c_int_ctr_2)
      END IF

      IF (ALLOCATED(ri_data%rho_ao_t)) THEN
         DO j = 1, SIZE(ri_data%rho_ao_t, 2)
            DO i = 1, SIZE(ri_data%rho_ao_t, 1)
               CALL dbt_destroy(ri_data%rho_ao_t(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%rho_ao_t)
      END IF

      IF (ALLOCATED(ri_data%ks_t)) THEN
         DO j = 1, SIZE(ri_data%ks_t, 2)
            DO i = 1, SIZE(ri_data%ks_t, 1)
               CALL dbt_destroy(ri_data%ks_t(i, j))
            END DO
         END DO
         DEALLOCATE (ri_data%ks_t)
      END IF

   END SUBROUTINE cleanup_kp

! **************************************************************************************************
!> \brief Update the KS matrices for each real-space image
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_kp(qs_env, ri_data, ks_matrix, ehfx, rho_ao, &
                                  geometry_did_change, nspins, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_kp'

      INTEGER :: b_img, batch_size, group_size, handle, handle2, i_batch, i_img, i_spin, iatom, &
         iblk, igroup, jatom, mb_img, n_batch_nze, natom, ngroups, nimg, nimg_nze
      INTEGER(int_8)                                     :: nflop, nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_ranges_at, batch_ranges_nze, &
                                                            idx_to_at_AO
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: iapc_pairs
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: sparsity_pattern
      LOGICAL                                            :: use_delta_p
      REAL(dp)                                           :: etmp, fac, occ, pfac, pref, t1, t2, t3, &
                                                            t4
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(dbcsr_type)                                   :: ks_desymm, rho_desymm, tmp
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_2c_pot
      TYPE(dbcsr_type), POINTER                          :: dbcsr_template
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: ks_t_split, t_2c_ao_tmp, t_2c_work, &
                                                            t_3c_int, t_3c_work_2, t_3c_work_3
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: ks_t, ks_t_sub, t_3c_apc, t_3c_apc_sub
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(section_vals_type), POINTER                   :: hfx_section

      NULLIFY (para_env, para_env_sub, blacs_env_sub, hfx_section, dbcsr_template)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, natom=natom)

      IF (nspins == 1) THEN
         fac = 0.5_dp*hf_fraction
      ELSE
         fac = 1.0_dp*hf_fraction
      END IF

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_kp(ks_matrix(1, 1)%matrix, ri_data, qs_env)
      END IF
      nimg = ri_data%nimg
      nimg_nze = ri_data%nimg_nze

      !We need to calculate the KS matrix for each periodic cell with index b: F_mu^0,nu^b
      !F_mu^0,nu^b = -0.5 sum_a,c P_sigma^0,lambda^c (mu^0, sigma^a| P^0) V_P^0,Q^b (Q^b| nu^b lambda^a+c)
      !with V_P^0,Q^b = (P^0|R^0)^-1 * (R^0|S^b) * (S^b|Q^b)^-1

      !We use a local RI basis set for each atom in the system, which inlcudes RI basis elements for
      !each neighboring atom standing within the KIND radius (decay of Gaussian with smallest exponent)

      !We also limit the number of periodic images we consider accorrding to the HFX potentail in the
      !RI basis, because if V_P^0,Q^b is zero everywhere, then image b can be ignored (RI basis less diffuse)

      !We manage to calculate each KS matrix doing a double loop on iamges, and a double loop on atoms
      !First, we pre-contract and store P_sigma^0,lambda^c (mu^0, sigma^a| P^0) (P^0|R^0)^-1 into T_mu^0,lambda^a+c,P^0
      !Then, we loop over b_img, iatom, jatom to get (R^0|S^b)
      !Finally, we do an additional loop over a+c images where we do (R^0|S^b) (S^b|Q^b)^-1 (Q^b| nu^b lambda^a+c)
      !and the final contraction with T_mu^0,lambda^a+c,P^0

      !Note that the 3-center integrals are pre-contracted with the RI metric, and that the same tensor can be used
      !(mu^0, sigma^a| P^0) (P^0|R^0)  <===> (S^b|Q^b)^-1 (Q^b| nu^b lambda^a+c) by relabelling the images

      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF%RI")
      CALL section_vals_val_get(hfx_section, "KP_USE_DELTA_P", l_val=use_delta_p)

      !By default, build the density tensor based on the difference of this SCF P and that of the prev. SCF
      pfac = -1.0_dp
      IF (.NOT. use_delta_p) pfac = 0.0_dp
      CALL get_pmat_images(ri_data%rho_ao_t, rho_ao, pfac, ri_data, qs_env)

      ALLOCATE (ks_t(nspins, nimg))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ri_data%ks_t(1, 1), ks_t(i_spin, i_img))
         END DO
      END DO

      ALLOCATE (idx_to_at_AO(SIZE(ri_data%bsizes_AO_split)))
      CALL get_idx_to_atom(idx_to_at_AO, ri_data%bsizes_AO_split, ri_data%bsizes_AO)

      !First we calculate and store T^1_mu^0,lambda^a+c,P = P_mu^0,lambda^c * (mu_0 sigma^a | P^0) (P^0|R^0)^-1
      !To avoid doing nimg**2 tiny contractions that do not scale well with a large number of CPUs,
      !we instead do a single loop over the a+c image index. For each a+c, we get a list of allowed
      !combination of a,c indices. Then we build TAS tensors P_mu^0,lambda^c with all concerned c's
      !and (mu^0 sigma^a | P^0)*(P^0|R^0)^-1 with all a's. Then we perform a single contraction with larger tensors,
      !were the sum over a,c is automatically taken care of
      ALLOCATE (t_3c_apc(nspins, nimg))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ri_data%t_3c_int_ctr_2(1, 1), t_3c_apc(i_spin, i_img))
         END DO
      END DO
      CALL contract_pmat_3c(t_3c_apc, ri_data%rho_ao_t, ri_data, qs_env)

      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF%RI")
      CALL section_vals_val_get(hfx_section, "KP_NGROUPS", i_val=ngroups)
      CALL section_vals_val_get(hfx_section, "KP_STACK_SIZE", i_val=batch_size)
      ri_data%kp_stack_size = batch_size

      IF (MOD(para_env%num_pe, ngroups) .NE. 0) THEN
         CPWARN("KP_NGROUPS must be an integer divisor of the total number of MPI ranks. It was set to 1.")
         ngroups = 1
         CALL section_vals_val_set(hfx_section, "KP_NGROUPS", i_val=ngroups)
      END IF
      IF ((MOD(ngroups, natom) .NE. 0) .AND. (MOD(natom, ngroups) .NE. 0) .AND. geometry_did_change) THEN
         IF (ngroups > 1) &
            CPWARN("Better load balancing is reached if NGROUPS is a multiple/divisor of the number of atoms")
      END IF
      group_size = para_env%num_pe/ngroups
      igroup = para_env%mepos/group_size

      ALLOCATE (para_env_sub)
      CALL para_env_sub%from_split(para_env, igroup)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      ! The sparsity pattern of each iatom, jatom pair, on each b_img, and on which subgroup
      ALLOCATE (sparsity_pattern(natom, natom, nimg))
      CALL get_sparsity_pattern(sparsity_pattern, ri_data, qs_env)
      CALL get_sub_dist(sparsity_pattern, ngroups, ri_data)

      !Get all the required tensors in the subgroups
      ALLOCATE (mat_2c_pot(nimg), ks_t_sub(nspins, nimg), t_2c_ao_tmp(1), ks_t_split(2), t_2c_work(3))
      CALL get_subgroup_2c_tensors(mat_2c_pot, t_2c_work, t_2c_ao_tmp, ks_t_split, ks_t_sub, &
                                   group_size, ngroups, para_env, para_env_sub, ri_data)

      ALLOCATE (t_3c_int(nimg), t_3c_apc_sub(nspins, nimg), t_3c_work_2(3), t_3c_work_3(3))
      CALL get_subgroup_3c_tensors(t_3c_int, t_3c_work_2, t_3c_work_3, t_3c_apc, t_3c_apc_sub, &
                                   group_size, ngroups, para_env, para_env_sub, ri_data)

      !We go atom by atom, therefore there is an automatic batching along that direction
      !Also, because we stack the 3c tensors nimg times, we naturally do some batching there too
      ALLOCATE (batch_ranges_at(natom + 1))
      batch_ranges_at(natom + 1) = SIZE(ri_data%bsizes_AO_split) + 1
      iatom = 0
      DO iblk = 1, SIZE(ri_data%bsizes_AO_split)
         IF (idx_to_at_AO(iblk) == iatom + 1) THEN
            iatom = iatom + 1
            batch_ranges_at(iatom) = iblk
         END IF
      END DO

      n_batch_nze = nimg_nze/batch_size
      IF (MODULO(nimg_nze, batch_size) .NE. 0) n_batch_nze = n_batch_nze + 1
      ALLOCATE (batch_ranges_nze(n_batch_nze + 1))
      DO i_batch = 1, n_batch_nze
         batch_ranges_nze(i_batch) = (i_batch - 1)*batch_size + 1
      END DO
      batch_ranges_nze(n_batch_nze + 1) = nimg_nze + 1

      CALL dbt_batched_contract_init(t_3c_work_3(1), batch_range_2=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_3(2), batch_range_2=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_2(1), batch_range_1=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_2(2), batch_range_1=batch_ranges_at)

      t1 = m_walltime()
      ri_data%kp_cost(:, :, :) = 0.0_dp
      ALLOCATE (iapc_pairs(nimg, 2))
      DO b_img = 1, nimg
         CALL dbt_batched_contract_init(ks_t_split(1))
         CALL dbt_batched_contract_init(ks_t_split(2))
         DO jatom = 1, natom
            DO iatom = 1, natom
               IF (.NOT. sparsity_pattern(iatom, jatom, b_img) == igroup) CYCLE
               pref = 1.0_dp
               IF (iatom == jatom .AND. b_img == 1) pref = 0.5_dp

               !measure the cost of the given i, j, b configuration
               t3 = m_walltime()

               !Get the proper HFX potential 2c integrals (R_i^0|S_j^b)
               CALL timeset(routineN//"_2c", handle2)
               CALL get_ext_2c_int(t_2c_work(1), mat_2c_pot, iatom, jatom, b_img, ri_data, qs_env, &
                                   blacs_env_ext=blacs_env_sub, para_env_ext=para_env_sub, &
                                   dbcsr_template=dbcsr_template)
               CALL dbt_copy(t_2c_work(1), t_2c_work(2), move_data=.TRUE.) !move to split blocks
               CALL dbt_filter(t_2c_work(2), ri_data%filter_eps)
               CALL timestop(handle2)

               CALL dbt_batched_contract_init(t_2c_work(2))
               CALL get_iapc_pairs(iapc_pairs, b_img, ri_data, qs_env)
               CALL timeset(routineN//"_3c", handle2)

               !Stack the (S^b|Q^b)^-1 * (Q^b| nu^b lambda^a+c) integrals over a+c and multiply by (R_i^0|S_j^b)
               DO i_batch = 1, n_batch_nze
                  CALL fill_3c_stack(t_3c_work_3(3), t_3c_int, iapc_pairs(:, 1), 3, ri_data, &
                                     filter_at=jatom, filter_dim=2, idx_to_at=idx_to_at_AO, &
                                     img_bounds=[batch_ranges_nze(i_batch), batch_ranges_nze(i_batch + 1)])
                  CALL dbt_copy(t_3c_work_3(3), t_3c_work_3(1), move_data=.TRUE.)

                  CALL dbt_contract(1.0_dp, t_2c_work(2), t_3c_work_3(1), &
                                    0.0_dp, t_3c_work_3(2), map_1=[1], map_2=[2, 3], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2, 3], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  CALL dbt_copy(t_3c_work_3(2), t_3c_work_2(2), order=[2, 1, 3], move_data=.TRUE.)
                  CALL dbt_copy(t_3c_work_3(3), t_3c_work_3(1))

                  !Stack the P_sigma^a,lambda^a+c * (mu^0 sigma^a | P^0)*(P^0|R^0)^-1 integrals over a+c and contract
                  !to get the final block of the KS matrix
                  DO i_spin = 1, nspins
                     CALL fill_3c_stack(t_3c_work_2(3), t_3c_apc_sub(i_spin, :), iapc_pairs(:, 2), 3, &
                                        ri_data, filter_at=iatom, filter_dim=1, idx_to_at=idx_to_at_AO, &
                                        img_bounds=[batch_ranges_nze(i_batch), batch_ranges_nze(i_batch + 1)])
                     CALL get_tensor_occupancy(t_3c_work_2(3), nze, occ)
                     IF (nze == 0) CYCLE
                     CALL dbt_copy(t_3c_work_2(3), t_3c_work_2(1), move_data=.TRUE.)
                     CALL dbt_contract(-pref*fac, t_3c_work_2(1), t_3c_work_2(2), &
                                       1.0_dp, ks_t_split(i_spin), map_1=[1], map_2=[2], &
                                       contract_1=[2, 3], notcontract_1=[1], &
                                       contract_2=[2, 3], notcontract_2=[1], &
                                       filter_eps=ri_data%filter_eps, &
                                       move_data=i_spin == nspins, flop=nflop)
                     ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  END DO
               END DO !i_batch
               CALL timestop(handle2)
               CALL dbt_batched_contract_finalize(t_2c_work(2))

               t4 = m_walltime()
               ri_data%kp_cost(iatom, jatom, b_img) = t4 - t3
            END DO !iatom
         END DO !jatom
         CALL dbt_batched_contract_finalize(ks_t_split(1))
         CALL dbt_batched_contract_finalize(ks_t_split(2))

         DO i_spin = 1, nspins
            CALL dbt_copy(ks_t_split(i_spin), t_2c_ao_tmp(1), move_data=.TRUE.)
            CALL dbt_copy(t_2c_ao_tmp(1), ks_t_sub(i_spin, b_img), summation=.TRUE.)
         END DO
      END DO !b_img
      CALL dbt_batched_contract_finalize(t_3c_work_3(1))
      CALL dbt_batched_contract_finalize(t_3c_work_3(2))
      CALL dbt_batched_contract_finalize(t_3c_work_2(1))
      CALL dbt_batched_contract_finalize(t_3c_work_2(2))
      CALL para_env%sync()
      CALL para_env%sum(ri_data%dbcsr_nflop)
      CALL para_env%sum(ri_data%kp_cost)
      t2 = m_walltime()
      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      !transfer KS tensor from subgroup to main group
      CALL gather_ks_matrix(ks_t, ks_t_sub, group_size, sparsity_pattern, para_env, ri_data)

      !Keep the 3c integrals on the subgroups to avoid communication at next SCF step
      DO i_img = 1, nimg
         CALL dbt_copy(t_3c_int(i_img), ri_data%kp_t_3c_int(i_img), move_data=.TRUE.)
      END DO

      !clean-up subgroup tensors
      CALL dbt_destroy(t_2c_ao_tmp(1))
      CALL dbt_destroy(ks_t_split(1))
      CALL dbt_destroy(ks_t_split(2))
      CALL dbt_destroy(t_2c_work(1))
      CALL dbt_destroy(t_2c_work(2))
      CALL dbt_destroy(t_3c_work_2(1))
      CALL dbt_destroy(t_3c_work_2(2))
      CALL dbt_destroy(t_3c_work_2(3))
      CALL dbt_destroy(t_3c_work_3(1))
      CALL dbt_destroy(t_3c_work_3(2))
      CALL dbt_destroy(t_3c_work_3(3))
      DO i_img = 1, nimg
         CALL dbt_destroy(t_3c_int(i_img))
         CALL dbcsr_release(mat_2c_pot(i_img))
         DO i_spin = 1, nspins
            CALL dbt_destroy(t_3c_apc_sub(i_spin, i_img))
            CALL dbt_destroy(ks_t_sub(i_spin, i_img))
         END DO
      END DO
      IF (ASSOCIATED(dbcsr_template)) THEN
         CALL dbcsr_release(dbcsr_template)
         DEALLOCATE (dbcsr_template)
      END IF

      !End of subgroup parallelization
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL para_env_sub%free()
      DEALLOCATE (para_env_sub)

      !Currently, rho_ao_t holds the density difference (wrt to pref SCF step).
      !ks_t also hold that diff, while only having half the blocks => need to add to prev ks_t and symmetrize
      !We need the full thing for the energy, on the next SCF step
      CALL get_pmat_images(ri_data%rho_ao_t, rho_ao, 0.0_dp, ri_data, qs_env)
      DO i_spin = 1, nspins
         DO b_img = 1, nimg
            CALL dbt_copy(ks_t(i_spin, b_img), ri_data%ks_t(i_spin, b_img), summation=.TRUE.)

            !desymmetrize
            mb_img = get_opp_index(b_img, qs_env)
            IF (mb_img > 0 .AND. mb_img .LE. nimg) THEN
               CALL dbt_copy(ks_t(i_spin, mb_img), ri_data%ks_t(i_spin, b_img), order=[2, 1], summation=.TRUE.)
            END IF
         END DO
      END DO
      DO b_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_destroy(ks_t(i_spin, b_img))
         END DO
      END DO

      !calculate the energy
      CALL dbt_create(ri_data%ks_t(1, 1), t_2c_ao_tmp(1))
      CALL dbcsr_create(tmp, template=ks_matrix(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(ks_desymm, template=ks_matrix(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(rho_desymm, template=ks_matrix(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      ehfx = 0.0_dp
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_filter(ri_data%ks_t(i_spin, i_img), ri_data%filter_eps)
            CALL dbt_copy(ri_data%ks_t(i_spin, i_img), t_2c_ao_tmp(1))
            CALL dbt_copy_tensor_to_matrix(t_2c_ao_tmp(1), ks_desymm)
            CALL dbt_copy_tensor_to_matrix(t_2c_ao_tmp(1), tmp)
            CALL dbcsr_add(ks_matrix(i_spin, i_img)%matrix, tmp, 1.0_dp, 1.0_dp)

            CALL dbt_copy(ri_data%rho_ao_t(i_spin, i_img), t_2c_ao_tmp(1))
            CALL dbt_copy_tensor_to_matrix(t_2c_ao_tmp(1), rho_desymm)

            CALL dbcsr_dot(ks_desymm, rho_desymm, etmp)
            ehfx = ehfx + 0.5_dp*etmp

            IF (.NOT. use_delta_p) CALL dbt_clear(ri_data%ks_t(i_spin, i_img))
         END DO
      END DO
      CALL dbcsr_release(rho_desymm)
      CALL dbcsr_release(ks_desymm)
      CALL dbcsr_release(tmp)
      CALL dbt_destroy(t_2c_ao_tmp(1))

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_update_ks_kp

! **************************************************************************************************
!> \brief Update the K-points RI-HFX forces
!> \param qs_env ...
!> \param ri_data ...
!> \param nspins ...
!> \param hf_fraction ...
!> \param rho_ao ...
!> \param use_virial ...
!> \note Because this routine uses stored quantities calculated in the energy calculation, they should
!>       always be called by pairs, and with the same input densities
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_forces_kp(qs_env, ri_data, nspins, hf_fraction, rho_ao, use_virial)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_forces_kp'

      INTEGER :: b_img, batch_size, group_size, handle, handle2, i_batch, i_img, i_loop, i_spin, &
         i_xyz, iatom, iblk, igroup, j_xyz, jatom, k_xyz, n_batch, natom, ngroups, nimg, nimg_nze
      INTEGER(int_8)                                     :: nflop, nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, batch_ranges_at, &
                                                            batch_ranges_nze, dist1, dist2, &
                                                            i_images, idx_to_at_AO, idx_to_at_RI, &
                                                            kind_of
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: iapc_pairs
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: force_pattern, sparsity_pattern
      INTEGER, DIMENSION(2, 1)                           :: bounds_iat, bounds_jat
      LOGICAL                                            :: use_virial_prv
      REAL(dp)                                           :: fac, occ, pref, t1, t2
      REAL(dp), DIMENSION(3, 3)                          :: work_virial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_2c_pot
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:, :)     :: mat_der_pot, mat_der_pot_sub
      TYPE(dbcsr_type), POINTER                          :: dbcsr_template
      TYPE(dbt_type)                                     :: t_2c_R, t_2c_R_split
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_2c_bint, t_2c_binv, t_2c_der_pot, &
                                                            t_2c_inv, t_2c_metric, t_2c_work, &
                                                            t_3c_der_stack, t_3c_work_2, &
                                                            t_3c_work_3
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :) :: rho_ao_t, rho_ao_t_sub, t_2c_der_metric, &
         t_2c_der_metric_sub, t_3c_apc, t_3c_apc_sub, t_3c_der_AO, t_3c_der_AO_sub, t_3c_der_RI, &
         t_3c_der_RI_sub
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (para_env, para_env_sub, hfx_section, blacs_env_sub, dbcsr_template, force, atomic_kind_set, &
               virial, particle_set, cell)

      CALL timeset(routineN, handle)

      use_virial_prv = .FALSE.
      IF (PRESENT(use_virial)) use_virial_prv = use_virial

      IF (nspins == 1) THEN
         fac = 0.5_dp*hf_fraction
      ELSE
         fac = 1.0_dp*hf_fraction
      END IF

      CALL get_qs_env(qs_env, natom=natom, para_env=para_env, force=force, cell=cell, virial=virial, &
                      atomic_kind_set=atomic_kind_set, particle_set=particle_set)
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ALLOCATE (idx_to_at_AO(SIZE(ri_data%bsizes_AO_split)))
      CALL get_idx_to_atom(idx_to_at_AO, ri_data%bsizes_AO_split, ri_data%bsizes_AO)

      ALLOCATE (idx_to_at_RI(SIZE(ri_data%bsizes_RI_split)))
      CALL get_idx_to_atom(idx_to_at_RI, ri_data%bsizes_RI_split, ri_data%bsizes_RI)

      nimg = ri_data%nimg
      ALLOCATE (t_3c_der_RI(nimg, 3), t_3c_der_AO(nimg, 3), mat_der_pot(nimg, 3), t_2c_der_metric(natom, 3))

      !We assume that the integrals are available from the SCF
      !pre-calculate the derivs. 3c tensors as (P^0| sigma^a mu^0), with t_3c_der_AO holding deriv wrt mu^0
      CALL precalc_derivatives(t_3c_der_RI, t_3c_der_AO, mat_der_pot, t_2c_der_metric, ri_data, qs_env)

      !Calculate the density matrix at each image
      ALLOCATE (rho_ao_t(nspins, nimg))
      CALL create_2c_tensor(rho_ao_t(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)
      IF (nspins == 2) CALL dbt_create(rho_ao_t(1, 1), rho_ao_t(2, 1))
      DO i_img = 2, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(rho_ao_t(1, 1), rho_ao_t(i_spin, i_img))
         END DO
      END DO
      CALL get_pmat_images(rho_ao_t, rho_ao, 0.0_dp, ri_data, qs_env)

      !Contract integrals with the density matrix
      ALLOCATE (t_3c_apc(nspins, nimg))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(ri_data%t_3c_int_ctr_2(1, 1), t_3c_apc(i_spin, i_img))
         END DO
      END DO
      CALL contract_pmat_3c(t_3c_apc, rho_ao_t, ri_data, qs_env)

      !Setup the subgroups
      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF%RI")
      CALL section_vals_val_get(hfx_section, "KP_NGROUPS", i_val=ngroups)
      group_size = para_env%num_pe/ngroups
      igroup = para_env%mepos/group_size

      ALLOCATE (para_env_sub)
      CALL para_env_sub%from_split(para_env, igroup)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      !Get the ususal sparsity pattern
      ALLOCATE (sparsity_pattern(natom, natom, nimg))
      CALL get_sparsity_pattern(sparsity_pattern, ri_data, qs_env)
      CALL get_sub_dist(sparsity_pattern, ngroups, ri_data)

      !Get the 2-center quantities in the subgroups (note: main group derivs are deleted wihtin)
      ALLOCATE (t_2c_inv(natom), mat_2c_pot(nimg), rho_ao_t_sub(nspins, nimg), t_2c_work(5), &
                t_2c_der_metric_sub(natom, 3), mat_der_pot_sub(nimg, 3), t_2c_bint(natom), &
                t_2c_metric(natom), t_2c_binv(natom))
      CALL get_subgroup_2c_derivs(t_2c_inv, t_2c_bint, t_2c_metric, mat_2c_pot, t_2c_work, rho_ao_t, &
                                  rho_ao_t_sub, t_2c_der_metric, t_2c_der_metric_sub, mat_der_pot, &
                                  mat_der_pot_sub, group_size, ngroups, para_env, para_env_sub, ri_data)
      CALL dbt_create(t_2c_work(1), t_2c_R) !nRI x nRI
      CALL dbt_create(t_2c_work(5), t_2c_R_split) !nRI x nRI with split blocks

      ALLOCATE (t_2c_der_pot(3))
      DO i_xyz = 1, 3
         CALL dbt_create(t_2c_R, t_2c_der_pot(i_xyz))
      END DO

      !Get the 3-center quantities in the subgroups. The integrals and t_3c_apc already there
      !TODO: also distribute the t_3c_apc tensors
      ALLOCATE (t_3c_work_2(3), t_3c_work_3(4), t_3c_der_stack(6), t_3c_der_AO_sub(nimg, 3), &
                t_3c_der_RI_sub(nimg, 3), t_3c_apc_sub(nspins, nimg))
      CALL get_subgroup_3c_derivs(t_3c_work_2, t_3c_work_3, t_3c_der_AO, t_3c_der_AO_sub, &
                                  t_3c_der_RI, t_3c_der_RI_sub, t_3c_apc, t_3c_apc_sub, t_3c_der_stack, &
                                  group_size, ngroups, para_env, para_env_sub, ri_data)

      !Set up batched contraction (go atom by atom)
      ALLOCATE (batch_ranges_at(natom + 1))
      batch_ranges_at(natom + 1) = SIZE(ri_data%bsizes_AO_split) + 1
      iatom = 0
      DO iblk = 1, SIZE(ri_data%bsizes_AO_split)
         IF (idx_to_at_AO(iblk) == iatom + 1) THEN
            iatom = iatom + 1
            batch_ranges_at(iatom) = iblk
         END IF
      END DO

      CALL dbt_batched_contract_init(t_3c_work_3(1), batch_range_2=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_3(2), batch_range_2=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_3(3), batch_range_2=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_2(1), batch_range_1=batch_ranges_at)
      CALL dbt_batched_contract_init(t_3c_work_2(2), batch_range_1=batch_ranges_at)

      !Preparing for the stacking of 3c tensors
      nimg_nze = ri_data%nimg_nze
      batch_size = ri_data%kp_stack_size
      n_batch = nimg_nze/batch_size
      IF (MODULO(nimg_nze, batch_size) .NE. 0) n_batch = n_batch + 1
      ALLOCATE (batch_ranges_nze(n_batch + 1))
      DO i_batch = 1, n_batch
         batch_ranges_nze(i_batch) = (i_batch - 1)*batch_size + 1
      END DO
      batch_ranges_nze(n_batch + 1) = nimg_nze + 1

      !Applying the external bump to ((P|Q)_D + B*(P|Q)_OD*B)^-1 from left and right
      !And keep the bump on LHS only version as well, with B*M^-1 = (M^-1*B)^T
      DO iatom = 1, natom
         CALL dbt_create(t_2c_inv(iatom), t_2c_binv(iatom))
         CALL dbt_copy(t_2c_inv(iatom), t_2c_binv(iatom))
         CALL apply_bump(t_2c_binv(iatom), iatom, ri_data, qs_env, from_left=.TRUE., from_right=.FALSE.)
         CALL apply_bump(t_2c_inv(iatom), iatom, ri_data, qs_env, from_left=.TRUE., from_right=.TRUE.)
      END DO

      t1 = m_walltime()
      work_virial = 0.0_dp
      ALLOCATE (iapc_pairs(nimg, 2), i_images(nimg))
      ALLOCATE (force_pattern(natom, natom, nimg))
      force_pattern(:, :, :) = -1
      !We proceed with 2 loops: one over the sparsity pattern from the SCF, one over the rest
      !We use the SCF cost model for the first loop, while we calculate the cost of the upcoming loop
      DO i_loop = 1, 2
         DO b_img = 1, nimg
            DO jatom = 1, natom
               DO iatom = 1, natom

                  pref = -0.5_dp*fac
                  IF (i_loop == 1 .AND. (.NOT. sparsity_pattern(iatom, jatom, b_img) == igroup)) CYCLE
                  IF (i_loop == 2 .AND. (.NOT. force_pattern(iatom, jatom, b_img) == igroup)) CYCLE

                  !Get the proper HFX potential 2c integrals (R_i^0|S_j^b), times (S_j^b|Q_j^b)^-1
                  CALL timeset(routineN//"_2c_1", handle2)
                  CALL get_ext_2c_int(t_2c_work(1), mat_2c_pot, iatom, jatom, b_img, ri_data, qs_env, &
                                      blacs_env_ext=blacs_env_sub, para_env_ext=para_env_sub, &
                                      dbcsr_template=dbcsr_template)
                  CALL dbt_contract(1.0_dp, t_2c_work(1), t_2c_inv(jatom), &
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  CALL dbt_copy(t_2c_work(2), t_2c_work(5), move_data=.TRUE.) !move to split blocks
                  CALL dbt_filter(t_2c_work(5), ri_data%filter_eps)
                  CALL timestop(handle2)

                  CALL timeset(routineN//"_3c", handle2)
                  bounds_iat(:, 1) = [SUM(ri_data%bsizes_AO(1:iatom - 1)) + 1, SUM(ri_data%bsizes_AO(1:iatom))]
                  bounds_jat(:, 1) = [SUM(ri_data%bsizes_AO(1:jatom - 1)) + 1, SUM(ri_data%bsizes_AO(1:jatom))]
                  CALL dbt_clear(t_2c_R_split)

                  DO i_spin = 1, nspins
                     CALL dbt_batched_contract_init(rho_ao_t_sub(i_spin, b_img))
                  END DO

                  CALL get_iapc_pairs(iapc_pairs, b_img, ri_data, qs_env, i_images) !i = a+c-b
                  DO i_batch = 1, n_batch

                     !Stack the 3c derivatives to take the trace later on
                     DO i_xyz = 1, 3
                        CALL dbt_clear(t_3c_der_stack(i_xyz))
                        CALL fill_3c_stack(t_3c_der_stack(i_xyz), t_3c_der_RI_sub(:, i_xyz), &
                                           iapc_pairs(:, 1), 3, ri_data, filter_at=jatom, &
                                           filter_dim=2, idx_to_at=idx_to_at_AO, &
                                           img_bounds=[batch_ranges_nze(i_batch), batch_ranges_nze(i_batch + 1)])

                        CALL dbt_clear(t_3c_der_stack(3 + i_xyz))
                        CALL fill_3c_stack(t_3c_der_stack(3 + i_xyz), t_3c_der_AO_sub(:, i_xyz), &
                                           iapc_pairs(:, 1), 3, ri_data, filter_at=jatom, &
                                           filter_dim=2, idx_to_at=idx_to_at_AO, &
                                           img_bounds=[batch_ranges_nze(i_batch), batch_ranges_nze(i_batch + 1)])
                     END DO

                     DO i_spin = 1, nspins
                        !stack the t_3c_apc tensors
                        CALL dbt_clear(t_3c_work_2(3))
                        CALL fill_3c_stack(t_3c_work_2(3), t_3c_apc_sub(i_spin, :), iapc_pairs(:, 2), 3, &
                                           ri_data, filter_at=iatom, filter_dim=1, idx_to_at=idx_to_at_AO, &
                                           img_bounds=[batch_ranges_nze(i_batch), batch_ranges_nze(i_batch + 1)])
                        CALL get_tensor_occupancy(t_3c_work_2(3), nze, occ)
                        IF (nze == 0) CYCLE
                        CALL dbt_copy(t_3c_work_2(3), t_3c_work_2(1), move_data=.TRUE.)

                        !Contract with the second density matrix: P_mu^0,nu^b * t_3c_apc,
                        !where t_3c_apc = P_sigma^a,lambda^a+c (mu^0 P^0 sigma^a) *(P^0|R^0)^-1 (stacked along a+c)
                        CALL dbt_contract(1.0_dp, rho_ao_t_sub(i_spin, b_img), t_3c_work_2(1), &
                                          0.0_dp, t_3c_work_2(2), map_1=[1], map_2=[2, 3], &
                                          contract_1=[1], notcontract_1=[2], &
                                          contract_2=[1], notcontract_2=[2, 3], &
                                          bounds_1=bounds_iat, bounds_2=bounds_jat, &
                                          filter_eps=ri_data%filter_eps, flop=nflop)
                        ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                        CALL get_tensor_occupancy(t_3c_work_2(2), nze, occ)
                        IF (nze == 0) CYCLE

                        !Contract with V_PQ so that we can take the trace with (Q^b|nu^b lmabda^a+c)^(x)
                        CALL dbt_copy(t_3c_work_2(2), t_3c_work_3(1), order=[2, 1, 3], move_data=.TRUE.)
                        CALL dbt_batched_contract_init(t_2c_work(5))
                        CALL dbt_contract(1.0_dp, t_2c_work(5), t_3c_work_3(1), &
                                          0.0_dp, t_3c_work_3(2), map_1=[1], map_2=[2, 3], &
                                          contract_1=[1], notcontract_1=[2], &
                                          contract_2=[1], notcontract_2=[2, 3], &
                                          filter_eps=ri_data%filter_eps, flop=nflop)
                        ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                        CALL dbt_batched_contract_finalize(t_2c_work(5))

                        !Contract with the 3c derivatives to get the force/virial
                        CALL dbt_copy(t_3c_work_3(2), t_3c_work_3(4), move_data=.TRUE.)
                        IF (use_virial_prv) THEN
                           CALL get_force_from_3c_trace(force, t_3c_work_3(4), t_3c_der_stack(1:3), &
                                                        t_3c_der_stack(4:6), atom_of_kind, kind_of, &
                                                        idx_to_at_RI, idx_to_at_AO, i_images, &
                                                        batch_ranges_nze(i_batch), 2.0_dp*pref, &
                                                        ri_data, qs_env, work_virial, cell, particle_set)
                        ELSE
                           CALL get_force_from_3c_trace(force, t_3c_work_3(4), t_3c_der_stack(1:3), &
                                                        t_3c_der_stack(4:6), atom_of_kind, kind_of, &
                                                        idx_to_at_RI, idx_to_at_AO, i_images, &
                                                        batch_ranges_nze(i_batch), 2.0_dp*pref, &
                                                        ri_data, qs_env)
                        END IF
                        CALL dbt_clear(t_3c_work_3(4))

                        !Contract with the 3-center integrals in order to have a matrix R_PQ such that
                        !we can take the trace sum_PQ R_PQ (P^0|Q^b)^(x)
                        IF (i_loop == 2) CYCLE

                        !Stack the 3c integrals
                        CALL fill_3c_stack(t_3c_work_3(4), ri_data%kp_t_3c_int, iapc_pairs(:, 1), 3, ri_data, &
                                           filter_at=jatom, filter_dim=2, idx_to_at=idx_to_at_AO, &
                                           img_bounds=[batch_ranges_nze(i_batch), batch_ranges_nze(i_batch + 1)])
                        CALL dbt_copy(t_3c_work_3(4), t_3c_work_3(3), move_data=.TRUE.)

                        CALL dbt_batched_contract_init(t_2c_R_split)
                        CALL dbt_contract(1.0_dp, t_3c_work_3(1), t_3c_work_3(3), &
                                          1.0_dp, t_2c_R_split, map_1=[1], map_2=[2], &
                                          contract_1=[2, 3], notcontract_1=[1], &
                                          contract_2=[2, 3], notcontract_2=[1], &
                                          filter_eps=ri_data%filter_eps, flop=nflop)
                        ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                        CALL dbt_batched_contract_finalize(t_2c_R_split)
                        CALL dbt_copy(t_3c_work_3(4), t_3c_work_3(1))
                     END DO
                  END DO
                  DO i_spin = 1, nspins
                     CALL dbt_batched_contract_finalize(rho_ao_t_sub(i_spin, b_img))
                  END DO
                  CALL timestop(handle2)

                  IF (i_loop == 2) CYCLE
                  pref = 2.0_dp*pref
                  IF (iatom == jatom .AND. b_img == 1) pref = 0.5_dp*pref

                  CALL timeset(routineN//"_2c_2", handle2)
                  !Note that the derivatives are in atomic block format (not split)
                  CALL dbt_copy(t_2c_R_split, t_2c_R, move_data=.TRUE.)

                  CALL get_ext_2c_int(t_2c_work(1), mat_2c_pot, iatom, jatom, b_img, ri_data, qs_env, &
                                      blacs_env_ext=blacs_env_sub, para_env_ext=para_env_sub, &
                                      dbcsr_template=dbcsr_template)

                  !We have to calculate: S^-1(iat) * R_PQ * S^-1(jat)    to trace with HFX pot der
                  !                      + R_PQ * S^-1(jat) * pot^T      to trace with S^(x) (iat)
                  !                      + pot^T * S^-1(iat) *R_PQ       to trace with S^(x) (jat)

                  !Because 3c tensors are all precontracted with the inverse RI metric,
                  !t_2c_R is currently implicitely multiplied by S^-1(iat) from the left
                  !and S^-1(jat) from the right, directly in the proper format for the trace
                  !with the HFX potential derivative

                  !Trace with HFX pot deriv, that we need to build first
                  DO i_xyz = 1, 3
                     CALL get_ext_2c_int(t_2c_der_pot(i_xyz), mat_der_pot_sub(:, i_xyz), iatom, jatom, &
                                         b_img, ri_data, qs_env, blacs_env_ext=blacs_env_sub, &
                                         para_env_ext=para_env_sub, dbcsr_template=dbcsr_template)
                  END DO

                  IF (use_virial_prv) THEN
                     CALL get_2c_der_force(force, t_2c_R, t_2c_der_pot, atom_of_kind, kind_of, &
                                           b_img, pref, ri_data, qs_env, work_virial, cell, particle_set)
                  ELSE
                     CALL get_2c_der_force(force, t_2c_R, t_2c_der_pot, atom_of_kind, kind_of, &
                                           b_img, pref, ri_data, qs_env)
                  END IF

                  DO i_xyz = 1, 3
                     CALL dbt_clear(t_2c_der_pot(i_xyz))
                  END DO

                  !R_PQ * S^-1(jat) * pot^T  (=A)
                  CALL dbt_contract(1.0_dp, t_2c_metric(iatom), t_2c_R, & !get rid of implicit S^-1(iat)
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  CALL dbt_contract(1.0_dp, t_2c_work(2), t_2c_work(1), &
                                    0.0_dp, t_2c_work(3), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[2], notcontract_2=[1], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  !With the RI bump function, things get more complex. M = (S|P)_D + B*(S|P)_OD*B
                  !Calculate M^-1*B*A + A*B*M^-1 to contract with B^x. A is in t_2c_work(3)
                  CALL dbt_contract(1.0_dp, t_2c_work(3), t_2c_binv(iatom), &
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_contract(1.0_dp, t_2c_binv(iatom), t_2c_work(3), & !use transpose of B*M^-1 = M^-1*B
                                    0.0_dp, t_2c_work(4), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_copy(t_2c_work(2), t_2c_work(4), summation=.TRUE.)
                  CALL get_2c_bump_forces(force, t_2c_work(4), iatom, atom_of_kind, kind_of, pref, &
                                          ri_data, qs_env, work_virial)

                  !Calculate -M^-1*B*A*B*M^-1 to contracte with diagonal RI metric deriv. t_2c_work(2) holds A*B*M^-1
                  CALL dbt_contract(1.0_dp, t_2c_binv(iatom), t_2c_work(2), &
                                    0.0_dp, t_2c_work(4), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  IF (use_virial_prv) THEN
                     CALL get_2c_der_force(force, t_2c_work(4), t_2c_der_metric_sub(iatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, work_virial, cell, particle_set, &
                                           diag=.TRUE., offdiag=.FALSE.)
                  ELSE
                     CALL get_2c_der_force(force, t_2c_work(4), t_2c_der_metric_sub(iatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, diag=.TRUE., offdiag=.FALSE.)
                  END IF

                  !Calculate -B*M^-1*B*A*B*M^-1*B to contract with off-diagonal RI metric derivs
                  CALL dbt_copy(t_2c_work(4), t_2c_work(2))
                  CALL apply_bump(t_2c_work(2), iatom, ri_data, qs_env, from_left=.TRUE., from_right=.TRUE.)

                  IF (use_virial_prv) THEN
                     CALL get_2c_der_force(force, t_2c_work(2), t_2c_der_metric_sub(iatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, work_virial, cell, particle_set, &
                                           diag=.FALSE., offdiag=.TRUE.)
                  ELSE
                     CALL get_2c_der_force(force, t_2c_work(2), t_2c_der_metric_sub(iatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, diag=.FALSE., offdiag=.TRUE.)
                  END IF

                  !Calculate -O*B*M^-1*B*A*B*M^-1 - M^-1*B*A*B*M^-1*B*O, where O is off-diagonal integrals
                  !t_2c_work(4) holds M^-1*B*A*B*M^-1, and exploit transpose of B*O (stored in t_2c_bint)
                  CALL dbt_contract(1.0_dp, t_2c_work(4), t_2c_bint(iatom), &
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_contract(1.0_dp, t_2c_bint(iatom), t_2c_work(4), &
                                    1.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL get_2c_bump_forces(force, t_2c_work(2), iatom, atom_of_kind, kind_of, -pref, &
                                          ri_data, qs_env, work_virial)

                  ! pot^T * S^-1(iat) * R_PQ (=A)
                  CALL dbt_contract(1.0_dp, t_2c_work(1), t_2c_R, &
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_contract(1.0_dp, t_2c_work(2), t_2c_metric(jatom), & !get rid of implicit S^-1(jat)
                                    0.0_dp, t_2c_work(3), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  !Do the same shenanigans with the S^(x) (jatom)
                  !Calculate M^-1*B*A + A*B*M^-1 to contract with B^x. A is in t_2c_work(3)
                  CALL dbt_contract(1.0_dp, t_2c_work(3), t_2c_binv(jatom), &
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_contract(1.0_dp, t_2c_binv(jatom), t_2c_work(3), & !use transpose of B*M^-1 = M^-1*B
                                    0.0_dp, t_2c_work(4), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_copy(t_2c_work(2), t_2c_work(4), summation=.TRUE.)
                  CALL get_2c_bump_forces(force, t_2c_work(4), jatom, atom_of_kind, kind_of, pref, &
                                          ri_data, qs_env, work_virial)

                  !Calculate -M^-1*B*A*B*M^-1 to contracte with diagonal RI metric deriv. t_2c_work(2) holds A*B*M^-1
                  CALL dbt_contract(1.0_dp, t_2c_binv(jatom), t_2c_work(2), &
                                    0.0_dp, t_2c_work(4), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  IF (use_virial_prv) THEN
                     CALL get_2c_der_force(force, t_2c_work(4), t_2c_der_metric_sub(jatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, work_virial, cell, particle_set, &
                                           diag=.TRUE., offdiag=.FALSE.)
                  ELSE
                     CALL get_2c_der_force(force, t_2c_work(4), t_2c_der_metric_sub(jatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, diag=.TRUE., offdiag=.FALSE.)
                  END IF

                  !Calculate -B*M^-1*B*A*B*M^-1*B to contract with off-diagonal RI metric derivs
                  CALL dbt_copy(t_2c_work(4), t_2c_work(2))
                  CALL apply_bump(t_2c_work(2), jatom, ri_data, qs_env, from_left=.TRUE., from_right=.TRUE.)

                  IF (use_virial_prv) THEN
                     CALL get_2c_der_force(force, t_2c_work(2), t_2c_der_metric_sub(jatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, work_virial, cell, particle_set, &
                                           diag=.FALSE., offdiag=.TRUE.)
                  ELSE
                     CALL get_2c_der_force(force, t_2c_work(2), t_2c_der_metric_sub(jatom, :), atom_of_kind, &
                                           kind_of, 1, -pref, ri_data, qs_env, diag=.FALSE., offdiag=.TRUE.)
                  END IF

                  !Calculate -O*B*M^-1*B*A*B*M^-1 - M^-1*B*A*B*M^-1*B*O, where O is off-diagonal integrals
                  !t_2c_work(4) holds M^-1*B*A*B*M^-1, and exploit transpose of B*O (stored in t_2c_bint)
                  CALL dbt_contract(1.0_dp, t_2c_work(4), t_2c_bint(jatom), &
                                    0.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL dbt_contract(1.0_dp, t_2c_bint(jatom), t_2c_work(4), &
                                    1.0_dp, t_2c_work(2), map_1=[1], map_2=[2], &
                                    contract_1=[1], notcontract_1=[2], &
                                    contract_2=[1], notcontract_2=[2], &
                                    filter_eps=ri_data%filter_eps, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

                  CALL get_2c_bump_forces(force, t_2c_work(2), jatom, atom_of_kind, kind_of, -pref, &
                                          ri_data, qs_env, work_virial)

                  CALL timestop(handle2)
               END DO !iatom
            END DO !jatom
         END DO !b_img

         IF (i_loop == 1) THEN
            CALL update_pattern_to_forces(force_pattern, sparsity_pattern, ngroups, ri_data, qs_env)
         END IF
      END DO !i_loop

      CALL dbt_batched_contract_finalize(t_3c_work_3(1))
      CALL dbt_batched_contract_finalize(t_3c_work_3(2))
      CALL dbt_batched_contract_finalize(t_3c_work_3(3))
      CALL dbt_batched_contract_finalize(t_3c_work_2(1))
      CALL dbt_batched_contract_finalize(t_3c_work_2(2))

      IF (use_virial_prv) THEN
         DO k_xyz = 1, 3
            DO j_xyz = 1, 3
               DO i_xyz = 1, 3
                  virial%pv_fock_4c(i_xyz, j_xyz) = virial%pv_fock_4c(i_xyz, j_xyz) &
                                                    + work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
               END DO
            END DO
         END DO
      END IF

      !End of subgroup parallelization
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL para_env_sub%free()
      DEALLOCATE (para_env_sub)

      CALL para_env%sync()
      t2 = m_walltime()
      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      !clean-up
      IF (ASSOCIATED(dbcsr_template)) THEN
         CALL dbcsr_release(dbcsr_template)
         DEALLOCATE (dbcsr_template)
      END IF
      CALL dbt_destroy(t_2c_R)
      CALL dbt_destroy(t_2c_R_split)
      CALL dbt_destroy(t_2c_work(1))
      CALL dbt_destroy(t_2c_work(2))
      CALL dbt_destroy(t_2c_work(3))
      CALL dbt_destroy(t_2c_work(4))
      CALL dbt_destroy(t_2c_work(5))
      CALL dbt_destroy(t_3c_work_2(1))
      CALL dbt_destroy(t_3c_work_2(2))
      CALL dbt_destroy(t_3c_work_2(3))
      CALL dbt_destroy(t_3c_work_3(1))
      CALL dbt_destroy(t_3c_work_3(2))
      CALL dbt_destroy(t_3c_work_3(3))
      CALL dbt_destroy(t_3c_work_3(4))
      CALL dbt_destroy(t_3c_der_stack(1))
      CALL dbt_destroy(t_3c_der_stack(2))
      CALL dbt_destroy(t_3c_der_stack(3))
      CALL dbt_destroy(t_3c_der_stack(4))
      CALL dbt_destroy(t_3c_der_stack(5))
      CALL dbt_destroy(t_3c_der_stack(6))
      DO i_xyz = 1, 3
         CALL dbt_destroy(t_2c_der_pot(i_xyz))
      END DO
      DO iatom = 1, natom
         CALL dbt_destroy(t_2c_inv(iatom))
         CALL dbt_destroy(t_2c_binv(iatom))
         CALL dbt_destroy(t_2c_bint(iatom))
         CALL dbt_destroy(t_2c_metric(iatom))
         DO i_xyz = 1, 3
            CALL dbt_destroy(t_2c_der_metric_sub(iatom, i_xyz))
         END DO
      END DO
      DO i_img = 1, nimg
         CALL dbcsr_release(mat_2c_pot(i_img))
         DO i_spin = 1, nspins
            CALL dbt_destroy(rho_ao_t_sub(i_spin, i_img))
            CALL dbt_destroy(t_3c_apc_sub(i_spin, i_img))
         END DO
      END DO
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbt_destroy(t_3c_der_RI_sub(i_img, i_xyz))
            CALL dbt_destroy(t_3c_der_AO_sub(i_img, i_xyz))
            CALL dbcsr_release(mat_der_pot_sub(i_img, i_xyz))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_update_forces_kp

! **************************************************************************************************
!> \brief A routine the applies the RI bump matrix from the left and/or the right, given an input
!>        matrix and the central RI atom. We assume atomic block sizes
!> \param t_2c_inout ...
!> \param atom_i ...
!> \param ri_data ...
!> \param qs_env ...
!> \param from_left ...
!> \param from_right ...
!> \param debump ...
! **************************************************************************************************
   SUBROUTINE apply_bump(t_2c_inout, atom_i, ri_data, qs_env, from_left, from_right, debump)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_inout
      INTEGER, INTENT(IN)                                :: atom_i
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: from_left, from_right, debump

      INTEGER                                            :: i_img, i_RI, iatom, ind(2), j_img, j_RI, &
                                                            jatom, natom, nblks(2), nimg, nkind
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, my_debump, my_left, my_right
      REAL(dp)                                           :: bval, r0, r1, ri(3), rj(3), rref(3), &
                                                            scoord(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, particle_set, kpoints, index_to_cell, cell_to_index, cell)

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, cell=cell, &
                      kpoints=kpoints, particle_set=particle_set)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      my_debump = .FALSE.
      IF (PRESENT(debump)) my_debump = debump

      my_left = .FALSE.
      IF (PRESENT(from_left)) my_left = from_left

      my_right = .FALSE.
      IF (PRESENT(from_right)) my_right = from_right
      CPASSERT(my_left .OR. my_right)

      CALL dbt_get_info(t_2c_inout, nblks_total=nblks)
      CPASSERT(nblks(1) == ri_data%ncell_RI*natom)
      CPASSERT(nblks(2) == ri_data%ncell_RI*natom)

      nimg = ri_data%nimg

      !Loop over the RI cells and atoms, and apply bump accordingly
      r1 = ri_data%kp_RI_range
      r0 = ri_data%kp_bump_rad
      rref = pbc(particle_set(atom_i)%r, cell)

!$OMP PARALLEL DEFAULT(NONE) SHARED(t_2c_inout,natom,ri_data,cell,particle_set,index_to_cell,my_left, &
!$OMP                               my_right,r0,r1,rref,my_debump) &
!$OMP PRIVATE(iter,ind,blk,found,i_RI,i_img,iatom,j_RI,j_img,jatom,scoord,ri,rj,bval)
      CALL dbt_iterator_start(iter, t_2c_inout)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_block(t_2c_inout, ind, blk, found)
         IF (.NOT. found) CYCLE

         i_RI = (ind(1) - 1)/natom + 1
         i_img = ri_data%RI_cell_to_img(i_RI)
         iatom = ind(1) - (i_RI - 1)*natom

         CALL real_to_scaled(scoord, pbc(particle_set(iatom)%r, cell), cell)
         CALL scaled_to_real(ri, scoord(:) + index_to_cell(:, i_img), cell)

         j_RI = (ind(2) - 1)/natom + 1
         j_img = ri_data%RI_cell_to_img(j_RI)
         jatom = ind(2) - (j_RI - 1)*natom

         CALL real_to_scaled(scoord, pbc(particle_set(jatom)%r, cell), cell)
         CALL scaled_to_real(rj, scoord(:) + index_to_cell(:, j_img), cell)

         IF (.NOT. my_debump) THEN
            IF (my_left) blk(:, :) = blk(:, :)*bump(NORM2(ri - rref), r0, r1)
            IF (my_right) blk(:, :) = blk(:, :)*bump(NORM2(rj - rref), r0, r1)
         ELSE
            !Note: by construction, the bump function is never quite zero, as its range is the same
            !      as that of the extended RI basis (but we are safe)
            bval = bump(NORM2(ri - rref), r0, r1)
            IF (my_left .AND. bval > EPSILON(1.0_dp)) blk(:, :) = blk(:, :)/bval
            bval = bump(NORM2(rj - rref), r0, r1)
            IF (my_right .AND. bval > EPSILON(1.0_dp)) blk(:, :) = blk(:, :)/bval
         END IF

         CALL dbt_put_block(t_2c_inout, ind, SHAPE(blk), blk)

         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      CALL dbt_filter(t_2c_inout, ri_data%filter_eps)

   END SUBROUTINE apply_bump

! **************************************************************************************************
!> \brief A routine that calculates the forces due to the derivative of the bump function
!> \param force ...
!> \param t_2c_in ...
!> \param atom_i ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param pref ...
!> \param ri_data ...
!> \param qs_env ...
!> \param work_virial ...
! **************************************************************************************************
   SUBROUTINE get_2c_bump_forces(force, t_2c_in, atom_i, atom_of_kind, kind_of, pref, ri_data, &
                                 qs_env, work_virial)
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_in
      INTEGER, INTENT(IN)                                :: atom_i
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of
      REAL(dp), INTENT(IN)                               :: pref
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT)           :: work_virial

      INTEGER :: i, i_img, i_RI, i_xyz, iat_of_kind, iatom, ikind, ind(2), j_img, j_RI, j_xyz, &
         jat_of_kind, jatom, jkind, natom, nblks(2), nimg, nkind
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(dp)                                           :: new_force, r0, r1, ri(3), rj(3), &
                                                            rref(3), scoord(3), x
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, particle_set, kpoints, index_to_cell, cell_to_index, cell)

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, cell=cell, &
                      kpoints=kpoints, particle_set=particle_set)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      CALL dbt_get_info(t_2c_in, nblks_total=nblks)
      CPASSERT(nblks(1) == ri_data%ncell_RI*natom)
      CPASSERT(nblks(2) == ri_data%ncell_RI*natom)

      nimg = ri_data%nimg

      !Loop over the RI cells and atoms, and apply bump accordingly
      r1 = ri_data%kp_RI_range
      r0 = ri_data%kp_bump_rad
      rref = pbc(particle_set(atom_i)%r, cell)

      iat_of_kind = atom_of_kind(atom_i)
      ikind = kind_of(atom_i)

!$OMP PARALLEL DEFAULT(NONE) SHARED(t_2c_in,natom,ri_data,cell,particle_set,index_to_cell,pref, &
!$OMP force,r0,r1,rref,atom_of_kind,kind_of,iat_of_kind,ikind,work_virial) &
!$OMP PRIVATE(iter,ind,blk,found,i_RI,i_img,iatom,j_RI,j_img,jatom,scoord,ri,rj,jkind,jat_of_kind, &
!$OMP         new_force,i_xyz,i,x,j_xyz)
      CALL dbt_iterator_start(iter, t_2c_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         IF (ind(1) .NE. ind(2)) CYCLE !bump matrix is diagonal

         CALL dbt_get_block(t_2c_in, ind, blk, found)
         IF (.NOT. found) CYCLE

         !bump is a function of x = SQRT((R - Rref)^2). We refer to R as jatom, and Rref as atom_i
         j_RI = (ind(2) - 1)/natom + 1
         j_img = ri_data%RI_cell_to_img(j_RI)
         jatom = ind(2) - (j_RI - 1)*natom
         jat_of_kind = atom_of_kind(jatom)
         jkind = kind_of(jatom)

         CALL real_to_scaled(scoord, pbc(particle_set(jatom)%r, cell), cell)
         CALL scaled_to_real(rj, scoord(:) + index_to_cell(:, j_img), cell)
         x = NORM2(rj - rref)
         IF (x < r0 .OR. x > r1) CYCLE

         new_force = 0.0_dp
         DO i = 1, SIZE(blk, 1)
            new_force = new_force + blk(i, i)
         END DO
         new_force = pref*new_force*dbump(x, r0, r1)

         !x = SQRT((R - Rref)^2), so we multiply by dx/dR and dx/dRref
         DO i_xyz = 1, 3
            !Force acting on second atom
!$OMP ATOMIC
            force(jkind)%fock_4c(i_xyz, jat_of_kind) = force(jkind)%fock_4c(i_xyz, jat_of_kind) + &
                                                       new_force*(rj(i_xyz) - rref(i_xyz))/x

            !virial acting on second atom
            CALL real_to_scaled(scoord, rj, cell)
            DO j_xyz = 1, 3
!$OMP ATOMIC
               work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) &
                                           + new_force*scoord(j_xyz)*(rj(i_xyz) - rref(i_xyz))/x
            END DO

            !Force acting on reference atom, defining the RI basis
!$OMP ATOMIC
            force(ikind)%fock_4c(i_xyz, iat_of_kind) = force(ikind)%fock_4c(i_xyz, iat_of_kind) - &
                                                       new_force*(rj(i_xyz) - rref(i_xyz))/x

            !virial of ref atom
            CALL real_to_scaled(scoord, rref, cell)
            DO j_xyz = 1, 3
!$OMP ATOMIC
               work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) &
                                           - new_force*scoord(j_xyz)*(rj(i_xyz) - rref(i_xyz))/x
            END DO
         END DO !i_xyz

         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

   END SUBROUTINE get_2c_bump_forces

! **************************************************************************************************
!> \brief The bumb function as defined by Juerg
!> \param x ...
!> \param r0 ...
!> \param r1 ...
!> \return ...
! **************************************************************************************************
   FUNCTION bump(x, r0, r1) RESULT(b)
      REAL(dp), INTENT(IN)                               :: x, r0, r1
      REAL(dp)                                           :: b

      REAL(dp)                                           :: r

      !Head-Gordon
      !b = 1.0_dp/(1.0_dp+EXP((r1-r0)/(r1-x)-(r1-r0)/(x-r0)))
      !Juerg
      r = (x - r0)/(r1 - r0)
      b = -6.0_dp*r**5 + 15.0_dp*r**4 - 10.0_dp*r**3 + 1.0_dp
      IF (x .GE. r1) b = 0.0_dp
      IF (x .LE. r0) b = 1.0_dp

   END FUNCTION bump

! **************************************************************************************************
!> \brief The derivative of the bump function
!> \param x ...
!> \param r0 ...
!> \param r1 ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbump(x, r0, r1) RESULT(b)
      REAL(dp), INTENT(IN)                               :: x, r0, r1
      REAL(dp)                                           :: b

      REAL(dp)                                           :: r

      r = (x - r0)/(r1 - r0)
      b = (-30.0_dp*r**4 + 60.0_dp*r**3 - 30.0_dp*r**2)/(r1 - r0)
      IF (x .GE. r1) b = 0.0_dp
      IF (x .LE. r0) b = 0.0_dp

   END FUNCTION dbump

! **************************************************************************************************
!> \brief return the cell index a+c corresponding to given cell index i and b, with i = a+c-b
!> \param i_index ...
!> \param b_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_apc_index_from_ib(i_index, b_index, qs_env) RESULT(apc_index)
      INTEGER, INTENT(IN)                                :: i_index, b_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: apc_index

      INTEGER, DIMENSION(3)                              :: cell_apc
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      !i = a+c-b => a+c = i+b
      cell_apc(:) = index_to_cell(:, i_index) + index_to_cell(:, b_index)

      IF (ANY([cell_apc(1), cell_apc(2), cell_apc(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([cell_apc(1), cell_apc(2), cell_apc(3)] > UBOUND(cell_to_index))) THEN

         apc_index = 0
      ELSE
         apc_index = cell_to_index(cell_apc(1), cell_apc(2), cell_apc(3))
      END IF

   END FUNCTION get_apc_index_from_ib

! **************************************************************************************************
!> \brief return the cell index i corresponding to the summ of cell_a and cell_c
!> \param a_index ...
!> \param c_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_apc_index(a_index, c_index, qs_env) RESULT(i_index)
      INTEGER, INTENT(IN)                                :: a_index, c_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: i_index

      INTEGER, DIMENSION(3)                              :: cell_i
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      cell_i(:) = index_to_cell(:, a_index) + index_to_cell(:, c_index)

      IF (ANY([cell_i(1), cell_i(2), cell_i(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([cell_i(1), cell_i(2), cell_i(3)] > UBOUND(cell_to_index))) THEN

         i_index = 0
      ELSE
         i_index = cell_to_index(cell_i(1), cell_i(2), cell_i(3))
      END IF

   END FUNCTION get_apc_index

! **************************************************************************************************
!> \brief return the cell index i corresponding to the summ of cell_a + cell_c - cell_b
!> \param apc_index ...
!> \param b_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_i_index(apc_index, b_index, qs_env) RESULT(i_index)
      INTEGER, INTENT(IN)                                :: apc_index, b_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: i_index

      INTEGER, DIMENSION(3)                              :: cell_i
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      cell_i(:) = index_to_cell(:, apc_index) - index_to_cell(:, b_index)

      IF (ANY([cell_i(1), cell_i(2), cell_i(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([cell_i(1), cell_i(2), cell_i(3)] > UBOUND(cell_to_index))) THEN

         i_index = 0
      ELSE
         i_index = cell_to_index(cell_i(1), cell_i(2), cell_i(3))
      END IF

   END FUNCTION get_i_index

! **************************************************************************************************
!> \brief A routine that returns all allowed a,c pairs such that a+c images corresponds to the value
!>        of the apc_index input. Takes into account that image a corresponds to 3c integrals, which
!>        are ordered in their own way
!> \param ac_pairs ...
!> \param apc_index ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_ac_pairs(ac_pairs, apc_index, ri_data, qs_env)
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: ac_pairs
      INTEGER, INTENT(IN)                                :: apc_index
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: a_index, actual_img, c_index, nimg

      nimg = SIZE(ac_pairs, 1)

      ac_pairs(:, :) = 0
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(ac_pairs,nimg,ri_data,qs_env,apc_index) &
!$OMP PRIVATE(a_index,actual_img,c_index)
      DO a_index = 1, nimg
         actual_img = ri_data%idx_to_img(a_index)
         !c = a+c - a
         c_index = get_i_index(apc_index, actual_img, qs_env)
         ac_pairs(a_index, 1) = a_index
         ac_pairs(a_index, 2) = c_index
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE get_ac_pairs

! **************************************************************************************************
!> \brief A routine that returns all allowed i,a+c pairs such that, for the given value of b, we have
!>        i = a+c-b. Takes into account that image i corrsponds to the 3c ints, which are ordered in
!>        their own way
!> \param iapc_pairs ...
!> \param b_index ...
!> \param ri_data ...
!> \param qs_env ...
!> \param actual_i_img ...
! **************************************************************************************************
   SUBROUTINE get_iapc_pairs(iapc_pairs, b_index, ri_data, qs_env, actual_i_img)
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: iapc_pairs
      INTEGER, INTENT(IN)                                :: b_index
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL     :: actual_i_img

      INTEGER                                            :: actual_img, apc_index, i_index, nimg

      nimg = SIZE(iapc_pairs, 1)
      IF (PRESENT(actual_i_img)) actual_i_img(:) = 0

      iapc_pairs(:, :) = 0
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(iapc_pairs,nimg,ri_data,qs_env,b_index,actual_i_img) &
!$OMP PRIVATE(i_index,actual_img,apc_index)
      DO i_index = 1, nimg
         actual_img = ri_data%idx_to_img(i_index)
         apc_index = get_apc_index_from_ib(actual_img, b_index, qs_env)
         IF (apc_index == 0) CYCLE
         iapc_pairs(i_index, 1) = i_index
         iapc_pairs(i_index, 2) = apc_index
         IF (PRESENT(actual_i_img)) actual_i_img(i_index) = actual_img
      END DO

   END SUBROUTINE get_iapc_pairs

! **************************************************************************************************
!> \brief A function that, given a cell index a, returun the index corresponding to -a, and zero if
!>        if out of bounds
!> \param a_index ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_opp_index(a_index, qs_env) RESULT(opp_index)
      INTEGER, INTENT(IN)                                :: a_index
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: opp_index

      INTEGER, DIMENSION(3)                              :: opp_cell
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (kpoints, cell_to_index, index_to_cell)

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      opp_cell(:) = -index_to_cell(:, a_index)

      IF (ANY([opp_cell(1), opp_cell(2), opp_cell(3)] < LBOUND(cell_to_index)) .OR. &
          ANY([opp_cell(1), opp_cell(2), opp_cell(3)] > UBOUND(cell_to_index))) THEN

         opp_index = 0
      ELSE
         opp_index = cell_to_index(opp_cell(1), opp_cell(2), opp_cell(3))
      END IF

   END FUNCTION get_opp_index

! **************************************************************************************************
!> \brief A routine that returns the actual non-symemtric density matrix for each image, by Fourier
!>        transforming the kpoint density matrix
!> \param rho_ao_t ...
!> \param rho_ao ...
!> \param scale_prev_p ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_pmat_images(rho_ao_t, rho_ao, scale_prev_p, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: rho_ao_t
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: rho_ao
      REAL(dp), INTENT(IN)                               :: scale_prev_p
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: cell_j(3), i_img, i_spin, iatom, icol, &
                                                            irow, j_img, jatom, mi_img, mj_img, &
                                                            nimg, nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(dp)                                           :: fac
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock, pblock_desymm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, rho_desymm
      TYPE(dbt_type)                                     :: tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl, sab_nl_nosym
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      NULLIFY (rho_desymm, kpoints, sab_nl_nosym, scf_env, matrix_ks, dft_control, &
               sab_nl, nl_iterator, cell_to_index, pblock, pblock_desymm)

      CALL get_qs_env(qs_env, kpoints=kpoints, scf_env=scf_env, matrix_ks_kp=matrix_ks, dft_control=dft_control)
      CALL get_kpoint_info(kpoints, sab_nl_nosym=sab_nl_nosym, cell_to_index=cell_to_index, sab_nl=sab_nl)

      IF (dft_control%do_admm) THEN
         CALL get_admm_env(qs_env%admm_env, matrix_ks_aux_fit_kp=matrix_ks)
      END IF

      nspins = SIZE(matrix_ks, 1)
      nimg = ri_data%nimg

      ALLOCATE (rho_desymm(nspins, nimg))
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            ALLOCATE (rho_desymm(i_spin, i_img)%matrix)
            CALL dbcsr_create(rho_desymm(i_spin, i_img)%matrix, template=matrix_ks(i_spin, i_img)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(rho_desymm(i_spin, i_img)%matrix, sab_nl_nosym)
         END DO
      END DO
      CALL dbt_create(rho_desymm(1, 1)%matrix, tmp)

      !We transfor the symmtric typed (but not actually symmetric: P_ab^i = P_ba^-i) real-spaced density
      !matrix into proper non-symemtric ones (using the same nl for consistency)
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_j)
         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         fac = 1.0_dp
         IF (iatom == jatom) fac = 0.5_dp
         mj_img = get_opp_index(j_img, qs_env)
         !if no opposite image, then no sum of P^j + P^-j => need full diag
         IF (mj_img == 0) fac = 1.0_dp

         irow = iatom
         icol = jatom
         IF (iatom > jatom) THEN
            !because symmetric nl. Value for atom pair i,j is actually stored in j,i if i > j
            irow = jatom
            icol = iatom
         END IF

         DO i_spin = 1, nspins
            CALL dbcsr_get_block_p(rho_ao(i_spin, j_img)%matrix, irow, icol, pblock, found)
            IF (.NOT. found) CYCLE

            !distribution of symm and non-symm matrix match in that way
            CALL dbcsr_get_block_p(rho_desymm(i_spin, j_img)%matrix, iatom, jatom, pblock_desymm, found)
            IF (.NOT. found) CYCLE

            IF (iatom > jatom) THEN
               pblock_desymm(:, :) = fac*TRANSPOSE(pblock(:, :))
            ELSE
               pblock_desymm(:, :) = fac*pblock(:, :)
            END IF
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_scale(rho_ao_t(i_spin, i_img), scale_prev_p)

            CALL dbt_copy_matrix_to_tensor(rho_desymm(i_spin, i_img)%matrix, tmp)
            CALL dbt_copy(tmp, rho_ao_t(i_spin, i_img), summation=.TRUE., move_data=.TRUE.)

            !symmetrize by addin transpose of opp img
            mi_img = get_opp_index(i_img, qs_env)
            IF (mi_img > 0 .AND. mi_img .LE. nimg) THEN
               CALL dbt_copy_matrix_to_tensor(rho_desymm(i_spin, mi_img)%matrix, tmp)
               CALL dbt_copy(tmp, rho_ao_t(i_spin, i_img), order=[2, 1], summation=.TRUE., move_data=.TRUE.)
            END IF
            CALL dbt_filter(rho_ao_t(i_spin, i_img), ri_data%filter_eps)
         END DO
      END DO

      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbcsr_release(rho_desymm(i_spin, i_img)%matrix)
            DEALLOCATE (rho_desymm(i_spin, i_img)%matrix)
         END DO
      END DO

      CALL dbt_destroy(tmp)
      DEALLOCATE (rho_desymm)

   END SUBROUTINE get_pmat_images

! **************************************************************************************************
!> \brief A routine that, given a cell index b and atom indices ij, returns a 2c tensor with the HFX
!>        potential (P_i^0|Q_j^b), within the extended RI basis
!> \param t_2c_pot ...
!> \param mat_orig ...
!> \param atom_i ...
!> \param atom_j ...
!> \param img_b ...
!> \param ri_data ...
!> \param qs_env ...
!> \param do_inverse ...
!> \param para_env_ext ...
!> \param blacs_env_ext ...
!> \param dbcsr_template ...
!> \param off_diagonal ...
!> \param skip_inverse ...
! **************************************************************************************************
   SUBROUTINE get_ext_2c_int(t_2c_pot, mat_orig, atom_i, atom_j, img_b, ri_data, qs_env, do_inverse, &
                             para_env_ext, blacs_env_ext, dbcsr_template, off_diagonal, skip_inverse)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_pot
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: mat_orig
      INTEGER, INTENT(IN)                                :: atom_i, atom_j, img_b
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_inverse
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_ext
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: dbcsr_template
      LOGICAL, INTENT(IN), OPTIONAL                      :: off_diagonal, skip_inverse

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_ext_2c_int'

      INTEGER :: blk, group, handle, handle2, i_img, i_RI, iatom, iblk, ikind, img_tot, j_img, &
         j_RI, jatom, jblk, jkind, n_dependent, natom, nblks_RI, nimg, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: present_atoms_i, present_atoms_j
      INTEGER, DIMENSION(3)                              :: cell_b, cell_i, cell_j, cell_tot
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_dist_ext, ri_blk_size_ext, &
                                                            row_dist, row_dist_ext
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell, pgrid
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_inverse_prv, found, my_offd, &
                                                            skip_inverse_prv, use_template
      REAL(dp)                                           :: bfac, dij, r0, r1, threshold
      REAL(dp), DIMENSION(3)                             :: ri, rij, rj, rref, scoord
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist, dbcsr_dist_ext
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter
      TYPE(dbcsr_type)                                   :: work, work_tight, work_tight_inv
      TYPE(dbt_type)                                     :: t_2c_tmp
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, nl_2c, nl_iterator, cell, kpoints, cell_to_index, index_to_cell, dist_2d, &
               para_env, pblock, blacs_env, particle_set, col_dist, row_dist, pgrid, &
               col_dist_ext, row_dist_ext)

      CALL timeset(routineN, handle)

      !Idea: run over the neighbor list once for i and once for j, and record in which cell the MIC
      !      atoms are. Then loop over the atoms and only take the pairs the we need

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, cell=cell, &
                      kpoints=kpoints, para_env=para_env, blacs_env=blacs_env, particle_set=particle_set)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      do_inverse_prv = .FALSE.
      IF (PRESENT(do_inverse)) do_inverse_prv = do_inverse
      IF (do_inverse_prv) THEN
         CPASSERT(atom_i == atom_j)
      END IF

      skip_inverse_prv = .FALSE.
      IF (PRESENT(skip_inverse)) skip_inverse_prv = skip_inverse

      my_offd = .FALSE.
      IF (PRESENT(off_diagonal)) my_offd = off_diagonal

      IF (PRESENT(para_env_ext)) para_env => para_env_ext
      IF (PRESENT(blacs_env_ext)) blacs_env => blacs_env_ext

      nimg = SIZE(mat_orig)

      CALL timeset(routineN//"_nl_iter", handle2)

      !create our own dist_2d in the subgroup
      ALLOCATE (dist1(natom), dist2(natom))
      DO iatom = 1, natom
         dist1(iatom) = MOD(iatom, blacs_env%num_pe(1))
         dist2(iatom) = MOD(iatom, blacs_env%num_pe(2))
      END DO
      CALL distribution_2d_create(dist_2d, dist1, dist2, nkind, particle_set, blacs_env_ext=blacs_env)

      ALLOCATE (basis_set_RI(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)

      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%ri_metric, &
                                   "HFX_2c_nl_RI", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)

      ALLOCATE (present_atoms_i(natom, nimg), present_atoms_j(natom, nimg))
      present_atoms_i = 0
      present_atoms_j = 0

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, r=rij, cell=cell_j, &
                                ikind=ikind, jkind=jkind)

         dij = NORM2(rij)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         IF (iatom == atom_i .AND. dij .LE. ri_data%kp_RI_range) present_atoms_i(jatom, j_img) = 1
         IF (iatom == atom_j .AND. dij .LE. ri_data%kp_RI_range) present_atoms_j(jatom, j_img) = 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(nl_2c)
      CALL distribution_2d_release(dist_2d)
      CALL timestop(handle2)

      CALL para_env%sum(present_atoms_i)
      CALL para_env%sum(present_atoms_j)

      !Need to build a work matrix with matching distribution to mat_orig
      !If template is provided, use it. If not, we create it.
      use_template = .FALSE.
      IF (PRESENT(dbcsr_template)) THEN
         IF (ASSOCIATED(dbcsr_template)) use_template = .TRUE.
      END IF

      IF (use_template) THEN
         CALL dbcsr_create(work, template=dbcsr_template)
      ELSE
         CALL dbcsr_get_info(mat_orig(1), distribution=dbcsr_dist)
         CALL dbcsr_distribution_get(dbcsr_dist, row_dist=row_dist, col_dist=col_dist, group=group, pgrid=pgrid)
         ALLOCATE (row_dist_ext(ri_data%ncell_RI*natom), col_dist_ext(ri_data%ncell_RI*natom))
         ALLOCATE (ri_blk_size_ext(ri_data%ncell_RI*natom))
         DO i_RI = 1, ri_data%ncell_RI
            row_dist_ext((i_RI - 1)*natom + 1:i_RI*natom) = row_dist(:)
            col_dist_ext((i_RI - 1)*natom + 1:i_RI*natom) = col_dist(:)
            RI_blk_size_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         END DO

         CALL dbcsr_distribution_new(dbcsr_dist_ext, group=group, pgrid=pgrid, &
                                     row_dist=row_dist_ext, col_dist=col_dist_ext)
         CALL dbcsr_create(work, dist=dbcsr_dist_ext, name="RI_ext", matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=RI_blk_size_ext, col_blk_size=RI_blk_size_ext)
         CALL dbcsr_distribution_release(dbcsr_dist_ext)
         DEALLOCATE (col_dist_ext, row_dist_ext, RI_blk_size_ext)

         IF (PRESENT(dbcsr_template)) THEN
            ALLOCATE (dbcsr_template)
            CALL dbcsr_create(dbcsr_template, template=work)
         END IF
      END IF !use_template

      cell_b(:) = index_to_cell(:, img_b)
      DO i_img = 1, nimg
         i_RI = ri_data%img_to_RI_cell(i_img)
         IF (i_RI == 0) CYCLE
         cell_i(:) = index_to_cell(:, i_img)
         DO j_img = 1, nimg
            j_RI = ri_data%img_to_RI_cell(j_img)
            IF (j_RI == 0) CYCLE
            cell_j(:) = index_to_cell(:, j_img)
            cell_tot = cell_j - cell_i + cell_b

            IF (ANY([cell_tot(1), cell_tot(2), cell_tot(3)] < LBOUND(cell_to_index)) .OR. &
                ANY([cell_tot(1), cell_tot(2), cell_tot(3)] > UBOUND(cell_to_index))) CYCLE
            img_tot = cell_to_index(cell_tot(1), cell_tot(2), cell_tot(3))
            IF (img_tot > nimg .OR. img_tot < 1) CYCLE

            CALL dbcsr_iterator_start(dbcsr_iter, mat_orig(img_tot))
            DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
               CALL dbcsr_iterator_next_block(dbcsr_iter, row=iatom, column=jatom, blk=blk)
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               IF (present_atoms_j(jatom, j_img) == 0) CYCLE
               IF (my_offd .AND. (i_RI - 1)*natom + iatom == (j_RI - 1)*natom + jatom) CYCLE

               CALL dbcsr_get_block_p(mat_orig(img_tot), iatom, jatom, pblock, found)
               IF (.NOT. found) CYCLE

               CALL dbcsr_put_block(work, (i_RI - 1)*natom + iatom, (j_RI - 1)*natom + jatom, pblock)

            END DO
            CALL dbcsr_iterator_stop(dbcsr_iter)

         END DO !j_img
      END DO !i_img
      CALL dbcsr_finalize(work)

      IF (do_inverse_prv) THEN

         r1 = ri_data%kp_RI_range
         r0 = ri_data%kp_bump_rad

         !Because there are a lot of empty rows/cols in work, we need to get rid of them for inversion
         nblks_RI = SUM(present_atoms_i)
         ALLOCATE (col_dist_ext(nblks_RI), row_dist_ext(nblks_RI), RI_blk_size_ext(nblks_RI))
         iblk = 0
         DO i_img = 1, nimg
            i_RI = ri_data%img_to_RI_cell(i_img)
            IF (i_RI == 0) CYCLE
            DO iatom = 1, natom
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               iblk = iblk + 1
               col_dist_ext(iblk) = col_dist(iatom)
               row_dist_ext(iblk) = row_dist(iatom)
               RI_blk_size_ext(iblk) = ri_data%bsizes_RI(iatom)
            END DO
         END DO

         CALL dbcsr_distribution_new(dbcsr_dist_ext, group=group, pgrid=pgrid, &
                                     row_dist=row_dist_ext, col_dist=col_dist_ext)
         CALL dbcsr_create(work_tight, dist=dbcsr_dist_ext, name="RI_ext", matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=RI_blk_size_ext, col_blk_size=RI_blk_size_ext)
         CALL dbcsr_create(work_tight_inv, dist=dbcsr_dist_ext, name="RI_ext", matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=RI_blk_size_ext, col_blk_size=RI_blk_size_ext)
         CALL dbcsr_distribution_release(dbcsr_dist_ext)
         DEALLOCATE (col_dist_ext, row_dist_ext, RI_blk_size_ext)

         !We apply a bump function to the RI metric inverse for smooth RI basis extension:
         ! S^-1 = B * ((P|Q)_D + B*(P|Q)_OD*B)^-1 * B, with D block-diagonal blocks and OD off-diagonal
         rref = pbc(particle_set(atom_i)%r, cell)

         iblk = 0
         DO i_img = 1, nimg
            i_RI = ri_data%img_to_RI_cell(i_img)
            IF (i_RI == 0) CYCLE
            DO iatom = 1, natom
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               iblk = iblk + 1

               CALL real_to_scaled(scoord, pbc(particle_set(iatom)%r, cell), cell)
               CALL scaled_to_real(ri, scoord(:) + index_to_cell(:, i_img), cell)

               jblk = 0
               DO j_img = 1, nimg
                  j_RI = ri_data%img_to_RI_cell(j_img)
                  IF (j_RI == 0) CYCLE
                  DO jatom = 1, natom
                     IF (present_atoms_j(jatom, j_img) == 0) CYCLE
                     jblk = jblk + 1

                     CALL real_to_scaled(scoord, pbc(particle_set(jatom)%r, cell), cell)
                     CALL scaled_to_real(rj, scoord(:) + index_to_cell(:, j_img), cell)

                     CALL dbcsr_get_block_p(work, (i_RI - 1)*natom + iatom, (j_RI - 1)*natom + jatom, pblock, found)
                     IF (.NOT. found) CYCLE

                     bfac = 1.0_dp
                     IF (iblk .NE. jblk) bfac = bump(NORM2(ri - rref), r0, r1)*bump(NORM2(rj - rref), r0, r1)
                     CALL dbcsr_put_block(work_tight, iblk, jblk, bfac*pblock(:, :))
                  END DO
               END DO
            END DO
         END DO
         CALL dbcsr_finalize(work_tight)
         CALL dbcsr_clear(work)

         IF (.NOT. skip_inverse_prv) THEN
            SELECT CASE (ri_data%t2c_method)
            CASE (hfx_ri_do_2c_iter)
               threshold = MAX(ri_data%filter_eps, 1.0e-12_dp)
               CALL invert_hotelling(work_tight_inv, work_tight, threshold=threshold, silent=.FALSE.)
            CASE (hfx_ri_do_2c_cholesky)
               CALL dbcsr_copy(work_tight_inv, work_tight)
               CALL cp_dbcsr_cholesky_decompose(work_tight_inv, para_env=para_env, blacs_env=blacs_env)
               CALL cp_dbcsr_cholesky_invert(work_tight_inv, para_env=para_env, blacs_env=blacs_env, &
                                             upper_to_full=.TRUE.)
            CASE (hfx_ri_do_2c_diag)
               CALL dbcsr_copy(work_tight_inv, work_tight)
               CALL cp_dbcsr_power(work_tight_inv, -1.0_dp, ri_data%eps_eigval, n_dependent, &
                                   para_env, blacs_env, verbose=ri_data%unit_nr_dbcsr > 0)
            END SELECT
         ELSE
            CALL dbcsr_copy(work_tight_inv, work_tight)
         END IF

         !move back data to standard extended RI pattern
         !Note: we apply the external bump to ((P|Q)_D + B*(P|Q)_OD*B)^-1 later, because this matrix
         !      is required for forces
         iblk = 0
         DO i_img = 1, nimg
            i_RI = ri_data%img_to_RI_cell(i_img)
            IF (i_RI == 0) CYCLE
            DO iatom = 1, natom
               IF (present_atoms_i(iatom, i_img) == 0) CYCLE
               iblk = iblk + 1

               jblk = 0
               DO j_img = 1, nimg
                  j_RI = ri_data%img_to_RI_cell(j_img)
                  IF (j_RI == 0) CYCLE
                  DO jatom = 1, natom
                     IF (present_atoms_j(jatom, j_img) == 0) CYCLE
                     jblk = jblk + 1

                     CALL dbcsr_get_block_p(work_tight_inv, iblk, jblk, pblock, found)
                     IF (.NOT. found) CYCLE

                     CALL dbcsr_put_block(work, (i_RI - 1)*natom + iatom, (j_RI - 1)*natom + jatom, pblock)
                  END DO
               END DO
            END DO
         END DO
         CALL dbcsr_finalize(work)

         CALL dbcsr_release(work_tight)
         CALL dbcsr_release(work_tight_inv)
      END IF

      CALL dbt_create(work, t_2c_tmp)
      CALL dbt_copy_matrix_to_tensor(work, t_2c_tmp)
      CALL dbt_copy(t_2c_tmp, t_2c_pot, move_data=.TRUE.)
      CALL dbt_filter(t_2c_pot, ri_data%filter_eps)

      CALL dbt_destroy(t_2c_tmp)
      CALL dbcsr_release(work)

      CALL timestop(handle)

   END SUBROUTINE get_ext_2c_int

! **************************************************************************************************
!> \brief Pre-contract the density matrices with the 3-center integrals:
!>        P_sigma^a,lambda^a+c (mu^0 sigma^a| P^0)
!> \param t_3c_apc ...
!> \param rho_ao_t ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE contract_pmat_3c(t_3c_apc, rho_ao_t, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_apc, rho_ao_t
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'contract_pmat_3c'

      INTEGER                                            :: apc_img, batch_size, handle, i_batch, &
                                                            i_img, i_spin, j_batch, n_batch_img, &
                                                            n_batch_nze, nimg, nimg_nze, nspins
      INTEGER(int_8)                                     :: nflop, nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_ranges_img, batch_ranges_nze, &
                                                            int_indices
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ac_pairs
      REAL(dp)                                           :: occ, t1, t2
      TYPE(dbt_type)                                     :: t_3c_tmp
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: ints_stack, res_stack, rho_stack
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      nimg = ri_data%nimg
      nimg_nze = ri_data%nimg_nze
      nspins = dft_control%nspins

      CALL dbt_create(t_3c_apc(1, 1), t_3c_tmp)

      batch_size = nimg/ri_data%n_mem

      !batching over all images
      n_batch_img = nimg/batch_size
      IF (MODULO(nimg, batch_size) .NE. 0) n_batch_img = n_batch_img + 1
      ALLOCATE (batch_ranges_img(n_batch_img + 1))
      DO i_batch = 1, n_batch_img
         batch_ranges_img(i_batch) = (i_batch - 1)*batch_size + 1
      END DO
      batch_ranges_img(n_batch_img + 1) = nimg + 1

      !batching over images with non-zero 3c integrals
      n_batch_nze = nimg_nze/batch_size
      IF (MODULO(nimg_nze, batch_size) .NE. 0) n_batch_nze = n_batch_nze + 1
      ALLOCATE (batch_ranges_nze(n_batch_nze + 1))
      DO i_batch = 1, n_batch_nze
         batch_ranges_nze(i_batch) = (i_batch - 1)*batch_size + 1
      END DO
      batch_ranges_nze(n_batch_nze + 1) = nimg_nze + 1

      !Create the stack tensors in the approriate distribution
      ALLOCATE (rho_stack(2), ints_stack(2), res_stack(2))
      CALL get_stack_tensors(res_stack, rho_stack, ints_stack, rho_ao_t(1, 1), &
                             ri_data%t_3c_int_ctr_1(1, 1), batch_size, ri_data, qs_env)

      ALLOCATE (ac_pairs(nimg, 2), int_indices(nimg_nze))
      DO i_img = 1, nimg_nze
         int_indices(i_img) = i_img
      END DO

      t1 = m_walltime()
      DO j_batch = 1, n_batch_nze
         !First batch is over the integrals. They are always in the same order, consistent with get_ac_pairs
         CALL fill_3c_stack(ints_stack(1), ri_data%t_3c_int_ctr_1(1, :), int_indices, 3, ri_data, &
                            img_bounds=[batch_ranges_nze(j_batch), batch_ranges_nze(j_batch + 1)])
         CALL dbt_copy(ints_stack(1), ints_stack(2), move_data=.TRUE.)

         DO i_spin = 1, nspins
            DO i_batch = 1, n_batch_img
               !Second batch is over the P matrix. Here we fill the stacked rho tensors col by col
               DO apc_img = batch_ranges_img(i_batch), batch_ranges_img(i_batch + 1) - 1
                  CALL get_ac_pairs(ac_pairs, apc_img, ri_data, qs_env)
                  CALL fill_2c_stack(rho_stack(1), rho_ao_t(i_spin, :), ac_pairs(:, 2), 1, ri_data, &
                                     img_bounds=[batch_ranges_nze(j_batch), batch_ranges_nze(j_batch + 1)], &
                                     shift=apc_img - batch_ranges_img(i_batch) + 1)

               END DO !apc_img
               CALL get_tensor_occupancy(rho_stack(1), nze, occ)
               IF (nze == 0) CYCLE
               CALL dbt_copy(rho_stack(1), rho_stack(2), move_data=.TRUE.)

               !The actual contraction
               CALL dbt_batched_contract_init(rho_stack(2))
               CALL dbt_contract(1.0_dp, ints_stack(2), rho_stack(2), &
                                 0.0_dp, res_stack(2), map_1=[1, 2], map_2=[3], &
                                 contract_1=[3], notcontract_1=[1, 2], &
                                 contract_2=[1], notcontract_2=[2], &
                                 filter_eps=ri_data%filter_eps, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               CALL dbt_batched_contract_finalize(rho_stack(2))
               CALL dbt_copy(res_stack(2), res_stack(1), move_data=.TRUE.)

               DO apc_img = batch_ranges_img(i_batch), batch_ranges_img(i_batch + 1) - 1
                  !Destack the resulting tensor and put it in t_3c_apc with correct apc_img
                  CALL unstack_t_3c_apc(t_3c_tmp, res_stack(1), apc_img - batch_ranges_img(i_batch) + 1)
                  CALL dbt_copy(t_3c_tmp, t_3c_apc(i_spin, apc_img), summation=.TRUE., move_data=.TRUE.)
               END DO

            END DO !i_batch
         END DO !i_spin
      END DO !j_batch
      DEALLOCATE (batch_ranges_img)
      DEALLOCATE (batch_ranges_nze)
      t2 = m_walltime()
      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      CALL dbt_destroy(rho_stack(1))
      CALL dbt_destroy(rho_stack(2))
      CALL dbt_destroy(ints_stack(1))
      CALL dbt_destroy(ints_stack(2))
      CALL dbt_destroy(res_stack(1))
      CALL dbt_destroy(res_stack(2))
      CALL dbt_destroy(t_3c_tmp)

      CALL timestop(handle)

   END SUBROUTINE contract_pmat_3c

! **************************************************************************************************
!> \brief Pre-contract 3-center integrals with the bumped invrse RI metric, for each atom
!> \param t_3c_int ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE precontract_3c_ints(t_3c_int, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_int
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'precontract_3c_ints'

      INTEGER                                            :: batch_size, handle, i_batch, i_img, &
                                                            i_RI, iatom, is, n_batch, natom, &
                                                            nblks, nblks_3c(3), nimg
      INTEGER(int_8)                                     :: nflop
      INTEGER, ALLOCATABLE, DIMENSION(:) :: batch_ranges, bsizes_RI_ext, bsizes_RI_ext_split, &
         bsizes_stack, dist1, dist2, dist3, dist_stack3, idx_to_at_AO, int_indices
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_type)                                     :: t_2c_RI_tmp(2), t_3c_tmp(3)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, natom=natom)

      nimg = ri_data%nimg
      ALLOCATE (int_indices(nimg))
      DO i_img = 1, nimg
         int_indices(i_img) = i_img
      END DO

      ALLOCATE (idx_to_at_AO(SIZE(ri_data%bsizes_AO_split)))
      CALL get_idx_to_atom(idx_to_at_AO, ri_data%bsizes_AO_split, ri_data%bsizes_AO)

      nblks = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*natom))
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         bsizes_RI_ext_split((i_RI - 1)*nblks + 1:i_RI*nblks) = ri_data%bsizes_RI_split(:)
      END DO
      CALL create_2c_tensor(t_2c_RI_tmp(1), dist1, dist2, ri_data%pgrid_2d, &
                            bsizes_RI_ext, bsizes_RI_ext, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)
      CALL create_2c_tensor(t_2c_RI_tmp(2), dist1, dist2, ri_data%pgrid_2d, &
                            bsizes_RI_ext_split, bsizes_RI_ext_split, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      !For more efficiency, we stack multiple images of the 3-center integrals into a single tensor
      batch_size = nimg/ri_data%n_mem
      n_batch = nimg/batch_size
      IF (MODULO(nimg, batch_size) .NE. 0) n_batch = n_batch + 1
      ALLOCATE (batch_ranges(n_batch + 1))
      DO i_batch = 1, n_batch
         batch_ranges(i_batch) = (i_batch - 1)*batch_size + 1
      END DO
      batch_ranges(n_batch + 1) = nimg + 1

      nblks = SIZE(ri_data%bsizes_AO_split)
      ALLOCATE (bsizes_stack(batch_size*nblks))
      DO is = 1, batch_size
         bsizes_stack((is - 1)*nblks + 1:is*nblks) = ri_data%bsizes_AO_split(:)
      END DO

      CALL dbt_get_info(t_3c_int(1, 1), nblks_total=nblks_3c)
      ALLOCATE (dist1(nblks_3c(1)), dist2(nblks_3c(2)), dist3(nblks_3c(3)), dist_stack3(batch_size*nblks_3c(3)))
      CALL dbt_get_info(t_3c_int(1, 1), proc_dist_1=dist1, proc_dist_2=dist2, proc_dist_3=dist3)
      DO is = 1, batch_size
         dist_stack3((is - 1)*nblks_3c(3) + 1:is*nblks_3c(3)) = dist3(:)
      END DO

      CALL dbt_distribution_new(t_dist, ri_data%pgrid, dist1, dist2, dist_stack3)
      CALL dbt_create(t_3c_tmp(1), "ints_stack", t_dist, [1], [2, 3], bsizes_RI_ext_split, &
                      ri_data%bsizes_AO_split, bsizes_stack)
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist3, dist_stack3)

      CALL dbt_create(t_3c_tmp(1), t_3c_tmp(2))
      CALL dbt_create(t_3c_int(1, 1), t_3c_tmp(3))

      DO iatom = 1, natom
         CALL dbt_copy(ri_data%t_2c_inv(1, iatom), t_2c_RI_tmp(1))
         CALL apply_bump(t_2c_RI_tmp(1), iatom, ri_data, qs_env, from_left=.TRUE., from_right=.TRUE.)
         CALL dbt_copy(t_2c_RI_tmp(1), t_2c_RI_tmp(2), move_data=.TRUE.)

         CALL dbt_batched_contract_init(t_2c_RI_tmp(2))
         DO i_batch = 1, n_batch

            CALL fill_3c_stack(t_3c_tmp(1), t_3c_int(1, :), int_indices, 3, ri_data, &
                               img_bounds=[batch_ranges(i_batch), batch_ranges(i_batch + 1)], &
                               filter_at=iatom, filter_dim=2, idx_to_at=idx_to_at_AO)

            CALL dbt_contract(1.0_dp, t_2c_RI_tmp(2), t_3c_tmp(1), &
                              0.0_dp, t_3c_tmp(2), map_1=[1], map_2=[2, 3], &
                              contract_1=[2], notcontract_1=[1], &
                              contract_2=[1], notcontract_2=[2, 3], &
                              filter_eps=ri_data%filter_eps, flop=nflop)
            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            DO i_img = batch_ranges(i_batch), batch_ranges(i_batch + 1) - 1
               CALL unstack_t_3c_apc(t_3c_tmp(3), t_3c_tmp(2), i_img - batch_ranges(i_batch) + 1)
               CALL dbt_copy(t_3c_tmp(3), ri_data%t_3c_int_ctr_1(1, i_img), summation=.TRUE., &
                             order=[2, 1, 3], move_data=.TRUE.)
            END DO
            CALL dbt_clear(t_3c_tmp(1))
         END DO
         CALL dbt_batched_contract_finalize(t_2c_RI_tmp(2))

      END DO
      CALL dbt_destroy(t_2c_RI_tmp(1))
      CALL dbt_destroy(t_2c_RI_tmp(2))
      CALL dbt_destroy(t_3c_tmp(1))
      CALL dbt_destroy(t_3c_tmp(2))
      CALL dbt_destroy(t_3c_tmp(3))

      DO i_img = 1, nimg
         CALL dbt_destroy(t_3c_int(1, i_img))
      END DO

      CALL timestop(handle)

   END SUBROUTINE precontract_3c_ints

! **************************************************************************************************
!> \brief Copy the data of a 2D tensor living in the main MPI group to a sub-group, given the proc
!>        mapping from one to the other (e.g. for a proc idx in the subgroup, we get the idx in the main)
!> \param t2c_sub ...
!> \param t2c_main ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE copy_2c_to_subgroup(t2c_sub, t2c_main, group_size, ngroups, para_env)
      TYPE(dbt_type), INTENT(INOUT)                      :: t2c_sub, t2c_main
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: batch_size, i, i_batch, i_msg, iblk, &
                                                            igroup, iproc, ir, is, jblk, n_batch, &
                                                            nocc, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes1, bsizes2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: block_dest, block_source
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: current_dest
      INTEGER, DIMENSION(2)                              :: ind, nblks
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_req, send_req

      !Stategy: we loop over the main tensor, and send all the data. Then we loop over the sub tensor
      !         and receive it. We do all of it with async MPI communication. The sub tensor needs
      !         to have blocks pre-reserved though

      CALL dbt_get_info(t2c_main, nblks_total=nblks)

      !Loop over the main tensor, count how many blocks are there, which ones, and on which proc
      ALLOCATE (block_source(nblks(1), nblks(2)))
      block_source = -1
      nocc = 0
!$OMP PARALLEL DEFAULT(NONE) SHARED(t2c_main,para_env,nocc,block_source) PRIVATE(iter,ind,blk,found)
      CALL dbt_iterator_start(iter, t2c_main)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_block(t2c_main, ind, blk, found)
         IF (.NOT. found) CYCLE

         block_source(ind(1), ind(2)) = para_env%mepos
!$OMP ATOMIC
         nocc = nocc + 1
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL para_env%sum(nocc)
      CALL para_env%sum(block_source)
      block_source = block_source + para_env%num_pe - 1

      !Loop over the sub tensor, get the block destination
      igroup = para_env%mepos/group_size
      ALLOCATE (block_dest(nblks(1), nblks(2)))
      block_dest = -1
      DO jblk = 1, nblks(2)
         DO iblk = 1, nblks(1)
            IF (block_source(iblk, jblk) == -1) CYCLE

            CALL dbt_get_stored_coordinates(t2c_sub, [iblk, jblk], iproc)
            block_dest(iblk, jblk) = igroup*group_size + iproc !mapping of iproc in subgroup to main group idx
         END DO
      END DO

      ALLOCATE (bsizes1(nblks(1)), bsizes2(nblks(2)))
      CALL dbt_get_info(t2c_main, blk_size_1=bsizes1, blk_size_2=bsizes2)

      ALLOCATE (current_dest(nblks(1), nblks(2), 0:ngroups - 1))
      DO igroup = 0, ngroups - 1
         !for a given subgroup, need to make the destination available to everyone in the main group
         current_dest(:, :, igroup) = block_dest(:, :)
         CALL para_env%bcast(current_dest(:, :, igroup), source=igroup*group_size) !bcast from first proc in sub-group
      END DO

      !Typically, the maximum value for a MPI tag is ~32000, so we go by batches to be sure
      batch_size = 32000
      n_batch = (nocc*ngroups)/batch_size
      IF (MODULO(nocc*ngroups, batch_size) .NE. 0) n_batch = n_batch + 1

      DO i_batch = 1, n_batch
         !Loop over groups, blocks and send/receive
         ALLOCATE (send_buff(batch_size), recv_buff(batch_size))
         ALLOCATE (send_req(batch_size), recv_req(batch_size))
         ir = 0
         is = 0
         i_msg = 0
         DO jblk = 1, nblks(2)
            DO iblk = 1, nblks(1)
               DO igroup = 0, ngroups - 1
                  IF (block_source(iblk, jblk) == -1) CYCLE

                  i_msg = i_msg + 1
                  IF (i_msg < (i_batch - 1)*batch_size + 1 .OR. i_msg > i_batch*batch_size) CYCLE

                  !a unique tag per block, within this batch
                  tag = i_msg - (i_batch - 1)*batch_size

                  found = .FALSE.
                  IF (para_env%mepos == block_source(iblk, jblk)) THEN
                     CALL dbt_get_block(t2c_main, [iblk, jblk], blk, found)
                  END IF

                  !If blocks live on same proc, simply copy. Else MPI send/recv
                  IF (block_source(iblk, jblk) == current_dest(iblk, jblk, igroup)) THEN
                     IF (found) CALL dbt_put_block(t2c_sub, [iblk, jblk], SHAPE(blk), blk)
                  ELSE
                     IF (para_env%mepos == block_source(iblk, jblk) .AND. found) THEN
                        ALLOCATE (send_buff(tag)%array(bsizes1(iblk), bsizes2(jblk)))
                        send_buff(tag)%array(:, :) = blk(:, :)
                        is = is + 1
                        CALL para_env%isend(msgin=send_buff(tag)%array, dest=current_dest(iblk, jblk, igroup), &
                                            request=send_req(is), tag=tag)
                     END IF

                     IF (para_env%mepos == current_dest(iblk, jblk, igroup)) THEN
                        ALLOCATE (recv_buff(tag)%array(bsizes1(iblk), bsizes2(jblk)))
                        ir = ir + 1
                        CALL para_env%irecv(msgout=recv_buff(tag)%array, source=block_source(iblk, jblk), &
                                            request=recv_req(ir), tag=tag)
                     END IF
                  END IF

                  IF (found) DEALLOCATE (blk)
               END DO
            END DO
         END DO

         CALL mp_waitall(send_req(1:is))
         CALL mp_waitall(recv_req(1:ir))
         !clean-up
         DO i = 1, batch_size
            IF (ASSOCIATED(send_buff(i)%array)) DEALLOCATE (send_buff(i)%array)
         END DO

         !Finally copy the data from the buffer to the sub-tensor
         i_msg = 0
         DO jblk = 1, nblks(2)
            DO iblk = 1, nblks(1)
               DO igroup = 0, ngroups - 1
                  IF (block_source(iblk, jblk) == -1) CYCLE

                  i_msg = i_msg + 1
                  IF (i_msg < (i_batch - 1)*batch_size + 1 .OR. i_msg > i_batch*batch_size) CYCLE

                  !a unique tag per block, within this batch
                  tag = i_msg - (i_batch - 1)*batch_size

                  IF (para_env%mepos == current_dest(iblk, jblk, igroup) .AND. &
                      block_source(iblk, jblk) .NE. current_dest(iblk, jblk, igroup)) THEN

                     ALLOCATE (blk(bsizes1(iblk), bsizes2(jblk)))
                     blk(:, :) = recv_buff(tag)%array(:, :)
                     CALL dbt_put_block(t2c_sub, [iblk, jblk], SHAPE(blk), blk)
                     DEALLOCATE (blk)
                  END IF
               END DO
            END DO
         END DO

         !clean-up
         DO i = 1, batch_size
            IF (ASSOCIATED(recv_buff(i)%array)) DEALLOCATE (recv_buff(i)%array)
         END DO
         DEALLOCATE (send_buff, recv_buff, send_req, recv_req)
      END DO !i_batch
      CALL dbt_finalize(t2c_sub)

   END SUBROUTINE copy_2c_to_subgroup

! **************************************************************************************************
!> \brief Copy the data of a 3D tensor living in the main MPI group to a sub-group, given the proc
!>        mapping from one to the other (e.g. for a proc idx in the subgroup, we get the idx in the main)
!> \param t3c_sub ...
!> \param t3c_main ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param iatom_to_subgroup ...
!> \param dim_at ...
!> \param idx_to_at ...
! **************************************************************************************************
   SUBROUTINE copy_3c_to_subgroup(t3c_sub, t3c_main, group_size, ngroups, para_env, iatom_to_subgroup, &
                                  dim_at, idx_to_at)
      TYPE(dbt_type), INTENT(INOUT)                      :: t3c_sub, t3c_main
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_1d_logical_p_type), DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: iatom_to_subgroup
      INTEGER, INTENT(IN), OPTIONAL                      :: dim_at
      INTEGER, DIMENSION(:), OPTIONAL                    :: idx_to_at

      INTEGER                                            :: batch_size, i, i_batch, i_msg, iatom, &
                                                            iblk, igroup, iproc, ir, is, jblk, &
                                                            kblk, n_batch, nocc, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes1, bsizes2, bsizes3
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: block_dest, block_source
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: current_dest
      INTEGER, DIMENSION(3)                              :: ind, nblks
      LOGICAL                                            :: filter_at, found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(cp_3d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_req, send_req

      !Stategy: we loop over the main tensor, and send all the data. Then we loop over the sub tensor
      !         and receive it. We do all of it with async MPI communication. The sub tensor needs
      !         to have blocks pre-reserved though

      CALL dbt_get_info(t3c_main, nblks_total=nblks)

      !in some cases, only copy a fraction of the 3c tensor to a given subgroup (corresponding to some atoms)
      filter_at = .FALSE.
      IF (PRESENT(iatom_to_subgroup) .AND. PRESENT(dim_at) .AND. PRESENT(idx_to_at)) THEN
         filter_at = .TRUE.
         CPASSERT(nblks(dim_at) == SIZE(idx_to_at))
      END IF

      !Loop over the main tensor, count how many blocks are there, which ones, and on which proc
      ALLOCATE (block_source(nblks(1), nblks(2), nblks(3)))
      block_source = -1
      nocc = 0
!$OMP PARALLEL DEFAULT(NONE) SHARED(t3c_main,para_env,nocc,block_source) PRIVATE(iter,ind,blk,found)
      CALL dbt_iterator_start(iter, t3c_main)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         CALL dbt_get_block(t3c_main, ind, blk, found)
         IF (.NOT. found) CYCLE

         block_source(ind(1), ind(2), ind(3)) = para_env%mepos
!$OMP ATOMIC
         nocc = nocc + 1
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL para_env%sum(nocc)
      CALL para_env%sum(block_source)
      block_source = block_source + para_env%num_pe - 1

      !Loop over the sub tensor, get the block destination
      igroup = para_env%mepos/group_size
      ALLOCATE (block_dest(nblks(1), nblks(2), nblks(3)))
      block_dest = -1
      DO kblk = 1, nblks(3)
         DO jblk = 1, nblks(2)
            DO iblk = 1, nblks(1)
               IF (block_source(iblk, jblk, kblk) == -1) CYCLE

               CALL dbt_get_stored_coordinates(t3c_sub, [iblk, jblk, kblk], iproc)
               block_dest(iblk, jblk, kblk) = igroup*group_size + iproc !mapping of iproc in subgroup to main group idx
            END DO
         END DO
      END DO

      ALLOCATE (bsizes1(nblks(1)), bsizes2(nblks(2)), bsizes3(nblks(3)))
      CALL dbt_get_info(t3c_main, blk_size_1=bsizes1, blk_size_2=bsizes2, blk_size_3=bsizes3)

      ALLOCATE (current_dest(nblks(1), nblks(2), nblks(3), 0:ngroups - 1))
      DO igroup = 0, ngroups - 1
         !for a given subgroup, need to make the destination available to everyone in the main group
         current_dest(:, :, :, igroup) = block_dest(:, :, :)
         CALL para_env%bcast(current_dest(:, :, :, igroup), source=igroup*group_size) !bcast from first proc in subgroup
      END DO

      !Typically, the maximum value for a MPI tag is ~32000, so we go by batches to be sure
      batch_size = 32000
      n_batch = (nocc*ngroups)/batch_size
      IF (MODULO(nocc*ngroups, batch_size) .NE. 0) n_batch = n_batch + 1

      DO i_batch = 1, n_batch
         !Loop over groups, blocks and send/receive
         ALLOCATE (send_buff(batch_size), recv_buff(batch_size))
         ALLOCATE (send_req(batch_size), recv_req(batch_size))
         ir = 0
         is = 0
         i_msg = 0
         DO kblk = 1, nblks(3)
            DO jblk = 1, nblks(2)
               DO iblk = 1, nblks(1)
                  DO igroup = 0, ngroups - 1
                     IF (block_source(iblk, jblk, kblk) == -1) CYCLE

                     i_msg = i_msg + 1
                     IF (i_msg < (i_batch - 1)*batch_size + 1 .OR. i_msg > i_batch*batch_size) CYCLE

                     !a unique tag per block, within this batch
                     tag = i_msg - (i_batch - 1)*batch_size

                     IF (filter_at) THEN
                        ind(:) = [iblk, jblk, kblk]
                        iatom = idx_to_at(ind(dim_at))
                        IF (.NOT. iatom_to_subgroup(iatom)%array(igroup + 1)) CYCLE
                     END IF

                     found = .FALSE.
                     IF (para_env%mepos == block_source(iblk, jblk, kblk)) THEN
                        CALL dbt_get_block(t3c_main, [iblk, jblk, kblk], blk, found)
                     END IF

                     !If blocks live on same proc, simply copy. Else MPI send/recv
                     IF (block_source(iblk, jblk, kblk) == current_dest(iblk, jblk, kblk, igroup)) THEN
                        IF (found) CALL dbt_put_block(t3c_sub, [iblk, jblk, kblk], SHAPE(blk), blk)
                     ELSE
                        IF (para_env%mepos == block_source(iblk, jblk, kblk) .AND. found) THEN
                           ALLOCATE (send_buff(tag)%array(bsizes1(iblk), bsizes2(jblk), bsizes3(kblk)))
                           send_buff(tag)%array(:, :, :) = blk(:, :, :)
                           is = is + 1
                           CALL para_env%isend(msgin=send_buff(tag)%array, &
                                               dest=current_dest(iblk, jblk, kblk, igroup), &
                                               request=send_req(is), tag=tag)
                        END IF

                        IF (para_env%mepos == current_dest(iblk, jblk, kblk, igroup)) THEN
                           ALLOCATE (recv_buff(tag)%array(bsizes1(iblk), bsizes2(jblk), bsizes3(kblk)))
                           ir = ir + 1
                           CALL para_env%irecv(msgout=recv_buff(tag)%array, source=block_source(iblk, jblk, kblk), &
                                               request=recv_req(ir), tag=tag)
                        END IF
                     END IF

                     IF (found) DEALLOCATE (blk)
                  END DO
               END DO
            END DO
         END DO

         CALL mp_waitall(send_req(1:is))
         CALL mp_waitall(recv_req(1:ir))
         !clean-up
         DO i = 1, batch_size
            IF (ASSOCIATED(send_buff(i)%array)) DEALLOCATE (send_buff(i)%array)
         END DO

         !Finally copy the data from the buffer to the sub-tensor
         i_msg = 0
         DO kblk = 1, nblks(3)
            DO jblk = 1, nblks(2)
               DO iblk = 1, nblks(1)
                  DO igroup = 0, ngroups - 1
                     IF (block_source(iblk, jblk, kblk) == -1) CYCLE

                     i_msg = i_msg + 1
                     IF (i_msg < (i_batch - 1)*batch_size + 1 .OR. i_msg > i_batch*batch_size) CYCLE

                     !a unique tag per block, within this batch
                     tag = i_msg - (i_batch - 1)*batch_size

                     IF (filter_at) THEN
                        ind(:) = [iblk, jblk, kblk]
                        iatom = idx_to_at(ind(dim_at))
                        IF (.NOT. iatom_to_subgroup(iatom)%array(igroup + 1)) CYCLE
                     END IF

                     IF (para_env%mepos == current_dest(iblk, jblk, kblk, igroup) .AND. &
                         block_source(iblk, jblk, kblk) .NE. current_dest(iblk, jblk, kblk, igroup)) THEN

                        ALLOCATE (blk(bsizes1(iblk), bsizes2(jblk), bsizes3(kblk)))
                        blk(:, :, :) = recv_buff(tag)%array(:, :, :)
                        CALL dbt_put_block(t3c_sub, [iblk, jblk, kblk], SHAPE(blk), blk)
                        DEALLOCATE (blk)
                     END IF
                  END DO
               END DO
            END DO
         END DO

         !clean-up
         DO i = 1, batch_size
            IF (ASSOCIATED(recv_buff(i)%array)) DEALLOCATE (recv_buff(i)%array)
         END DO
         DEALLOCATE (send_buff, recv_buff, send_req, recv_req)
      END DO !i_batch
      CALL dbt_finalize(t3c_sub)

   END SUBROUTINE copy_3c_to_subgroup

! **************************************************************************************************
!> \brief A routine that gather the pieces of the KS matrix accross the subgroup and puts it in the
!>        main group. Each b_img, iatom, jatom tuple is one a single CPU
!> \param ks_t ...
!> \param ks_t_sub ...
!> \param group_size ...
!> \param sparsity_pattern ...
!> \param para_env ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE gather_ks_matrix(ks_t, ks_t_sub, group_size, sparsity_pattern, para_env, ri_data)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: ks_t, ks_t_sub
      INTEGER, INTENT(IN)                                :: group_size
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: sparsity_pattern
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gather_ks_matrix'

      INTEGER                                            :: b_img, dest, handle, i, i_spin, iatom, &
                                                            igroup, ir, is, jatom, n_mess, natom, &
                                                            nimg, nspins, source, tag
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_req, send_req

      CALL timeset(routineN, handle)

      nimg = SIZE(sparsity_pattern, 3)
      natom = SIZE(sparsity_pattern, 2)
      nspins = SIZE(ks_t, 1)

      DO b_img = 1, nimg
         n_mess = 0
         DO i_spin = 1, nspins
            DO jatom = 1, natom
               DO iatom = 1, natom
                  IF (sparsity_pattern(iatom, jatom, b_img) > -1) n_mess = n_mess + 1
               END DO
            END DO
         END DO

         ALLOCATE (send_buff(n_mess), recv_buff(n_mess))
         ALLOCATE (send_req(n_mess), recv_req(n_mess))
         ir = 0
         is = 0
         n_mess = 0
         tag = 0

         DO i_spin = 1, nspins
            DO jatom = 1, natom
               DO iatom = 1, natom
                  IF (sparsity_pattern(iatom, jatom, b_img) < 0) CYCLE
                  n_mess = n_mess + 1
                  tag = tag + 1

                  !sending the message
                  CALL dbt_get_stored_coordinates(ks_t(i_spin, b_img), [iatom, jatom], dest)
                  CALL dbt_get_stored_coordinates(ks_t_sub(i_spin, b_img), [iatom, jatom], source) !source within sub
                  igroup = sparsity_pattern(iatom, jatom, b_img)
                  source = source + igroup*group_size
                  IF (para_env%mepos == source) THEN
                     CALL dbt_get_block(ks_t_sub(i_spin, b_img), [iatom, jatom], blk, found)
                     IF (source == dest) THEN
                        IF (found) CALL dbt_put_block(ks_t(i_spin, b_img), [iatom, jatom], SHAPE(blk), blk)
                     ELSE
                        ALLOCATE (send_buff(n_mess)%array(ri_data%bsizes_AO(iatom), ri_data%bsizes_AO(jatom)))
                        send_buff(n_mess)%array(:, :) = 0.0_dp
                        IF (found) THEN
                           send_buff(n_mess)%array(:, :) = blk(:, :)
                        END IF
                        is = is + 1
                        CALL para_env%isend(msgin=send_buff(n_mess)%array, dest=dest, &
                                            request=send_req(is), tag=tag)
                     END IF
                     DEALLOCATE (blk)
                  END IF

                  !receiving the message
                  IF (para_env%mepos == dest .AND. source .NE. dest) THEN
                     ALLOCATE (recv_buff(n_mess)%array(ri_data%bsizes_AO(iatom), ri_data%bsizes_AO(jatom)))
                     ir = ir + 1
                     CALL para_env%irecv(msgout=recv_buff(n_mess)%array, source=source, &
                                         request=recv_req(ir), tag=tag)
                  END IF
               END DO !iatom
            END DO !jatom
         END DO !ispin

         CALL mp_waitall(send_req(1:is))
         CALL mp_waitall(recv_req(1:ir))

         !Copy the messages received into the KS matrix
         n_mess = 0
         DO i_spin = 1, nspins
            DO jatom = 1, natom
               DO iatom = 1, natom
                  IF (sparsity_pattern(iatom, jatom, b_img) < 0) CYCLE
                  n_mess = n_mess + 1

                  CALL dbt_get_stored_coordinates(ks_t(i_spin, b_img), [iatom, jatom], dest)
                  IF (para_env%mepos == dest) THEN
                     IF (.NOT. ASSOCIATED(recv_buff(n_mess)%array)) CYCLE
                     ALLOCATE (blk(ri_data%bsizes_AO(iatom), ri_data%bsizes_AO(jatom)))
                     blk(:, :) = recv_buff(n_mess)%array(:, :)
                     CALL dbt_put_block(ks_t(i_spin, b_img), [iatom, jatom], SHAPE(blk), blk)
                     DEALLOCATE (blk)
                  END IF
               END DO
            END DO
         END DO

         !clean-up
         DO i = 1, n_mess
            IF (ASSOCIATED(send_buff(i)%array)) DEALLOCATE (send_buff(i)%array)
            IF (ASSOCIATED(recv_buff(i)%array)) DEALLOCATE (recv_buff(i)%array)
         END DO
         DEALLOCATE (send_buff, recv_buff, send_req, recv_req)
      END DO !b_img

      CALL timestop(handle)

   END SUBROUTINE gather_ks_matrix

! **************************************************************************************************
!> \brief copy all required 2c tensors from the main MPI group to the subgroups
!> \param mat_2c_pot ...
!> \param t_2c_work ...
!> \param t_2c_ao_tmp ...
!> \param ks_t_split ...
!> \param ks_t_sub ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE get_subgroup_2c_tensors(mat_2c_pot, t_2c_work, t_2c_ao_tmp, ks_t_split, ks_t_sub, &
                                      group_size, ngroups, para_env, para_env_sub, ri_data)
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: mat_2c_pot
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_work, t_2c_ao_tmp, ks_t_split
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: ks_t_sub
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER :: routineN = 'get_subgroup_2c_tensors'

      INTEGER                                            :: handle, i, i_img, i_RI, i_spin, iproc, &
                                                            j, natom, nblks, nimg, nspins
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_RI_ext_split, &
                                                            dist1, dist2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, RI_blk_size, row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: dbcsr_pgrid
      REAL(dp)                                           :: occ
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist_sub
      TYPE(dbt_pgrid_type)                               :: pgrid_2d
      TYPE(dbt_type)                                     :: work, work_sub

      CALL timeset(routineN, handle)

      !Create the 2d pgrid
      pdims_2d = 0
      CALL dbt_pgrid_create(para_env_sub, pdims_2d, pgrid_2d)

      natom = SIZE(ri_data%bsizes_RI)
      nblks = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*natom))
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         bsizes_RI_ext_split((i_RI - 1)*nblks + 1:i_RI*nblks) = ri_data%bsizes_RI_split(:)
      END DO

      !nRI x nRI 2c tensors
      CALL create_2c_tensor(t_2c_work(1), dist1, dist2, pgrid_2d, &
                            bsizes_RI_ext, bsizes_RI_ext, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      CALL create_2c_tensor(t_2c_work(2), dist1, dist2, pgrid_2d, &
                            bsizes_RI_ext_split, bsizes_RI_ext_split, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      !the AO based tensors
      CALL create_2c_tensor(ks_t_split(1), dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_create(ks_t_split(1), ks_t_split(2))

      CALL create_2c_tensor(t_2c_ao_tmp(1), dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_AO, ri_data%bsizes_AO, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      nspins = SIZE(ks_t_sub, 1)
      nimg = SIZE(ks_t_sub, 2)
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(t_2c_ao_tmp(1), ks_t_sub(i_spin, i_img))
         END DO
      END DO

      !Finally the HFX potential matrices
      !For now, we do a convoluted things where we go to tensors first, then back to matrices.
      CALL create_2c_tensor(work_sub, dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_RI, ri_data%bsizes_RI, &
                            name="(RI | RI)")
      CALL dbt_create(ri_data%kp_mat_2c_pot(1, 1), work)

      ALLOCATE (dbcsr_pgrid(0:pdims_2d(1) - 1, 0:pdims_2d(2) - 1))
      iproc = 0
      DO i = 0, pdims_2d(1) - 1
         DO j = 0, pdims_2d(2) - 1
            dbcsr_pgrid(i, j) = iproc
            iproc = iproc + 1
         END DO
      END DO

      !We need to have the same exact 2d block dist as the tensors
      ALLOCATE (col_dist(natom), row_dist(natom))
      row_dist(:) = dist1(:)
      col_dist(:) = dist2(:)

      ALLOCATE (RI_blk_size(natom))
      RI_blk_size(:) = ri_data%bsizes_RI(:)

      CALL dbcsr_distribution_new(dbcsr_dist_sub, group=para_env_sub%get_handle(), pgrid=dbcsr_pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(mat_2c_pot(1), dist=dbcsr_dist_sub, name="sub", matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=RI_blk_size, col_blk_size=RI_blk_size)

      DO i_img = 1, nimg
         IF (i_img > 1) CALL dbcsr_create(mat_2c_pot(i_img), template=mat_2c_pot(1))
         CALL dbt_copy_matrix_to_tensor(ri_data%kp_mat_2c_pot(1, i_img), work)
         CALL get_tensor_occupancy(work, nze, occ)
         IF (nze == 0) CYCLE

         CALL copy_2c_to_subgroup(work_sub, work, group_size, ngroups, para_env)
         CALL dbt_copy_tensor_to_matrix(work_sub, mat_2c_pot(i_img))
         CALL dbcsr_filter(mat_2c_pot(i_img), ri_data%filter_eps)
         CALL dbt_clear(work_sub)
      END DO

      CALL dbt_destroy(work)
      CALL dbt_destroy(work_sub)
      CALL dbt_pgrid_destroy(pgrid_2d)
      CALL dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE (col_dist, row_dist, RI_blk_size, dbcsr_pgrid)
      CALL timestop(handle)

   END SUBROUTINE get_subgroup_2c_tensors

! **************************************************************************************************
!> \brief copy all required 3c tensors from the main MPI group to the subgroups
!> \param t_3c_int ...
!> \param t_3c_work_2 ...
!> \param t_3c_work_3 ...
!> \param t_3c_apc ...
!> \param t_3c_apc_sub ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE get_subgroup_3c_tensors(t_3c_int, t_3c_work_2, t_3c_work_3, t_3c_apc, t_3c_apc_sub, &
                                      group_size, ngroups, para_env, para_env_sub, ri_data)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_int, t_3c_work_2, t_3c_work_3
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_apc, t_3c_apc_sub
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER :: routineN = 'get_subgroup_3c_tensors'

      INTEGER                                            :: batch_size, handle, handle2, i_img, &
                                                            i_RI, i_spin, ib, natom, nblks_AO, &
                                                            nblks_RI, nimg, nspins
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_RI_ext_split, &
                                                            bsizes_stack, dist1, dist2, dist3, &
                                                            dist_stack, idx_to_at
      INTEGER, DIMENSION(3)                              :: pdims
      REAL(dp)                                           :: occ
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(dbt_type)                                     :: tmp, work_atom_block, work_atom_block_sub

      CALL timeset(routineN, handle)

      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext_split((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI_split(:)
      END DO

      !Preparing larger block sizes for efficient communication (less, bigger messages)
      nblks_RI = SIZE(ri_data%bsizes_RI)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*nblks_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI(:)
      END DO

      batch_size = ri_data%kp_stack_size
      nblks_AO = SIZE(ri_data%bsizes_AO_split)
      ALLOCATE (bsizes_stack(batch_size*nblks_AO))
      DO ib = 1, batch_size
         bsizes_stack((ib - 1)*nblks_AO + 1:ib*nblks_AO) = ri_data%bsizes_AO_split(:)
      END DO

      !Create the pgrid for the configuration correspoinding to ri_data%t_3c_int_ctr_3
      natom = SIZE(ri_data%bsizes_RI)
      pdims = 0
      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid, &
                            tensor_dims=[SIZE(bsizes_RI_ext_split), natom, batch_size*SIZE(ri_data%bsizes_AO_split)])

      !Create all required 3c tensors in that configuration
      CALL create_3c_tensor(t_3c_int(1), dist1, dist2, dist3, &
                            pgrid, bsizes_RI_ext_split, ri_data%bsizes_AO_split, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(RI | AO AO)")
      nimg = SIZE(t_3c_int)
      DO i_img = 2, nimg
         CALL dbt_create(t_3c_int(1), t_3c_int(i_img))
      END DO

      !The stacked work tensors, in a distribution that matches that of t_3c_int
      ALLOCATE (dist_stack(batch_size*nblks_AO))
      DO ib = 1, batch_size
         dist_stack((ib - 1)*nblks_AO + 1:ib*nblks_AO) = dist3(:)
      END DO

      CALL dbt_distribution_new(t_dist, pgrid, dist1, dist2, dist_stack)
      CALL dbt_create(t_3c_work_3(1), "work_3_stack", t_dist, [1], [2, 3], &
                      bsizes_RI_ext_split, ri_data%bsizes_AO_split, bsizes_stack)
      CALL dbt_create(t_3c_work_3(1), t_3c_work_3(2))
      CALL dbt_create(t_3c_work_3(1), t_3c_work_3(3))
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist3, dist_stack)

      !For more efficient communication, we use intermediate tensors with larger block size
      CALL create_3c_tensor(work_atom_block_sub, dist1, dist2, dist3, &
                            pgrid, bsizes_RI_ext, ri_data%bsizes_AO, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(work_atom_block, dist1, dist2, dist3, &
                            ri_data%pgrid, bsizes_RI_ext, ri_data%bsizes_AO, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      !Finally copy the integrals into the subgroups (if not there already)
      CALL timeset(routineN//"_ints", handle2)
      IF (ALLOCATED(ri_data%kp_t_3c_int)) THEN
         DO i_img = 1, nimg
            CALL dbt_copy(ri_data%kp_t_3c_int(i_img), t_3c_int(i_img), move_data=.TRUE.)
         END DO
      ELSE
         ALLOCATE (ri_data%kp_t_3c_int(nimg))
         DO i_img = 1, nimg
            CALL dbt_create(t_3c_int(i_img), ri_data%kp_t_3c_int(i_img))
            CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_1(1, i_img), nze, occ)
            IF (nze == 0) CYCLE
            CALL dbt_copy(ri_data%t_3c_int_ctr_1(1, i_img), work_atom_block, order=[2, 1, 3])
            CALL copy_3c_to_subgroup(work_atom_block_sub, work_atom_block, group_size, ngroups, para_env)
            CALL dbt_copy(work_atom_block_sub, t_3c_int(i_img), move_data=.TRUE.)
         END DO
      END IF
      CALL timestop(handle2)
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_destroy(work_atom_block)
      CALL dbt_destroy(work_atom_block_sub)

      !Do the same for the t_3c_ctr_2 configuration
      pdims = 0
      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid, &
                            tensor_dims=[1, SIZE(bsizes_RI_ext_split), batch_size*SIZE(ri_data%bsizes_AO_split)])

      !For more efficient communication, we use intermediate tensors with larger block size
      CALL create_3c_tensor(work_atom_block_sub, dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO, bsizes_RI_ext, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(work_atom_block, dist1, dist2, dist3, &
                            ri_data%pgrid_1, ri_data%bsizes_AO, bsizes_RI_ext, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      !template for t_3c_apc_sub
      CALL create_3c_tensor(tmp, dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO_split, bsizes_RI_ext_split, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(AO RI | AO)")

      !create t_3c_work_2 tensors in a distribution that matches the above
      ALLOCATE (dist_stack(batch_size*nblks_AO))
      DO ib = 1, batch_size
         dist_stack((ib - 1)*nblks_AO + 1:ib*nblks_AO) = dist3(:)
      END DO

      CALL dbt_distribution_new(t_dist, pgrid, dist1, dist2, dist_stack)
      CALL dbt_create(t_3c_work_2(1), "work_2_stack", t_dist, [1], [2, 3], &
                      ri_data%bsizes_AO_split, bsizes_RI_ext_split, bsizes_stack)
      CALL dbt_create(t_3c_work_2(1), t_3c_work_2(2))
      CALL dbt_create(t_3c_work_2(1), t_3c_work_2(3))
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist3, dist_stack)

      !Finally copy data from t_3c_apc to the subgroups
      ALLOCATE (idx_to_at(SIZE(ri_data%bsizes_AO)))
      CALL get_idx_to_atom(idx_to_at, ri_data%bsizes_AO, ri_data%bsizes_AO)
      nspins = SIZE(t_3c_apc, 1)
      CALL timeset(routineN//"_apc", handle2)
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(tmp, t_3c_apc_sub(i_spin, i_img))
            CALL get_tensor_occupancy(t_3c_apc(i_spin, i_img), nze, occ)
            IF (nze == 0) CYCLE
            CALL dbt_copy(t_3c_apc(i_spin, i_img), work_atom_block, move_data=.TRUE.)
            CALL copy_3c_to_subgroup(work_atom_block_sub, work_atom_block, group_size, &
                                     ngroups, para_env, ri_data%iatom_to_subgroup, 1, idx_to_at)
            CALL dbt_copy(work_atom_block_sub, t_3c_apc_sub(i_spin, i_img), move_data=.TRUE.)
         END DO
         DO i_spin = 1, nspins
            CALL dbt_destroy(t_3c_apc(i_spin, i_img))
         END DO
      END DO
      CALL timestop(handle2)
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_destroy(tmp)
      CALL dbt_destroy(work_atom_block)
      CALL dbt_destroy(work_atom_block_sub)

      CALL timestop(handle)

   END SUBROUTINE get_subgroup_3c_tensors

! **************************************************************************************************
!> \brief copy all required 2c force tensors from the main MPI group to the subgroups
!> \param t_2c_inv ...
!> \param t_2c_bint ...
!> \param t_2c_metric ...
!> \param mat_2c_pot ...
!> \param t_2c_work ...
!> \param rho_ao_t ...
!> \param rho_ao_t_sub ...
!> \param t_2c_der_metric ...
!> \param t_2c_der_metric_sub ...
!> \param mat_der_pot ...
!> \param mat_der_pot_sub ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ri_data ...
!> \note Main MPI group tensors are deleted within this routine, for memory optimization
! **************************************************************************************************
   SUBROUTINE get_subgroup_2c_derivs(t_2c_inv, t_2c_bint, t_2c_metric, mat_2c_pot, t_2c_work, rho_ao_t, &
                                     rho_ao_t_sub, t_2c_der_metric, t_2c_der_metric_sub, mat_der_pot, &
                                     mat_der_pot_sub, group_size, ngroups, para_env, para_env_sub, ri_data)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_inv, t_2c_bint, t_2c_metric
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: mat_2c_pot
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_work
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: rho_ao_t, rho_ao_t_sub, t_2c_der_metric, &
                                                            t_2c_der_metric_sub
      TYPE(dbcsr_type), DIMENSION(:, :), INTENT(INOUT)   :: mat_der_pot, mat_der_pot_sub
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER :: routineN = 'get_subgroup_2c_derivs'

      INTEGER                                            :: handle, i, i_img, i_RI, i_spin, i_xyz, &
                                                            iatom, iproc, j, natom, nblks, nimg, &
                                                            nspins
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_RI_ext_split, &
                                                            dist1, dist2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, RI_blk_size, row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: dbcsr_pgrid
      REAL(dp)                                           :: occ
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist_sub
      TYPE(dbt_pgrid_type)                               :: pgrid_2d
      TYPE(dbt_type)                                     :: work, work_sub

      CALL timeset(routineN, handle)

      !Note: a fair portion of this routine is copied from the energy version of it
      !Create the 2d pgrid
      pdims_2d = 0
      CALL dbt_pgrid_create(para_env_sub, pdims_2d, pgrid_2d)

      natom = SIZE(ri_data%bsizes_RI)
      nblks = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*natom))
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         bsizes_RI_ext_split((i_RI - 1)*nblks + 1:i_RI*nblks) = ri_data%bsizes_RI_split(:)
      END DO

      !nRI x nRI 2c tensors
      CALL create_2c_tensor(t_2c_inv(1), dist1, dist2, pgrid_2d, &
                            bsizes_RI_ext, bsizes_RI_ext, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      CALL dbt_create(t_2c_inv(1), t_2c_bint(1))
      CALL dbt_create(t_2c_inv(1), t_2c_metric(1))
      DO iatom = 2, natom
         CALL dbt_create(t_2c_inv(1), t_2c_inv(iatom))
         CALL dbt_create(t_2c_inv(1), t_2c_bint(iatom))
         CALL dbt_create(t_2c_inv(1), t_2c_metric(iatom))
      END DO
      CALL dbt_create(t_2c_inv(1), t_2c_work(1))
      CALL dbt_create(t_2c_inv(1), t_2c_work(2))
      CALL dbt_create(t_2c_inv(1), t_2c_work(3))
      CALL dbt_create(t_2c_inv(1), t_2c_work(4))

      CALL create_2c_tensor(t_2c_work(5), dist1, dist2, pgrid_2d, &
                            bsizes_RI_ext_split, bsizes_RI_ext_split, &
                            name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      !copy the data from the main group.
      DO iatom = 1, natom
         CALL copy_2c_to_subgroup(t_2c_inv(iatom), ri_data%t_2c_inv(1, iatom), group_size, ngroups, para_env)
         CALL copy_2c_to_subgroup(t_2c_bint(iatom), ri_data%t_2c_int(1, iatom), group_size, ngroups, para_env)
         CALL copy_2c_to_subgroup(t_2c_metric(iatom), ri_data%t_2c_pot(1, iatom), group_size, ngroups, para_env)
      END DO

      !This includes the derivatives of the RI metric, for which there is one per atom
      DO i_xyz = 1, 3
         DO iatom = 1, natom
            CALL dbt_create(t_2c_inv(1), t_2c_der_metric_sub(iatom, i_xyz))
            CALL copy_2c_to_subgroup(t_2c_der_metric_sub(iatom, i_xyz), t_2c_der_metric(iatom, i_xyz), &
                                     group_size, ngroups, para_env)
            CALL dbt_destroy(t_2c_der_metric(iatom, i_xyz))
         END DO
      END DO

      !AO x AO 2c tensors
      CALL create_2c_tensor(rho_ao_t_sub(1, 1), dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)
      nspins = SIZE(rho_ao_t, 1)
      nimg = SIZE(rho_ao_t, 2)

      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            IF (.NOT. (i_img == 1 .AND. i_spin == 1)) &
               CALL dbt_create(rho_ao_t_sub(1, 1), rho_ao_t_sub(i_spin, i_img))
            CALL copy_2c_to_subgroup(rho_ao_t_sub(i_spin, i_img), rho_ao_t(i_spin, i_img), &
                                     group_size, ngroups, para_env)
            CALL dbt_destroy(rho_ao_t(i_spin, i_img))
         END DO
      END DO

      !The RIxRI matrices, going through tensors
      CALL create_2c_tensor(work_sub, dist1, dist2, pgrid_2d, &
                            ri_data%bsizes_RI, ri_data%bsizes_RI, &
                            name="(RI | RI)")
      CALL dbt_create(ri_data%kp_mat_2c_pot(1, 1), work)

      ALLOCATE (dbcsr_pgrid(0:pdims_2d(1) - 1, 0:pdims_2d(2) - 1))
      iproc = 0
      DO i = 0, pdims_2d(1) - 1
         DO j = 0, pdims_2d(2) - 1
            dbcsr_pgrid(i, j) = iproc
            iproc = iproc + 1
         END DO
      END DO

      !We need to have the same exact 2d block dist as the tensors
      ALLOCATE (col_dist(natom), row_dist(natom))
      row_dist(:) = dist1(:)
      col_dist(:) = dist2(:)

      ALLOCATE (RI_blk_size(natom))
      RI_blk_size(:) = ri_data%bsizes_RI(:)

      CALL dbcsr_distribution_new(dbcsr_dist_sub, group=para_env_sub%get_handle(), pgrid=dbcsr_pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(mat_2c_pot(1), dist=dbcsr_dist_sub, name="sub", matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=RI_blk_size, col_blk_size=RI_blk_size)

      !The HFX potential
      DO i_img = 1, nimg
         IF (i_img > 1) CALL dbcsr_create(mat_2c_pot(i_img), template=mat_2c_pot(1))
         CALL dbt_copy_matrix_to_tensor(ri_data%kp_mat_2c_pot(1, i_img), work)
         CALL get_tensor_occupancy(work, nze, occ)
         IF (nze == 0) CYCLE

         CALL copy_2c_to_subgroup(work_sub, work, group_size, ngroups, para_env)
         CALL dbt_copy_tensor_to_matrix(work_sub, mat_2c_pot(i_img))
         CALL dbcsr_filter(mat_2c_pot(i_img), ri_data%filter_eps)
         CALL dbt_clear(work_sub)
      END DO

      !The derivatives of the HFX potential
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbcsr_create(mat_der_pot_sub(i_img, i_xyz), template=mat_2c_pot(1))
            CALL dbt_copy_matrix_to_tensor(mat_der_pot(i_img, i_xyz), work)
            CALL dbcsr_release(mat_der_pot(i_img, i_xyz))
            CALL get_tensor_occupancy(work, nze, occ)
            IF (nze == 0) CYCLE

            CALL copy_2c_to_subgroup(work_sub, work, group_size, ngroups, para_env)
            CALL dbt_copy_tensor_to_matrix(work_sub, mat_der_pot_sub(i_img, i_xyz))
            CALL dbcsr_filter(mat_der_pot_sub(i_img, i_xyz), ri_data%filter_eps)
            CALL dbt_clear(work_sub)
         END DO
      END DO

      CALL dbt_destroy(work)
      CALL dbt_destroy(work_sub)
      CALL dbt_pgrid_destroy(pgrid_2d)
      CALL dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE (col_dist, row_dist, RI_blk_size, dbcsr_pgrid)

      CALL timestop(handle)

   END SUBROUTINE get_subgroup_2c_derivs

! **************************************************************************************************
!> \brief copy all required 3c derivative tensors from the main MPI group to the subgroups
!> \param t_3c_work_2 ...
!> \param t_3c_work_3 ...
!> \param t_3c_der_AO ...
!> \param t_3c_der_AO_sub ...
!> \param t_3c_der_RI ...
!> \param t_3c_der_RI_sub ...
!> \param t_3c_apc ...
!> \param t_3c_apc_sub ...
!> \param t_3c_der_stack ...
!> \param group_size ...
!> \param ngroups ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ri_data ...
!> \note the tensor containing the derivatives in the main MPI group are deleted for memory
! **************************************************************************************************
   SUBROUTINE get_subgroup_3c_derivs(t_3c_work_2, t_3c_work_3, t_3c_der_AO, t_3c_der_AO_sub, &
                                     t_3c_der_RI, t_3c_der_RI_sub, t_3c_apc, t_3c_apc_sub, &
                                     t_3c_der_stack, group_size, ngroups, para_env, para_env_sub, &
                                     ri_data)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_work_2, t_3c_work_3
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_der_AO, t_3c_der_AO_sub, &
                                                            t_3c_der_RI, t_3c_der_RI_sub, &
                                                            t_3c_apc, t_3c_apc_sub
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_der_stack
      INTEGER, INTENT(IN)                                :: group_size, ngroups
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(len=*), PARAMETER :: routineN = 'get_subgroup_3c_derivs'

      INTEGER                                            :: batch_size, handle, i_img, i_RI, i_spin, &
                                                            i_xyz, ib, nblks_AO, nblks_RI, nimg, &
                                                            nspins, pdims(3)
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_RI_ext_split, &
                                                            bsizes_stack, dist1, dist2, dist3, &
                                                            dist_stack, idx_to_at
      REAL(dp)                                           :: occ
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(dbt_type)                                     :: tmp, work_atom_block, work_atom_block_sub

      CALL timeset(routineN, handle)

      !We use intermediate tensors with larger block size for more optimized communication
      nblks_RI = SIZE(ri_data%bsizes_RI)
      ALLOCATE (bsizes_RI_ext(ri_data%ncell_RI*nblks_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI(:)
      END DO

      CALL dbt_get_info(ri_data%kp_t_3c_int(1), pdims=pdims)
      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid)

      CALL create_3c_tensor(work_atom_block_sub, dist1, dist2, dist3, &
                            pgrid, bsizes_RI_ext, ri_data%bsizes_AO, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(work_atom_block, dist1, dist2, dist3, &
                            ri_data%pgrid_2, bsizes_RI_ext, ri_data%bsizes_AO, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)
      CALL dbt_pgrid_destroy(pgrid)

      !We use the 3c integrals on the subgroup as template for the derivatives
      nimg = ri_data%nimg
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbt_create(ri_data%kp_t_3c_int(1), t_3c_der_AO_sub(i_img, i_xyz))
            CALL get_tensor_occupancy(t_3c_der_AO(i_img, i_xyz), nze, occ)
            IF (nze == 0) CYCLE

            CALL dbt_copy(t_3c_der_AO(i_img, i_xyz), work_atom_block, move_data=.TRUE.)
            CALL copy_3c_to_subgroup(work_atom_block_sub, work_atom_block, &
                                     group_size, ngroups, para_env)
            CALL dbt_copy(work_atom_block_sub, t_3c_der_AO_sub(i_img, i_xyz), move_data=.TRUE.)
         END DO

         DO i_img = 1, nimg
            CALL dbt_create(ri_data%kp_t_3c_int(1), t_3c_der_RI_sub(i_img, i_xyz))
            CALL get_tensor_occupancy(t_3c_der_RI(i_img, i_xyz), nze, occ)
            IF (nze == 0) CYCLE

            CALL dbt_copy(t_3c_der_RI(i_img, i_xyz), work_atom_block, move_data=.TRUE.)
            CALL copy_3c_to_subgroup(work_atom_block_sub, work_atom_block, &
                                     group_size, ngroups, para_env)
            CALL dbt_copy(work_atom_block_sub, t_3c_der_RI_sub(i_img, i_xyz), move_data=.TRUE.)
         END DO

         DO i_img = 1, nimg
            CALL dbt_destroy(t_3c_der_RI(i_img, i_xyz))
            CALL dbt_destroy(t_3c_der_AO(i_img, i_xyz))
         END DO
      END DO
      CALL dbt_destroy(work_atom_block_sub)
      CALL dbt_destroy(work_atom_block)

      !Deal with t_3c_apc
      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (bsizes_RI_ext_split(ri_data%ncell_RI*nblks_RI))
      DO i_RI = 1, ri_data%ncell_RI
         bsizes_RI_ext_split((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI_split(:)
      END DO

      pdims = 0
      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid, &
                            tensor_dims=[1, SIZE(bsizes_RI_ext_split), batch_size*SIZE(ri_data%bsizes_AO_split)])

      CALL create_3c_tensor(work_atom_block_sub, dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO, bsizes_RI_ext, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(work_atom_block, dist1, dist2, dist3, &
                            ri_data%pgrid_1, ri_data%bsizes_AO, bsizes_RI_ext, &
                            ri_data%bsizes_AO, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(tmp, dist1, dist2, dist3, &
                            pgrid, ri_data%bsizes_AO_split, bsizes_RI_ext_split, &
                            ri_data%bsizes_AO_split, [1], [2, 3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      ALLOCATE (idx_to_at(SIZE(ri_data%bsizes_AO)))
      CALL get_idx_to_atom(idx_to_at, ri_data%bsizes_AO, ri_data%bsizes_AO)
      nspins = SIZE(t_3c_apc, 1)
      DO i_img = 1, nimg
         DO i_spin = 1, nspins
            CALL dbt_create(tmp, t_3c_apc_sub(i_spin, i_img))
            CALL get_tensor_occupancy(t_3c_apc(i_spin, i_img), nze, occ)
            IF (nze == 0) CYCLE
            CALL dbt_copy(t_3c_apc(i_spin, i_img), work_atom_block, move_data=.TRUE.)
            CALL copy_3c_to_subgroup(work_atom_block_sub, work_atom_block, group_size, &
                                     ngroups, para_env, ri_data%iatom_to_subgroup, 1, idx_to_at)
            CALL dbt_copy(work_atom_block_sub, t_3c_apc_sub(i_spin, i_img), move_data=.TRUE.)
         END DO
         DO i_spin = 1, nspins
            CALL dbt_destroy(t_3c_apc(i_spin, i_img))
         END DO
      END DO
      CALL dbt_destroy(tmp)
      CALL dbt_destroy(work_atom_block)
      CALL dbt_destroy(work_atom_block_sub)
      CALL dbt_pgrid_destroy(pgrid)

      !t_3c_work_3 based on structure of 3c integrals/derivs
      batch_size = ri_data%kp_stack_size
      nblks_AO = SIZE(ri_data%bsizes_AO_split)
      ALLOCATE (bsizes_stack(batch_size*nblks_AO))
      DO ib = 1, batch_size
         bsizes_stack((ib - 1)*nblks_AO + 1:ib*nblks_AO) = ri_data%bsizes_AO_split(:)
      END DO

      ALLOCATE (dist1(ri_data%ncell_RI*nblks_RI), dist2(nblks_AO), dist3(nblks_AO))
      CALL dbt_get_info(ri_data%kp_t_3c_int(1), proc_dist_1=dist1, proc_dist_2=dist2, &
                        proc_dist_3=dist3, pdims=pdims)

      ALLOCATE (dist_stack(batch_size*nblks_AO))
      DO ib = 1, batch_size
         dist_stack((ib - 1)*nblks_AO + 1:ib*nblks_AO) = dist3(:)
      END DO

      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid)
      CALL dbt_distribution_new(t_dist, pgrid, dist1, dist2, dist_stack)
      CALL dbt_create(t_3c_work_3(1), "work_3_stack", t_dist, [1], [2, 3], &
                      bsizes_RI_ext_split, ri_data%bsizes_AO_split, bsizes_stack)
      CALL dbt_create(t_3c_work_3(1), t_3c_work_3(2))
      CALL dbt_create(t_3c_work_3(1), t_3c_work_3(3))
      CALL dbt_create(t_3c_work_3(1), t_3c_work_3(4))
      CALL dbt_distribution_destroy(t_dist)
      CALL dbt_pgrid_destroy(pgrid)
      DEALLOCATE (dist1, dist2, dist3, dist_stack)

      !the derivatives are stacked in the same way
      CALL dbt_create(t_3c_work_3(1), t_3c_der_stack(1))
      CALL dbt_create(t_3c_work_3(1), t_3c_der_stack(2))
      CALL dbt_create(t_3c_work_3(1), t_3c_der_stack(3))
      CALL dbt_create(t_3c_work_3(1), t_3c_der_stack(4))
      CALL dbt_create(t_3c_work_3(1), t_3c_der_stack(5))
      CALL dbt_create(t_3c_work_3(1), t_3c_der_stack(6))

      !t_3c_work_2 based on structure of t_3c_apc
      ALLOCATE (dist1(nblks_AO), dist2(ri_data%ncell_RI*nblks_RI), dist3(nblks_AO))
      CALL dbt_get_info(t_3c_apc_sub(1, 1), proc_dist_1=dist1, proc_dist_2=dist2, &
                        proc_dist_3=dist3, pdims=pdims)

      ALLOCATE (dist_stack(batch_size*nblks_AO))
      DO ib = 1, batch_size
         dist_stack((ib - 1)*nblks_AO + 1:ib*nblks_AO) = dist3(:)
      END DO

      CALL dbt_pgrid_create(para_env_sub, pdims, pgrid)
      CALL dbt_distribution_new(t_dist, pgrid, dist1, dist2, dist_stack)
      CALL dbt_create(t_3c_work_2(1), "work_3_stack", t_dist, [1], [2, 3], &
                      ri_data%bsizes_AO_split, bsizes_RI_ext_split, bsizes_stack)
      CALL dbt_create(t_3c_work_2(1), t_3c_work_2(2))
      CALL dbt_create(t_3c_work_2(1), t_3c_work_2(3))
      CALL dbt_distribution_destroy(t_dist)
      CALL dbt_pgrid_destroy(pgrid)
      DEALLOCATE (dist1, dist2, dist3, dist_stack)

      CALL timestop(handle)

   END SUBROUTINE get_subgroup_3c_derivs

! **************************************************************************************************
!> \brief A routine that reorders the t_3c_int tensors such that all items which are fully empty
!>        are bunched together. This way, we can get much more efficient screening based on NZE
!> \param t_3c_ints ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE reorder_3c_ints(t_3c_ints, ri_data)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_ints
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reorder_3c_ints'

      INTEGER                                            :: handle, i_img, idx, idx_empty, idx_full, &
                                                            nimg
      INTEGER(int_8)                                     :: nze
      REAL(dp)                                           :: occ
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_3c_tmp

      CALL timeset(routineN, handle)

      nimg = ri_data%nimg
      ALLOCATE (t_3c_tmp(nimg))
      DO i_img = 1, nimg
         CALL dbt_create(t_3c_ints(i_img), t_3c_tmp(i_img))
         CALL dbt_copy(t_3c_ints(i_img), t_3c_tmp(i_img), move_data=.TRUE.)
      END DO

      !Loop over the images, check if ints have NZE == 0, and put them at the start or end of the
      !initial tensor array. Keep the mapping in an array
      ALLOCATE (ri_data%idx_to_img(nimg))
      idx_full = 0
      idx_empty = nimg + 1

      DO i_img = 1, nimg
         CALL get_tensor_occupancy(t_3c_tmp(i_img), nze, occ)
         IF (nze == 0) THEN
            idx_empty = idx_empty - 1
            CALL dbt_copy(t_3c_tmp(i_img), t_3c_ints(idx_empty), move_data=.TRUE.)
            ri_data%idx_to_img(idx_empty) = i_img
         ELSE
            idx_full = idx_full + 1
            CALL dbt_copy(t_3c_tmp(i_img), t_3c_ints(idx_full), move_data=.TRUE.)
            ri_data%idx_to_img(idx_full) = i_img
         END IF
         CALL dbt_destroy(t_3c_tmp(i_img))
      END DO

      !store the highest image index with non-zero integrals
      ri_data%nimg_nze = idx_full

      ALLOCATE (ri_data%img_to_idx(nimg))
      DO idx = 1, nimg
         ri_data%img_to_idx(ri_data%idx_to_img(idx)) = idx
      END DO

      CALL timestop(handle)

   END SUBROUTINE reorder_3c_ints

! **************************************************************************************************
!> \brief A routine that reorders the 3c derivatives, the same way that the integrals are, also to
!>        increase efficiency of screening
!> \param t_3c_derivs ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE reorder_3c_derivs(t_3c_derivs, ri_data)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_derivs
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reorder_3c_derivs'

      INTEGER                                            :: handle, i_img, i_xyz, idx, nimg
      INTEGER(int_8)                                     :: nze
      REAL(dp)                                           :: occ
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_3c_tmp

      CALL timeset(routineN, handle)

      nimg = ri_data%nimg
      ALLOCATE (t_3c_tmp(nimg))
      DO i_img = 1, nimg
         CALL dbt_create(t_3c_derivs(1, 1), t_3c_tmp(i_img))
      END DO

      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbt_copy(t_3c_derivs(i_img, i_xyz), t_3c_tmp(i_img), move_data=.TRUE.)
         END DO
         DO i_img = 1, nimg
            idx = ri_data%img_to_idx(i_img)
            CALL dbt_copy(t_3c_tmp(i_img), t_3c_derivs(idx, i_xyz), move_data=.TRUE.)
            CALL get_tensor_occupancy(t_3c_derivs(idx, i_xyz), nze, occ)
            IF (nze > 0) ri_data%nimg_nze = MAX(idx, ri_data%nimg_nze)
         END DO
      END DO

      DO i_img = 1, nimg
         CALL dbt_destroy(t_3c_tmp(i_img))
      END DO

      CALL timestop(handle)

   END SUBROUTINE reorder_3c_derivs

! **************************************************************************************************
!> \brief Get the sparsity pattern related to the non-symmetric AO basis overlap neighbor list
!> \param pattern ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_sparsity_pattern(pattern, ri_data, qs_env)
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: pattern
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: iatom, j_img, jatom, mj_img, natom, nimg
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bins
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: tmp_pattern
      INTEGER, DIMENSION(3)                              :: cell_j
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c

      NULLIFY (nl_2c, nl_iterator, kpoints, cell_to_index, dft_control, index_to_cell, para_env)

      CALL get_qs_env(qs_env, kpoints=kpoints, dft_control=dft_control, para_env=para_env, natom=natom)
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell, sab_nl=nl_2c)

      nimg = ri_data%nimg
      pattern(:, :, :) = 0

      !We use the symmetric nl for all images that have an opposite cell
      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_j)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         mj_img = get_opp_index(j_img, qs_env)
         IF (mj_img > nimg .OR. mj_img < 1) CYCLE

         IF (ri_data%present_images(j_img) == 0) CYCLE

         pattern(iatom, jatom, j_img) = 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      !If there is no opposite cell present, then we take into account the non-symmetric nl
      CALL get_kpoint_info(kpoints, sab_nl_nosym=nl_2c)

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_j)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         mj_img = get_opp_index(j_img, qs_env)
         IF (mj_img .LE. nimg .AND. mj_img > 0) CYCLE

         IF (ri_data%present_images(j_img) == 0) CYCLE

         pattern(iatom, jatom, j_img) = 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL para_env%sum(pattern)

      !If the opposite image is considered, then there is no need to compute diagonal twice
      DO j_img = 2, nimg
         DO iatom = 1, natom
            IF (pattern(iatom, iatom, j_img) .NE. 0) THEN
               mj_img = get_opp_index(j_img, qs_env)
               IF (mj_img > nimg .OR. mj_img < 1) CYCLE
               pattern(iatom, iatom, mj_img) = 0
            END IF
         END DO
      END DO

      ! We want to equilibrate the sparsity pattern such that there are same amount of blocks
      ! for each atom i of i,j pairs
      ALLOCATE (bins(natom))
      bins(:) = 0

      ALLOCATE (tmp_pattern(natom, natom, nimg))
      tmp_pattern(:, :, :) = 0
      DO j_img = 1, nimg
         DO jatom = 1, natom
            DO iatom = 1, natom
               IF (pattern(iatom, jatom, j_img) == 0) CYCLE
               mj_img = get_opp_index(j_img, qs_env)

               !Should we take the i,j,b or th j,i,-b atomic block?
               IF (mj_img > nimg .OR. mj_img < 1) THEN
                  !No opposite image, no choice
                  bins(iatom) = bins(iatom) + 1
                  tmp_pattern(iatom, jatom, j_img) = 1
               ELSE

                  IF (bins(iatom) > bins(jatom)) THEN
                     bins(jatom) = bins(jatom) + 1
                     tmp_pattern(jatom, iatom, mj_img) = 1
                  ELSE
                     bins(iatom) = bins(iatom) + 1
                     tmp_pattern(iatom, jatom, j_img) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO

      ! -1 => unoccupied, 0 => occupied
      pattern(:, :, :) = tmp_pattern(:, :, :) - 1

   END SUBROUTINE get_sparsity_pattern

! **************************************************************************************************
!> \brief Distribute the iatom, jatom, b_img triplet over the subgroupd to spread the load
!>        the group id for each triplet is passed as the value of sparsity_pattern(i, j, b),
!>        with -1 being an unoccupied block
!> \param sparsity_pattern ...
!> \param ngroups ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE get_sub_dist(sparsity_pattern, ngroups, ri_data)
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: sparsity_pattern
      INTEGER, INTENT(IN)                                :: ngroups
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: b_img, ctr, iat, iatom, igroup, jatom, &
                                                            natom, nimg, ub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: max_at_per_group
      REAL(dp)                                           :: cost
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: bins

      natom = SIZE(sparsity_pattern, 2)
      nimg = SIZE(sparsity_pattern, 3)

      !To avoid unnecessary data replication accross the subgroups, we want to have a limited number
      !of subgroup with the data of a given iatom. At the minimum, all groups have 1 atom
      !We assume that the cost associated to each iatom is roughly the same
      IF (.NOT. ALLOCATED(ri_data%iatom_to_subgroup)) THEN
         ALLOCATE (ri_data%iatom_to_subgroup(natom), max_at_per_group(ngroups))
         DO iatom = 1, natom
            NULLIFY (ri_data%iatom_to_subgroup(iatom)%array)
            ALLOCATE (ri_data%iatom_to_subgroup(iatom)%array(ngroups))
            ri_data%iatom_to_subgroup(iatom)%array(:) = .FALSE.
         END DO

         ub = natom/ngroups
         IF (ub*ngroups < natom) ub = ub + 1
         max_at_per_group(:) = MAX(1, ub)

         !We want each atom to be present the same amount of times. Some groups might have more atoms
         !than other to achieve this.
         ctr = 0
         DO WHILE (MODULO(SUM(max_at_per_group), natom) .NE. 0)
            igroup = MODULO(ctr, ngroups) + 1
            max_at_per_group(igroup) = max_at_per_group(igroup) + 1
            ctr = ctr + 1
         END DO

         ctr = 0
         DO igroup = 1, ngroups
            DO iat = 1, max_at_per_group(igroup)
               iatom = MODULO(ctr, natom) + 1
               ri_data%iatom_to_subgroup(iatom)%array(igroup) = .TRUE.
               ctr = ctr + 1
            END DO
         END DO
      END IF

      ALLOCATE (bins(ngroups))
      bins = 0.0_dp
      DO b_img = 1, nimg
         DO jatom = 1, natom
            DO iatom = 1, natom
               IF (sparsity_pattern(iatom, jatom, b_img) == -1) CYCLE
               igroup = MINLOC(bins, 1, MASK=ri_data%iatom_to_subgroup(iatom)%array) - 1

               !Use cost information from previous SCF if available
               IF (ANY(ri_data%kp_cost > EPSILON(0.0_dp))) THEN
                  cost = ri_data%kp_cost(iatom, jatom, b_img)
               ELSE
                  cost = REAL(ri_data%bsizes_AO(iatom)*ri_data%bsizes_AO(jatom), dp)
               END IF
               bins(igroup + 1) = bins(igroup + 1) + cost
               sparsity_pattern(iatom, jatom, b_img) = igroup
            END DO
         END DO
      END DO

   END SUBROUTINE get_sub_dist

! **************************************************************************************************
!> \brief A rouine that updates the sparsity pattern for force calculation, where all i,j,b combinations
!>        are visited.
!> \param force_pattern ...
!> \param scf_pattern ...
!> \param ngroups ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE update_pattern_to_forces(force_pattern, scf_pattern, ngroups, ri_data, qs_env)
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: force_pattern, scf_pattern
      INTEGER, INTENT(IN)                                :: ngroups
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: b_img, iatom, igroup, jatom, mb_img, &
                                                            natom, nimg
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: bins

      natom = SIZE(scf_pattern, 2)
      nimg = SIZE(scf_pattern, 3)

      ALLOCATE (bins(ngroups))
      bins = 0.0_dp

      DO b_img = 1, nimg
         mb_img = get_opp_index(b_img, qs_env)
         DO jatom = 1, natom
            DO iatom = 1, natom
               !Important: same distribution as KS matrix, because reuse t_3c_apc
               igroup = MINLOC(bins, 1, MASK=ri_data%iatom_to_subgroup(iatom)%array) - 1

               !check that block not already treated
               IF (scf_pattern(iatom, jatom, b_img) > -1) CYCLE

               !If not, take the cost of block j, i, -b (same energy contribution)
               IF (mb_img > 0 .AND. mb_img .LE. nimg) THEN
                  IF (scf_pattern(jatom, iatom, mb_img) == -1) CYCLE
                  bins(igroup + 1) = bins(igroup + 1) + ri_data%kp_cost(jatom, iatom, mb_img)
                  force_pattern(iatom, jatom, b_img) = igroup
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE update_pattern_to_forces

! **************************************************************************************************
!> \brief A routine that determines the extend of the KP RI-HFX periodic images, including for the
!>        extension of the RI basis
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_kp_and_ri_images(ri_data, qs_env)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kp_and_ri_images'

      INTEGER :: cell_j(3), cell_k(3), handle, i_img, iatom, ikind, j_img, jatom, jcell, katom, &
         kcell, kp_index_lbounds(3), kp_index_ubounds(3), natom, ngroups, nimg, nkind, pcoord(3), &
         pdims(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_AO_1, dist_AO_2, dist_RI, &
                                                            nRI_per_atom, present_img, RI_cells
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(dp)                                           :: bump_fact, dij, dik, fact, image_range, &
                                                            RI_range, rij(3), rik(3)
      TYPE(dbt_type)                                     :: t_dummy
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_AO, basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_cart_type)                                 :: mp_comm_t3c
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: hfx_section

      NULLIFY (qs_kind_set, dist_2d, nl_2c, nl_iterator, dft_control, &
               particle_set, kpoints, para_env, cell_to_index, hfx_section)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, distribution_2d=dist_2d, &
                      dft_control=dft_control, particle_set=particle_set, kpoints=kpoints, &
                      para_env=para_env, natom=natom)
      nimg = dft_control%nimages
      CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      kp_index_lbounds = LBOUND(cell_to_index)
      kp_index_ubounds = UBOUND(cell_to_index)

      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF%RI")
      CALL section_vals_val_get(hfx_section, "KP_RI_EXTENSION_FACTOR", r_val=fact)
      CALL section_vals_val_get(hfx_section, "KP_NGROUPS", i_val=ngroups)

      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)

      !In case of shortrange HFX potential, it is imprtant to be consistent with the rest of the KP
      !code, and use EPS_SCHWARZ to determine the range (rather than eps_filter_2c in normal RI-HFX)
      IF (ri_data%hfx_pot%potential_type == do_potential_short) THEN
         CALL erfc_cutoff(ri_data%eps_schwarz, ri_data%hfx_pot%omega, ri_data%hfx_pot%cutoff_radius)
      END IF

      !Determine the range for contributing periodic images, and for the RI basis extension
      ri_data%kp_RI_range = 0.0_dp
      ri_data%kp_image_range = 0.0_dp
      DO ikind = 1, nkind

         CALL init_interaction_radii_orb_basis(basis_set_AO(ikind)%gto_basis_set, fact*ri_data%eps_pgf_orb)
         CALL get_gto_basis_set(basis_set_AO(ikind)%gto_basis_set, kind_radius=RI_range)
         ri_data%kp_RI_range = MAX(RI_range, ri_data%kp_RI_range)

         CALL init_interaction_radii_orb_basis(basis_set_AO(ikind)%gto_basis_set, ri_data%eps_pgf_orb)
         CALL init_interaction_radii_orb_basis(basis_set_RI(ikind)%gto_basis_set, ri_data%eps_pgf_orb)
         CALL get_gto_basis_set(basis_set_RI(ikind)%gto_basis_set, kind_radius=image_range)

         image_range = 2.0_dp*image_range + cutoff_screen_factor*ri_data%hfx_pot%cutoff_radius
         ri_data%kp_image_range = MAX(image_range, ri_data%kp_image_range)
      END DO

      CALL section_vals_val_get(hfx_section, "KP_RI_BUMP_FACTOR", r_val=bump_fact)
      ri_data%kp_bump_rad = bump_fact*ri_data%kp_RI_range

      !For the extent of the KP RI-HFX images, we are limited by the RI-HFX potential in
      !(mu^0 sigma^a|P^0) (P^0|Q^b) (Q^b|nu^b lambda^a+c), if there is no contact between
      !any P^0 and Q^b, then image b does not contribute
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%hfx_pot, &
                                   "HFX_2c_nl_RI", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)

      ALLOCATE (present_img(nimg))
      present_img = 0
      ri_data%nimg = 0
      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, r=rij, cell=cell_j)

         dij = NORM2(rij)

         j_img = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (j_img > nimg .OR. j_img < 1) CYCLE

         IF (dij > ri_data%kp_image_range) CYCLE

         ri_data%nimg = MAX(j_img, ri_data%nimg)
         present_img(j_img) = 1

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL release_neighbor_list_sets(nl_2c)
      CALL para_env%max(ri_data%nimg)
      IF (ri_data%nimg > nimg) &
         CPABORT("Make sure the smallest exponent of the RI-HFX basis is larger than that of the ORB basis.")

      !Keep track of which images will not contribute, so that can be ignored before calculation
      CALL para_env%sum(present_img)
      ALLOCATE (ri_data%present_images(ri_data%nimg))
      ri_data%present_images = 0
      DO i_img = 1, ri_data%nimg
         IF (present_img(i_img) > 0) ri_data%present_images(i_img) = 1
      END DO

      CALL create_3c_tensor(t_dummy, dist_AO_1, dist_AO_2, dist_RI, &
                            ri_data%pgrid, ri_data%bsizes_AO, ri_data%bsizes_AO, ri_data%bsizes_RI, &
                            map1=[1, 2], map2=[3], name="(AO AO | RI)")

      CALL dbt_mp_environ_pgrid(ri_data%pgrid, pdims, pcoord)
      CALL mp_comm_t3c%create(ri_data%pgrid%mp_comm_2d, 3, pdims)
      CALL distribution_3d_create(dist_3d, dist_AO_1, dist_AO_2, dist_RI, &
                                  nkind, particle_set, mp_comm_t3c, own_comm=.TRUE.)
      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
      CALL dbt_destroy(t_dummy)

      !For the extension of the RI basis P in (mu^0 sigma^a |P^i), we consider an atom if the distance,
      !between mu^0 and P^i if smaller or equal to the kind radius of mu^0
      CALL build_3c_neighbor_lists(nl_3c, basis_set_AO, basis_set_AO, basis_set_RI, dist_3d, &
                                   ri_data%ri_metric, "HFX_3c_nl", qs_env, op_pos=2, sym_ij=.FALSE., &
                                   own_dist=.TRUE.)

      ALLOCATE (RI_cells(nimg))
      RI_cells = 0

      ALLOCATE (nRI_per_atom(natom))
      nRI_per_atom = 0

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, cell_k=cell_k, rik=rik, cell_j=cell_j, &
                                   iatom=iatom, jatom=jatom, katom=katom)
         dik = NORM2(rik)

         IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
             ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

         jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
         IF (jcell > nimg .OR. jcell < 1) CYCLE

         IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
             ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

         kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
         IF (kcell > nimg .OR. kcell < 1) CYCLE

         IF (dik > ri_data%kp_RI_range) CYCLE
         RI_cells(kcell) = 1

         IF (jcell == 1 .AND. iatom == jatom) nRI_per_atom(iatom) = nRI_per_atom(iatom) + ri_data%bsizes_RI(katom)
      END DO
      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)
      CALL neighbor_list_3c_destroy(nl_3c)
      CALL para_env%sum(RI_cells)
      CALL para_env%sum(nRI_per_atom)

      ALLOCATE (ri_data%img_to_RI_cell(nimg))
      ri_data%ncell_RI = 0
      ri_data%img_to_RI_cell = 0
      DO i_img = 1, nimg
         IF (RI_cells(i_img) > 0) THEN
            ri_data%ncell_RI = ri_data%ncell_RI + 1
            ri_data%img_to_RI_cell(i_img) = ri_data%ncell_RI
         END IF
      END DO

      ALLOCATE (ri_data%RI_cell_to_img(ri_data%ncell_RI))
      DO i_img = 1, nimg
         IF (ri_data%img_to_RI_cell(i_img) > 0) ri_data%RI_cell_to_img(ri_data%img_to_RI_cell(i_img)) = i_img
      END DO

      !Print some info
      IF (ri_data%unit_nr > 0) THEN
         WRITE (ri_data%unit_nr, FMT="(/T3,A,I29)") &
            "KP-HFX_RI_INFO| Number of RI-KP parallel groups:", ngroups
         WRITE (ri_data%unit_nr, FMT="(T3,A,F35.3)") &
            "KP-HFX_RI_INFO| RI basis extension factor:", fact
         WRITE (ri_data%unit_nr, FMT="(T3,A,F31.3,A)") &
            "KP-HFX_RI_INFO| RI basis extension radius:", ri_data%kp_RI_range*angstrom, " Ang"
         WRITE (ri_data%unit_nr, FMT="(T3,A,F12.3,A, F6.3, A)") &
            "KP-HFX_RI_INFO| RI basis bump factor and bump radius:", bump_fact, " /", &
            ri_data%kp_bump_rad*angstrom, " Ang"
         WRITE (ri_data%unit_nr, FMT="(T3,A,I16,A)") &
            "KP-HFX_RI_INFO| The extended RI bases cover up to ", ri_data%ncell_RI, " unit cells"
         WRITE (ri_data%unit_nr, FMT="(T3,A,I18)") &
            "KP-HFX_RI_INFO| Average number of sgf in extended RI bases:", SUM(nRI_per_atom)/natom
         WRITE (ri_data%unit_nr, FMT="(T3,A,F13.3,A)") &
            "KP-HFX_RI_INFO| Consider all image cells within a radius of ", ri_data%kp_image_range*angstrom, " Ang"
         WRITE (ri_data%unit_nr, FMT="(T3,A,I27/)") &
            "KP-HFX_RI_INFO| Number of image cells considered: ", ri_data%nimg
         CALL m_flush(ri_data%unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_kp_and_ri_images

! **************************************************************************************************
!> \brief A routine that creates tensors structure for rho_ao and 3c_ints in a stacked format for
!>        the efficient contractions of rho_sigma^0,lambda^c * (mu^0 sigam^a | P) => TAS tensors
!> \param res_stack ...
!> \param rho_stack ...
!> \param ints_stack ...
!> \param rho_template ...
!> \param ints_template ...
!> \param stack_size ...
!> \param ri_data ...
!> \param qs_env ...
!> \note The result tensor has the exact same shape and distribution as the integral tensor
! **************************************************************************************************
   SUBROUTINE get_stack_tensors(res_stack, rho_stack, ints_stack, rho_template, ints_template, &
                                stack_size, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: res_stack, rho_stack, ints_stack
      TYPE(dbt_type), INTENT(INOUT)                      :: rho_template, ints_template
      INTEGER, INTENT(IN)                                :: stack_size
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: is, nblks, nblks_3c(3), pdims_3d(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_RI_ext, bsizes_stack, dist1, &
                                                            dist2, dist3, dist_stack1, &
                                                            dist_stack2, dist_stack3
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env)

      CALL get_qs_env(qs_env, para_env=para_env)

      nblks = SIZE(ri_data%bsizes_AO_split)
      ALLOCATE (bsizes_stack(stack_size*nblks))
      DO is = 1, stack_size
         bsizes_stack((is - 1)*nblks + 1:is*nblks) = ri_data%bsizes_AO_split(:)
      END DO

      ALLOCATE (dist1(nblks), dist2(nblks), dist_stack1(stack_size*nblks), dist_stack2(stack_size*nblks))
      CALL dbt_get_info(rho_template, proc_dist_1=dist1, proc_dist_2=dist2)
      DO is = 1, stack_size
         dist_stack1((is - 1)*nblks + 1:is*nblks) = dist1(:)
         dist_stack2((is - 1)*nblks + 1:is*nblks) = dist2(:)
      END DO

      !First 2c tensor matches the distribution of template
      !It is stacked in both directions
      CALL dbt_distribution_new(t_dist, ri_data%pgrid_2d, dist_stack1, dist_stack2)
      CALL dbt_create(rho_stack(1), "RHO_stack", t_dist, [1], [2], bsizes_stack, bsizes_stack)
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist_stack1, dist_stack2)

      !Second 2c tensor has optimal distribution on the 2d pgrid
      CALL create_2c_tensor(rho_stack(2), dist1, dist2, ri_data%pgrid_2d, bsizes_stack, bsizes_stack, name="RHO_stack")
      DEALLOCATE (dist1, dist2)

      CALL dbt_get_info(ints_template, nblks_total=nblks_3c)
      ALLOCATE (dist1(nblks_3c(1)), dist2(nblks_3c(2)), dist3(nblks_3c(3)))
      ALLOCATE (dist_stack3(stack_size*nblks_3c(3)), bsizes_RI_ext(nblks_3c(2)))
      CALL dbt_get_info(ints_template, proc_dist_1=dist1, proc_dist_2=dist2, &
                        proc_dist_3=dist3, blk_size_2=bsizes_RI_ext)
      DO is = 1, stack_size
         dist_stack3((is - 1)*nblks_3c(3) + 1:is*nblks_3c(3)) = dist3(:)
      END DO

      !First 3c tensor matches the distribution of template
      CALL dbt_distribution_new(t_dist, ri_data%pgrid_1, dist1, dist2, dist_stack3)
      CALL dbt_create(ints_stack(1), "ints_stack", t_dist, [1, 2], [3], ri_data%bsizes_AO_split, &
                      bsizes_RI_ext, bsizes_stack)
      CALL dbt_distribution_destroy(t_dist)
      DEALLOCATE (dist1, dist2, dist3, dist_stack3)

      !Second 3c tensor has optimal pgrid
      pdims_3d = 0
      CALL dbt_pgrid_create(para_env, pdims_3d, pgrid, tensor_dims=[nblks_3c(1), nblks_3c(2), stack_size*nblks_3c(3)])
      CALL create_3c_tensor(ints_stack(2), dist1, dist2, dist3, pgrid, ri_data%bsizes_AO_split, &
                            bsizes_RI_ext, bsizes_stack, [1, 2], [3], name="ints_stack")
      DEALLOCATE (dist1, dist2, dist3)
      CALL dbt_pgrid_destroy(pgrid)

      !The result tensor has the same shape and dist as the integral tensor
      CALL dbt_create(ints_stack(1), res_stack(1))
      CALL dbt_create(ints_stack(2), res_stack(2))

   END SUBROUTINE get_stack_tensors

! **************************************************************************************************
!> \brief Fill the stack of 3c tensors accrding to the order in the images input
!> \param t_3c_stack ...
!> \param t_3c_in ...
!> \param images ...
!> \param stack_dim ...
!> \param ri_data ...
!> \param filter_at ...
!> \param filter_dim ...
!> \param idx_to_at ...
!> \param img_bounds ...
! **************************************************************************************************
   SUBROUTINE fill_3c_stack(t_3c_stack, t_3c_in, images, stack_dim, ri_data, filter_at, filter_dim, &
                            idx_to_at, img_bounds)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_stack
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_3c_in
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: images
      INTEGER, INTENT(IN)                                :: stack_dim
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN), OPTIONAL                      :: filter_at, filter_dim
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL     :: idx_to_at
      INTEGER, INTENT(IN), OPTIONAL                      :: img_bounds(2)

      INTEGER                                            :: dest(3), i_img, idx, ind(3), lb, nblks, &
                                                            nimg, offset, ub
      LOGICAL                                            :: do_filter, found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(dbt_iterator_type)                            :: iter

      !We loop over the a images from the ac_pairs, then copy the 3c ints to the correct spot in
      !in the stack tensor (corresponding to pair index). Distributions match by construction
      nimg = ri_data%nimg
      nblks = SIZE(ri_data%bsizes_AO_split)

      do_filter = .FALSE.
      IF (PRESENT(filter_at) .AND. PRESENT(filter_dim) .AND. PRESENT(idx_to_at)) do_filter = .TRUE.

      lb = 1
      ub = nimg
      offset = 0
      IF (PRESENT(img_bounds)) THEN
         lb = img_bounds(1)
         ub = img_bounds(2) - 1
         offset = lb - 1
      END IF

      DO idx = lb, ub
         i_img = images(idx)
         IF (i_img == 0 .OR. i_img > nimg) CYCLE

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(idx,i_img,t_3c_in,t_3c_stack,nblks,stack_dim,filter_at,filter_dim,idx_to_at,do_filter,offset) &
!$OMP PRIVATE(iter,ind,blk,found,dest)
         CALL dbt_iterator_start(iter, t_3c_in(i_img))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)
            CALL dbt_get_block(t_3c_in(i_img), ind, blk, found)
            IF (.NOT. found) CYCLE

            IF (do_filter) THEN
               IF (.NOT. idx_to_at(ind(filter_dim)) == filter_at) CYCLE
            END IF

            IF (stack_dim == 1) THEN
               dest = [(idx - offset - 1)*nblks + ind(1), ind(2), ind(3)]
            ELSE IF (stack_dim == 2) THEN
               dest = [ind(1), (idx - offset - 1)*nblks + ind(2), ind(3)]
            ELSE
               dest = [ind(1), ind(2), (idx - offset - 1)*nblks + ind(3)]
            END IF

            CALL dbt_put_block(t_3c_stack, dest, SHAPE(blk), blk)
            DEALLOCATE (blk)
         END DO
         CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      END DO !i_img
      CALL dbt_finalize(t_3c_stack)

   END SUBROUTINE fill_3c_stack

! **************************************************************************************************
!> \brief Fill the stack of 2c tensors based on the content of images input
!> \param t_2c_stack ...
!> \param t_2c_in ...
!> \param images ...
!> \param stack_dim ...
!> \param ri_data ...
!> \param img_bounds ...
!> \param shift ...
! **************************************************************************************************
   SUBROUTINE fill_2c_stack(t_2c_stack, t_2c_in, images, stack_dim, ri_data, img_bounds, shift)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_stack
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_in
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: images
      INTEGER, INTENT(IN)                                :: stack_dim
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN), OPTIONAL                      :: img_bounds(2), shift

      INTEGER                                            :: dest(2), i_img, idx, ind(2), lb, &
                                                            my_shift, nblks, nimg, offset, ub
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: blk
      TYPE(dbt_iterator_type)                            :: iter

      !We loop over the a images from the ac_pairs, then copy the 3c ints to the correct spot in
      !in the stack tensor (corresponding to pair index). Distributions match by construction
      nimg = ri_data%nimg
      nblks = SIZE(ri_data%bsizes_AO_split)

      lb = 1
      ub = nimg
      offset = 0
      IF (PRESENT(img_bounds)) THEN
         lb = img_bounds(1)
         ub = img_bounds(2) - 1
         offset = lb - 1
      END IF

      my_shift = 1
      IF (PRESENT(shift)) my_shift = shift

      DO idx = lb, ub
         i_img = images(idx)
         IF (i_img == 0 .OR. i_img > nimg) CYCLE

!$OMP PARALLEL DEFAULT(NONE) SHARED(idx,i_img,t_2c_in,t_2c_stack,nblks,stack_dim,offset,my_shift) &
!$OMP PRIVATE(iter,ind,blk,found,dest)
         CALL dbt_iterator_start(iter, t_2c_in(i_img))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)
            CALL dbt_get_block(t_2c_in(i_img), ind, blk, found)
            IF (.NOT. found) CYCLE

            IF (stack_dim == 1) THEN
               dest = [(idx - offset - 1)*nblks + ind(1), (my_shift - 1)*nblks + ind(2)]
            ELSE
               dest = [(my_shift - 1)*nblks + ind(1), (idx - offset - 1)*nblks + ind(2)]
            END IF

            CALL dbt_put_block(t_2c_stack, dest, SHAPE(blk), blk)
            DEALLOCATE (blk)
         END DO
         CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      END DO !idx
      CALL dbt_finalize(t_2c_stack)

   END SUBROUTINE fill_2c_stack

! **************************************************************************************************
!> \brief Unstacks a stacked 3c tensor containing t_3c_apc
!> \param t_3c_apc ...
!> \param t_stacked ...
!> \param idx ...
! **************************************************************************************************
   SUBROUTINE unstack_t_3c_apc(t_3c_apc, t_stacked, idx)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_apc, t_stacked
      INTEGER, INTENT(IN)                                :: idx

      INTEGER                                            :: current_idx
      INTEGER, DIMENSION(3)                              :: ind, nblks_3c
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(dbt_iterator_type)                            :: iter

      !Note: t_3c_apc and t_stacked must have the same ditribution
      CALL dbt_get_info(t_3c_apc, nblks_total=nblks_3c)

!$OMP PARALLEL DEFAULT(NONE) SHARED(t_3c_apc,t_stacked,idx,nblks_3c) PRIVATE(iter,ind,blk,found,current_idx)
      CALL dbt_iterator_start(iter, t_stacked)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)

         !tensor is stacked along the 3rd dimension
         current_idx = (ind(3) - 1)/nblks_3c(3) + 1
         IF (.NOT. idx == current_idx) CYCLE

         CALL dbt_get_block(t_stacked, ind, blk, found)
         IF (.NOT. found) CYCLE

         CALL dbt_put_block(t_3c_apc, [ind(1), ind(2), ind(3) - (idx - 1)*nblks_3c(3)], SHAPE(blk), blk)
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

   END SUBROUTINE unstack_t_3c_apc

! **************************************************************************************************
!> \brief copies the 3c integrals correspoinding to a single atom mu from the general (P^0| mu^0 sigam^a)
!> \param t_3c_at ...
!> \param t_3c_ints ...
!> \param iatom ...
!> \param dim_at ...
!> \param idx_to_at ...
! **************************************************************************************************
   SUBROUTINE get_atom_3c_ints(t_3c_at, t_3c_ints, iatom, dim_at, idx_to_at)
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_at, t_3c_ints
      INTEGER, INTENT(IN)                                :: iatom, dim_at
      INTEGER, DIMENSION(:), INTENT(IN)                  :: idx_to_at

      INTEGER, DIMENSION(3)                              :: ind
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: blk
      TYPE(dbt_iterator_type)                            :: iter

!$OMP PARALLEL DEFAULT(NONE) SHARED(t_3c_ints,t_3c_at,iatom,idx_to_at,dim_at) PRIVATE(iter,ind,blk,found)
      CALL dbt_iterator_start(iter, t_3c_ints)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)
         IF (.NOT. idx_to_at(ind(dim_at)) == iatom) CYCLE

         CALL dbt_get_block(t_3c_ints, ind, blk, found)
         IF (.NOT. found) CYCLE

         CALL dbt_put_block(t_3c_at, ind, SHAPE(blk), blk)
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      CALL dbt_finalize(t_3c_at)

   END SUBROUTINE get_atom_3c_ints

! **************************************************************************************************
!> \brief Precalculate the 3c and 2c derivatives tensors
!> \param t_3c_der_RI ...
!> \param t_3c_der_AO ...
!> \param mat_der_pot ...
!> \param t_2c_der_metric ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE precalc_derivatives(t_3c_der_RI, t_3c_der_AO, mat_der_pot, t_2c_der_metric, ri_data, qs_env)
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_der_RI, t_3c_der_AO
      TYPE(dbcsr_type), DIMENSION(:, :), INTENT(INOUT)   :: mat_der_pot
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_2c_der_metric
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precalc_derivatives'

      INTEGER                                            :: handle, handle2, i_img, i_mem, i_RI, &
                                                            i_xyz, iatom, n_mem, natom, nblks_RI, &
                                                            ncell_RI, nimg, nkind, nthreads
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bsizes_RI_ext, bsizes_RI_ext_split, dist_AO_1, &
         dist_AO_2, dist_RI, dist_RI_ext, dummy_end, dummy_start, end_blocks, start_blocks
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      REAL(dp)                                           :: occ
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type)                                   :: dbcsr_template
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:, :)     :: mat_der_metric
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: pgrid
      TYPE(dbt_type)                                     :: t_3c_template
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :, :)    :: t_3c_der_AO_prv, t_3c_der_RI_prv
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_AO, basis_set_RI
      TYPE(mp_cart_type)                                 :: mp_comm_t3c
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, dist_2d, nl_2c, particle_set, dft_control, para_env, row_bsize, col_bsize)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, distribution_2d=dist_2d, natom=natom, &
                      particle_set=particle_set, dft_control=dft_control, para_env=para_env)

      nimg = ri_data%nimg
      ncell_RI = ri_data%ncell_RI

      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, basis=basis_set_AO)

      !Dealing with the 3c derivatives
      nthreads = 1
!$    nthreads = omp_get_num_threads()
      pdims = 0
      CALL dbt_pgrid_create(para_env, pdims, pgrid, tensor_dims=[MAX(1, natom/(ri_data%n_mem*nthreads)), natom, natom])

      CALL create_3c_tensor(t_3c_template, dist_AO_1, dist_AO_2, dist_RI, pgrid, &
                            ri_data%bsizes_AO, ri_data%bsizes_AO, ri_data%bsizes_RI, &
                            map1=[1, 2], map2=[3], name="tmp")
      CALL dbt_destroy(t_3c_template)

      !We stack the RI basis images. Keep consistent distribution
      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      ALLOCATE (dist_RI_ext(natom*ncell_RI))
      ALLOCATE (bsizes_RI_ext(natom*ncell_RI))
      ALLOCATE (bsizes_RI_ext_split(nblks_RI*ncell_RI))
      DO i_RI = 1, ncell_RI
         bsizes_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = ri_data%bsizes_RI(:)
         dist_RI_ext((i_RI - 1)*natom + 1:i_RI*natom) = dist_RI(:)
         bsizes_RI_ext_split((i_RI - 1)*nblks_RI + 1:i_RI*nblks_RI) = ri_data%bsizes_RI_split(:)
      END DO

      CALL dbt_distribution_new(t_dist, pgrid, dist_AO_1, dist_AO_2, dist_RI_ext)
      CALL dbt_create(t_3c_template, "KP_3c_der", t_dist, [1, 2], [3], &
                      ri_data%bsizes_AO, ri_data%bsizes_AO, bsizes_RI_ext)
      CALL dbt_distribution_destroy(t_dist)

      ALLOCATE (t_3c_der_RI_prv(nimg, 1, 3), t_3c_der_AO_prv(nimg, 1, 3))
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbt_create(t_3c_template, t_3c_der_RI_prv(i_img, 1, i_xyz))
            CALL dbt_create(t_3c_template, t_3c_der_AO_prv(i_img, 1, i_xyz))
         END DO
      END DO
      CALL dbt_destroy(t_3c_template)

      CALL dbt_mp_environ_pgrid(pgrid, pdims, pcoord)
      CALL mp_comm_t3c%create(pgrid%mp_comm_2d, 3, pdims)
      CALL distribution_3d_create(dist_3d, dist_AO_1, dist_AO_2, dist_RI, &
                                  nkind, particle_set, mp_comm_t3c, own_comm=.TRUE.)
      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
      CALL dbt_pgrid_destroy(pgrid)

      CALL build_3c_neighbor_lists(nl_3c, basis_set_AO, basis_set_AO, basis_set_RI, dist_3d, ri_data%ri_metric, &
                                   "HFX_3c_nl", qs_env, op_pos=2, sym_jk=.FALSE., own_dist=.TRUE.)

      n_mem = ri_data%n_mem
      CALL create_tensor_batches(ri_data%bsizes_RI, n_mem, dummy_start, dummy_end, &
                                 start_blocks, end_blocks)
      DEALLOCATE (dummy_start, dummy_end)

      CALL create_3c_tensor(t_3c_template, dist_RI, dist_AO_1, dist_AO_2, ri_data%pgrid_2, &
                            bsizes_RI_ext_split, ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            map1=[1], map2=[2, 3], name="der (RI | AO AO)")
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbt_create(t_3c_template, t_3c_der_RI(i_img, i_xyz))
            CALL dbt_create(t_3c_template, t_3c_der_AO(i_img, i_xyz))
         END DO
      END DO

      DO i_mem = 1, n_mem
         CALL build_3c_derivatives(t_3c_der_AO_prv, t_3c_der_RI_prv, ri_data%filter_eps, qs_env, &
                                   nl_3c, basis_set_AO, basis_set_AO, basis_set_RI, &
                                   ri_data%ri_metric, der_eps=ri_data%eps_schwarz_forces, op_pos=2, &
                                   do_kpoints=.TRUE., do_hfx_kpoints=.TRUE., &
                                   bounds_k=[start_blocks(i_mem), end_blocks(i_mem)], &
                                   RI_range=ri_data%kp_RI_range, img_to_RI_cell=ri_data%img_to_RI_cell)

         CALL timeset(routineN//"_cpy", handle2)
         !We go from (mu^0 sigma^i | P^j) to (P^i| sigma^j mu^0) and finally to (P^i| mu^0 sigma^j)
         DO i_img = 1, nimg
            DO i_xyz = 1, 3
               !derivative wrt to mu^0
               CALL get_tensor_occupancy(t_3c_der_AO_prv(i_img, 1, i_xyz), nze, occ)
               IF (nze > 0) THEN
                  CALL dbt_copy(t_3c_der_AO_prv(i_img, 1, i_xyz), t_3c_template, &
                                order=[3, 2, 1], move_data=.TRUE.)
                  CALL dbt_filter(t_3c_template, ri_data%filter_eps)
                  CALL dbt_copy(t_3c_template, t_3c_der_AO(i_img, i_xyz), &
                                order=[1, 3, 2], move_data=.TRUE., summation=.TRUE.)
               END IF

               !derivative wrt to P^i
               CALL get_tensor_occupancy(t_3c_der_RI_prv(i_img, 1, i_xyz), nze, occ)
               IF (nze > 0) THEN
                  CALL dbt_copy(t_3c_der_RI_prv(i_img, 1, i_xyz), t_3c_template, &
                                order=[3, 2, 1], move_data=.TRUE.)
                  CALL dbt_filter(t_3c_template, ri_data%filter_eps)
                  CALL dbt_copy(t_3c_template, t_3c_der_RI(i_img, i_xyz), &
                                order=[1, 3, 2], move_data=.TRUE., summation=.TRUE.)
               END IF
            END DO
         END DO
         CALL timestop(handle2)
      END DO
      CALL dbt_destroy(t_3c_template)

      CALL neighbor_list_3c_destroy(nl_3c)
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbt_destroy(t_3c_der_RI_prv(i_img, 1, i_xyz))
            CALL dbt_destroy(t_3c_der_AO_prv(i_img, 1, i_xyz))
         END DO
      END DO
      DEALLOCATE (t_3c_der_RI_prv, t_3c_der_AO_prv)

      !Reorder 3c derivatives to be consistant with ints
      CALL reorder_3c_derivs(t_3c_der_RI, ri_data)
      CALL reorder_3c_derivs(t_3c_der_AO, ri_data)

      CALL timeset(routineN//"_2c", handle2)
      !The 2-center derivatives
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(ri_data%bsizes_RI)))
      ALLOCATE (col_bsize(SIZE(ri_data%bsizes_RI)))
      row_bsize(:) = ri_data%bsizes_RI
      col_bsize(:) = ri_data%bsizes_RI

      CALL dbcsr_create(dbcsr_template, "2c_der", dbcsr_dist, dbcsr_type_no_symmetry, &
                        row_bsize, col_bsize)
      CALL dbcsr_distribution_release(dbcsr_dist)
      DEALLOCATE (col_bsize, row_bsize)

      ALLOCATE (mat_der_metric(nimg, 3))
      DO i_xyz = 1, 3
         DO i_img = 1, nimg
            CALL dbcsr_create(mat_der_pot(i_img, i_xyz), template=dbcsr_template)
            CALL dbcsr_create(mat_der_metric(i_img, i_xyz), template=dbcsr_template)
         END DO
      END DO
      CALL dbcsr_release(dbcsr_template)

      !HFX potential derivatives
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%hfx_pot, &
                                   "HFX_2c_nl_pot", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)
      CALL build_2c_derivatives(mat_der_pot, ri_data%filter_eps_2c, qs_env, nl_2c, &
                                basis_set_RI, basis_set_RI, ri_data%hfx_pot, do_kpoints=.TRUE.)
      CALL release_neighbor_list_sets(nl_2c)

      !RI metric derivatives
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%ri_metric, &
                                   "HFX_2c_nl_pot", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)
      CALL build_2c_derivatives(mat_der_metric, ri_data%filter_eps_2c, qs_env, nl_2c, &
                                basis_set_RI, basis_set_RI, ri_data%ri_metric, do_kpoints=.TRUE.)
      CALL release_neighbor_list_sets(nl_2c)

      !Get into extended RI basis and tensor format
      DO i_xyz = 1, 3
         DO iatom = 1, natom
            CALL dbt_create(ri_data%t_2c_inv(1, 1), t_2c_der_metric(iatom, i_xyz))
            CALL get_ext_2c_int(t_2c_der_metric(iatom, i_xyz), mat_der_metric(:, i_xyz), &
                                iatom, iatom, 1, ri_data, qs_env)
         END DO
         DO i_img = 1, nimg
            CALL dbcsr_release(mat_der_metric(i_img, i_xyz))
         END DO
      END DO
      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE precalc_derivatives

! **************************************************************************************************
!> \brief Update the forces due to the derivative of the a 2-center product d/dR (Q|R)
!> \param force ...
!> \param t_2c_contr A precontracted tensor containing sum_abcdPS (ab|P)(P|Q)^-1 (R|S)^-1 (S|cd) P_ac P_bd
!> \param t_2c_der the d/dR (Q|R) tensor, in all 3 cartesian directions
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param img in which periodic image the second center of the tensor is
!> \param pref ...
!> \param ri_data ...
!> \param qs_env ...
!> \param work_virial ...
!> \param cell ...
!> \param particle_set ...
!> \param diag ...
!> \param offdiag ...
!> \note IMPORTANT: t_tc_contr and t_2c_der need to have the same distribution. Atomic block sizes are
!>                  assumed
! **************************************************************************************************
   SUBROUTINE get_2c_der_force(force, t_2c_contr, t_2c_der, atom_of_kind, kind_of, img, pref, &
                               ri_data, qs_env, work_virial, cell, particle_set, diag, offdiag)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_contr
      TYPE(dbt_type), DIMENSION(:), INTENT(INOUT)        :: t_2c_der
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of
      INTEGER, INTENT(IN)                                :: img
      REAL(dp), INTENT(IN)                               :: pref
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: work_virial
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      LOGICAL, INTENT(IN), OPTIONAL                      :: diag, offdiag

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_2c_der_force'

      INTEGER                                            :: handle, i_img, i_RI, i_xyz, iat, &
                                                            iat_of_kind, ikind, j_img, j_RI, &
                                                            j_xyz, jat, jat_of_kind, jkind, natom
      INTEGER, DIMENSION(2)                              :: ind
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL                                            :: found, my_diag, my_offdiag, use_virial
      REAL(dp)                                           :: new_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), TARGET     :: contr_blk, der_blk
      REAL(dp), DIMENSION(3)                             :: scoord
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (kpoints, index_to_cell)

      !Loop over the blocks of d/dR (Q|R), contract with the corresponding block of t_2c_contr and
      !update the relevant force

      CALL timeset(routineN, handle)

      use_virial = .FALSE.
      IF (PRESENT(work_virial) .AND. PRESENT(cell) .AND. PRESENT(particle_set)) use_virial = .TRUE.

      my_diag = .FALSE.
      IF (PRESENT(diag)) my_diag = diag

      my_offdiag = .FALSE.
      IF (PRESENT(diag)) my_offdiag = offdiag

      CALL get_qs_env(qs_env, kpoints=kpoints, natom=natom)
      CALL get_kpoint_info(kpoints, index_to_cell=index_to_cell)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_2c_der,t_2c_contr,work_virial,force,use_virial,natom,index_to_cell,ri_data,img) &
!$OMP SHARED(pref,atom_of_kind,kind_of,particle_set,cell,my_diag,my_offdiag) &
!$OMP PRIVATE(i_xyz,j_xyz,iter,ind,der_blk,contr_blk,found,new_force,i_RI,i_img,j_RI,j_img) &
!$OMP PRIVATE(iat,jat,iat_of_kind,jat_of_kind,ikind,jkind,scoord)
      DO i_xyz = 1, 3
         CALL dbt_iterator_start(iter, t_2c_der(i_xyz))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)

            !Only take forecs due to block diagonal or block off-diagonal, depending on arguments
            IF ((my_diag .AND. .NOT. my_offdiag) .OR. (.NOT. my_diag .AND. my_offdiag)) THEN
               IF (my_diag .AND. (ind(1) .NE. ind(2))) CYCLE
               IF (my_offdiag .AND. (ind(1) == ind(2))) CYCLE
            END IF

            CALL dbt_get_block(t_2c_der(i_xyz), ind, der_blk, found)
            CPASSERT(found)
            CALL dbt_get_block(t_2c_contr, ind, contr_blk, found)

            IF (found) THEN

               !an element of d/dR (Q|R) corresponds to 2 things because of translational invariance
               !(Q'| R) = - (Q| R'), once wrt the center on Q, and once on R
               new_force = pref*SUM(der_blk(:, :)*contr_blk(:, :))

               i_RI = (ind(1) - 1)/natom + 1
               i_img = ri_data%RI_cell_to_img(i_RI)
               iat = ind(1) - (i_RI - 1)*natom
               iat_of_kind = atom_of_kind(iat)
               ikind = kind_of(iat)

               j_RI = (ind(2) - 1)/natom + 1
               j_img = ri_data%RI_cell_to_img(j_RI)
               jat = ind(2) - (j_RI - 1)*natom
               jat_of_kind = atom_of_kind(jat)
               jkind = kind_of(jat)

               !Force on iatom (first center)
!$OMP ATOMIC
               force(ikind)%fock_4c(i_xyz, iat_of_kind) = force(ikind)%fock_4c(i_xyz, iat_of_kind) &
                                                          + new_force

               IF (use_virial) THEN

                  CALL real_to_scaled(scoord, pbc(particle_set(iat)%r, cell), cell)
                  scoord(:) = scoord(:) + REAL(index_to_cell(:, i_img), dp)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
                  END DO
               END IF

               !Force on jatom (second center)
!$OMP ATOMIC
               force(jkind)%fock_4c(i_xyz, jat_of_kind) = force(jkind)%fock_4c(i_xyz, jat_of_kind) &
                                                          - new_force

               IF (use_virial) THEN

                  CALL real_to_scaled(scoord, pbc(particle_set(jat)%r, cell), cell)
                  scoord(:) = scoord(:) + REAL(index_to_cell(:, j_img) + index_to_cell(:, img), dp)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) - new_force*scoord(j_xyz)
                  END DO
               END IF

               DEALLOCATE (contr_blk)
            END IF

            DEALLOCATE (der_blk)
         END DO !iter
         CALL dbt_iterator_stop(iter)

      END DO !i_xyz
!$OMP END PARALLEL
      CALL timestop(handle)

   END SUBROUTINE get_2c_der_force

! **************************************************************************************************
!> \brief This routines calculates the force contribution from a trace over 3D tensors, i.e.
!>        force = sum_ijk A_ijk B_ijk., the B tensor is (P^0| sigma^0 lambda^img), with P in the
!>        extended RI basis. Note that all tensors are stacked along the 3rd dimension
!> \param force ...
!> \param t_3c_contr ...
!> \param t_3c_der_1 ...
!> \param t_3c_der_2 ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param idx_to_at_RI ...
!> \param idx_to_at_AO ...
!> \param i_images ...
!> \param lb_img ...
!> \param pref ...
!> \param ri_data ...
!> \param qs_env ...
!> \param work_virial ...
!> \param cell ...
!> \param particle_set ...
! **************************************************************************************************
   SUBROUTINE get_force_from_3c_trace(force, t_3c_contr, t_3c_der_1, t_3c_der_2, atom_of_kind, kind_of, &
                                      idx_to_at_RI, idx_to_at_AO, i_images, lb_img, pref, &
                                      ri_data, qs_env, work_virial, cell, particle_set)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_contr
      TYPE(dbt_type), DIMENSION(3), INTENT(INOUT)        :: t_3c_der_1, t_3c_der_2
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of, idx_to_at_RI, &
                                                            idx_to_at_AO, i_images
      INTEGER, INTENT(IN)                                :: lb_img
      REAL(dp), INTENT(IN)                               :: pref
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: work_virial
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_force_from_3c_trace'

      INTEGER :: handle, i_RI, i_xyz, iat, iat_of_kind, idx, ikind, j_xyz, jat, jat_of_kind, &
         jkind, kat, kat_of_kind, kkind, nblks_AO, nblks_RI, RI_img
      INTEGER, DIMENSION(3)                              :: ind
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL                                            :: found, found_1, found_2, use_virial
      REAL(dp)                                           :: new_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET  :: contr_blk, der_blk_1, der_blk_2, &
                                                            der_blk_3
      REAL(dp), DIMENSION(3)                             :: scoord
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(kpoint_type), POINTER                         :: kpoints

      NULLIFY (kpoints, index_to_cell)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, index_to_cell=index_to_cell)

      nblks_RI = SIZE(ri_data%bsizes_RI_split)
      nblks_AO = SIZE(ri_data%bsizes_AO_split)

      use_virial = .FALSE.
      IF (PRESENT(work_virial) .AND. PRESENT(cell) .AND. PRESENT(particle_set)) use_virial = .TRUE.

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_3c_der_1, t_3c_der_2,t_3c_contr,work_virial,force,use_virial,index_to_cell,i_images,lb_img) &
!$OMP SHARED(pref,idx_to_at_AO,atom_of_kind,kind_of,particle_set,cell,idx_to_at_RI,ri_data,nblks_RI,nblks_AO) &
!$OMP PRIVATE(i_xyz,j_xyz,iter,ind,der_blk_1,contr_blk,found,new_force,iat,iat_of_kind,ikind,scoord) &
!$OMP PRIVATE(jat,kat,jat_of_kind,kat_of_kind,jkind,kkind,i_RI,RI_img,der_blk_2,der_blk_3,found_1,found_2,idx)
      CALL dbt_iterator_start(iter, t_3c_contr)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)

         CALL dbt_get_block(t_3c_contr, ind, contr_blk, found)
         IF (found) THEN

            DO i_xyz = 1, 3
               CALL dbt_get_block(t_3c_der_1(i_xyz), ind, der_blk_1, found_1)
               IF (.NOT. found_1) THEN
                  DEALLOCATE (der_blk_1)
                  ALLOCATE (der_blk_1(SIZE(contr_blk, 1), SIZE(contr_blk, 2), SIZE(contr_blk, 3)))
                  der_blk_1(:, :, :) = 0.0_dp
               END IF
               CALL dbt_get_block(t_3c_der_2(i_xyz), ind, der_blk_2, found_2)
               IF (.NOT. found_2) THEN
                  DEALLOCATE (der_blk_2)
                  ALLOCATE (der_blk_2(SIZE(contr_blk, 1), SIZE(contr_blk, 2), SIZE(contr_blk, 3)))
                  der_blk_2(:, :, :) = 0.0_dp
               END IF

               ALLOCATE (der_blk_3(SIZE(contr_blk, 1), SIZE(contr_blk, 2), SIZE(contr_blk, 3)))
               der_blk_3(:, :, :) = -(der_blk_1(:, :, :) + der_blk_2(:, :, :))

               !We assume the tensors are in the format (P^0| sigma^0 mu^a+c-b), with P a member of the
               !extended RI basis set

               !Force for the first center (RI extended basis, zero cell)
               new_force = pref*SUM(der_blk_1(:, :, :)*contr_blk(:, :, :))

               i_RI = (ind(1) - 1)/nblks_RI + 1
               RI_img = ri_data%RI_cell_to_img(i_RI)
               iat = idx_to_at_RI(ind(1) - (i_RI - 1)*nblks_RI)
               iat_of_kind = atom_of_kind(iat)
               ikind = kind_of(iat)

!$OMP ATOMIC
               force(ikind)%fock_4c(i_xyz, iat_of_kind) = force(ikind)%fock_4c(i_xyz, iat_of_kind) &
                                                          + new_force

               IF (use_virial) THEN

                  CALL real_to_scaled(scoord, pbc(particle_set(iat)%r, cell), cell)
                  scoord(:) = scoord(:) + REAL(index_to_cell(:, RI_img), dp)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
                  END DO
               END IF

               !Force with respect to the second center (AO basis, zero cell)
               new_force = pref*SUM(der_blk_2(:, :, :)*contr_blk(:, :, :))
               jat = idx_to_at_AO(ind(2))
               jat_of_kind = atom_of_kind(jat)
               jkind = kind_of(jat)

!$OMP ATOMIC
               force(jkind)%fock_4c(i_xyz, jat_of_kind) = force(jkind)%fock_4c(i_xyz, jat_of_kind) &
                                                          + new_force

               IF (use_virial) THEN

                  CALL real_to_scaled(scoord, pbc(particle_set(jat)%r, cell), cell)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
                  END DO
               END IF

               !Force with respect to the third center (AO basis, apc_img - b_img)
               !Note: tensors are stacked along the 3rd direction
               new_force = pref*SUM(der_blk_3(:, :, :)*contr_blk(:, :, :))
               idx = (ind(3) - 1)/nblks_AO + 1
               kat = idx_to_at_AO(ind(3) - (idx - 1)*nblks_AO)
               kat_of_kind = atom_of_kind(kat)
               kkind = kind_of(kat)

!$OMP ATOMIC
               force(kkind)%fock_4c(i_xyz, kat_of_kind) = force(kkind)%fock_4c(i_xyz, kat_of_kind) &
                                                          + new_force

               IF (use_virial) THEN
                  CALL real_to_scaled(scoord, pbc(particle_set(kat)%r, cell), cell)
                  scoord(:) = scoord(:) + REAL(index_to_cell(:, i_images(lb_img - 1 + idx)), dp)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
                  END DO
               END IF

               DEALLOCATE (der_blk_1, der_blk_2, der_blk_3)
            END DO !i_xyz
            DEALLOCATE (contr_blk)
         END IF !found
      END DO !iter
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      CALL timestop(handle)

   END SUBROUTINE get_force_from_3c_trace

END MODULE

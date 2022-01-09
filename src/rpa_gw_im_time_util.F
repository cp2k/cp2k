!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for GW with imaginary time
!> \par History
!>      06.2019 split from rpa_im_time.F [Frederick Stein]
! **************************************************************************************************
MODULE rpa_gw_im_time_util

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_row_template
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_get_diag, &
        dbcsr_get_info, dbcsr_get_stored_coordinates, dbcsr_init_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_release_p, dbcsr_reserve_all_blocks, &
        dbcsr_set_diag, dbcsr_type, dbcsr_type_no_symmetry
   USE dbt_api,                         ONLY: &
        dbt_contract, dbt_copy, dbt_copy_matrix_to_tensor, dbt_create, dbt_default_distvec, &
        dbt_destroy, dbt_get_info, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE hfx_types,                       ONLY: alloc_containers,&
                                              block_ind_type,&
                                              hfx_compression_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_dims_create,&
                                              mp_sum
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_tensors,                      ONLY: compress_tensor,&
                                              decompress_tensor,&
                                              get_tensor_occupancy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              pgf_block_sizes,&
                                              split_block_sizes
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_im_time_util'

   PUBLIC :: get_tensor_3c_overl_int_gw, compute_weight_re_im, get_atom_index_from_basis_function_index, &
             compute_shortest_lattice_vector
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param t_3c_overl_int ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param t_3c_overl_int_ao_mo ...
!> \param t_3c_O_mo_compressed ...
!> \param t_3c_O_mo_ind ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param mo_coeff ...
!> \param matrix_s ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param para_env ...
!> \param do_ic_model ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE get_tensor_3c_overl_int_gw(t_3c_overl_int, t_3c_O_compressed, t_3c_O_ind, &
                                         t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                         t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                         starts_array_mc, ends_array_mc, &
                                         mo_coeff, matrix_s, &
                                         gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, &
                                         para_env, &
                                         do_ic_model, &
                                         t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                         qs_env, unit_nr, do_beta)

      TYPE(dbt_type), DIMENSION(:, :)                    :: t_3c_overl_int
      TYPE(hfx_compression_type), DIMENSION(:, :, :)     :: t_3c_O_compressed
      TYPE(block_ind_type), DIMENSION(:, :, :)           :: t_3c_O_ind
      TYPE(dbt_type)                                     :: t_3c_overl_int_ao_mo
      TYPE(hfx_compression_type)                         :: t_3c_O_mo_compressed
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: t_3c_O_mo_ind
      TYPE(dbt_type)                                     :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                                            nmo
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: do_ic_model
      TYPE(dbt_type)                                     :: t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_tensor_3c_overl_int_gw'

      INTEGER :: cut_memory, handle, i_mem, icol_global, imo, irow_global, min_bsize, &
         min_bsize_mo, nkind, nmo_blk_gw, npcols, nprows, size_MO, unit_nr_prv
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2, dist3, sizes_AO, &
                                                            sizes_AO_split, sizes_MO, sizes_MO_1, &
                                                            sizes_RI, sizes_RI_split, tmp
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(2, 1)                           :: bounds
      INTEGER, DIMENSION(2, 3)                           :: ibounds
      INTEGER, DIMENSION(3)                              :: bounds_3c, pdims
      INTEGER, DIMENSION(:), POINTER                     :: distp_1, distp_2, sizes_MO_blocked, &
                                                            sizes_MO_p1, sizes_MO_p2
      LOGICAL                                            :: memory_info, my_do_beta
      REAL(dp)                                           :: compression_factor, memory_3c, occ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: norm
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), POINTER                          :: fm_mat_mo_coeff_gw
      TYPE(dbcsr_distribution_type)                      :: dist, dist_templ
      TYPE(dbcsr_type)                                   :: mat_mo_coeff_gw_reflected_norm, &
                                                            mat_norm, mat_norm_diag, mat_work
      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_gw, &
                                                            mat_mo_coeff_gw_reflected
      TYPE(dbt_pgrid_type)                               :: pgrid_2d, pgrid_AO, pgrid_ic, pgrid_MO
      TYPE(dbt_type)                                     :: mo_coeff_gw_t, mo_coeff_gw_t_tmp, &
                                                            t_3c_overl_int_ao_ao, &
                                                            t_3c_overl_int_mo_ao, &
                                                            t_3c_overl_int_mo_mo
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ao, basis_set_ri_aux
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      memory_info = qs_env%mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_prv = unit_nr
      ELSE
         unit_nr_prv = 0
      END IF

      my_do_beta = .FALSE.
      IF (PRESENT(do_beta)) my_do_beta = do_beta

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, particle_set=particle_set, atomic_kind_set=atomic_kind_set)

      NULLIFY (fm_mat_mo_coeff_gw)
      CALL cp_fm_create(fm_mat_mo_coeff_gw, mo_coeff%matrix_struct)
      CALL cp_fm_to_fm(mo_coeff, fm_mat_mo_coeff_gw)

      ! set MO coeffs to zero where
      DO irow_global = 1, nmo
         DO icol_global = 1, homo - gw_corr_lev_occ
            CALL cp_fm_set_element(fm_mat_mo_coeff_gw, irow_global, icol_global, 0.0_dp)
         END DO
         DO icol_global = homo + gw_corr_lev_virt + 1, nmo
            CALL cp_fm_set_element(fm_mat_mo_coeff_gw, irow_global, icol_global, 0.0_dp)
         END DO
      END DO

      NULLIFY (mat_mo_coeff_gw)
      CALL dbcsr_init_p(mat_mo_coeff_gw)

      CALL cp_dbcsr_m_by_n_from_row_template(mat_mo_coeff_gw, template=matrix_s(1)%matrix, n=nmo, &
                                             sym=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mat_mo_coeff_gw, &
                            mat_mo_coeff_gw, &
                            keep_sparsity=.FALSE.)

      ! just remove the blocks which have been set to zero
      CALL dbcsr_filter(mat_mo_coeff_gw, 1.0E-20_dp)

      min_bsize = qs_env%mp2_env%ri_rpa_im_time%min_bsize
      min_bsize_mo = qs_env%mp2_env%ri_rpa_im_time%min_bsize_mo

      CALL split_block_sizes([gw_corr_lev_occ + gw_corr_lev_virt], sizes_MO, min_bsize_mo)
      ALLOCATE (sizes_MO_1(nmo))
      sizes_MO_1(:) = 1

      nmo_blk_gw = SIZE(sizes_MO)
      CALL move_alloc(sizes_MO, tmp)
      ALLOCATE (sizes_MO(nmo_blk_gw + 2))
      sizes_MO(1) = homo - gw_corr_lev_occ
      sizes_MO(2:SIZE(tmp) + 1) = tmp(:)
      sizes_MO(SIZE(tmp) + 2) = nmo - (homo + gw_corr_lev_virt)

      ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
      CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_ri_aux)
      CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_ao)

      CALL pgf_block_sizes(atomic_kind_set, basis_set_ao, min_bsize, sizes_AO_split)
      CALL pgf_block_sizes(atomic_kind_set, basis_set_ri_aux, min_bsize, sizes_RI_split)

      DEALLOCATE (basis_set_ao, basis_set_ri_aux)

      pdims = 0
      CALL dbt_pgrid_create(para_env%group, pdims, pgrid_AO, &
                            tensor_dims=[SIZE(sizes_RI_split), SIZE(sizes_AO_split), SIZE(sizes_AO_split)])

      pdims_2d = 0
      CALL mp_dims_create(para_env%num_pe, pdims_2d)

      ! we iterate over MO blocks for saving memory during contraction, thus we should not parallelize over MO dimension
      pdims = [pdims_2d(1), pdims_2d(2), 1]
      CALL dbt_pgrid_create(para_env%group, pdims, pgrid_MO, &
                            tensor_dims=[SIZE(sizes_RI_split), SIZE(sizes_AO_split), 1])

      pdims_2d = 0
      CALL dbt_pgrid_create(para_env%group, pdims_2d, pgrid_2d, &
                            tensor_dims=[SIZE(sizes_AO_split), nmo])

      CALL create_3c_tensor(t_3c_overl_int_ao_ao, dist1, dist2, dist3, pgrid_AO, &
                            sizes_RI_split, sizes_AO_split, sizes_AO_split, [1, 2], [3], name="(RI AO | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      MARK_USED(t_3c_O_compressed)
      MARK_USED(t_3c_O_ind)

      CALL create_3c_tensor(t_3c_overl_int_ao_mo, dist1, dist2, dist3, pgrid_AO, &
                            sizes_RI_split, sizes_AO_split, sizes_MO_1, [1, 2], [3], name="(RI AO | MO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(t_3c_overl_int_gw_RI, dist1, dist2, dist3, pgrid_MO, &
                            sizes_RI_split, sizes_AO_split, sizes_MO, [1], [2, 3], name="(RI | AO MO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(t_3c_overl_int_gw_AO, dist1, dist2, dist3, pgrid_MO, &
                            sizes_AO_split, sizes_RI_split, sizes_MO, [1], [2, 3], name="(AO | RI MO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL dbt_pgrid_destroy(pgrid_AO)
      CALL dbt_pgrid_destroy(pgrid_MO)

      CALL create_2c_tensor(mo_coeff_gw_t, dist1, dist2, pgrid_2d, sizes_AO_split, sizes_MO_1, name="(AO|MO)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_pgrid_destroy(pgrid_2d)

      CALL dbt_create(mat_mo_coeff_gw, mo_coeff_gw_t_tmp, name="MO coeffs")
      CALL dbt_copy_matrix_to_tensor(mat_mo_coeff_gw, mo_coeff_gw_t_tmp)

      CALL dbt_copy(mo_coeff_gw_t_tmp, mo_coeff_gw_t)

      bounds(1, 1) = homo - gw_corr_lev_occ + 1
      bounds(2, 1) = homo + gw_corr_lev_virt

      CALL dbt_get_info(t_3c_overl_int_ao_ao, nfull_total=bounds_3c)

      ibounds(:, 1) = [1, bounds_3c(1)]
      ibounds(:, 3) = [1, bounds_3c(3)]

      cut_memory = SIZE(starts_array_mc)

      DO i_mem = 1, cut_memory
         CALL decompress_tensor(t_3c_overl_int(1, 1), t_3c_O_ind(1, 1, i_mem)%ind, t_3c_O_compressed(1, 1, i_mem), &
                                qs_env%mp2_env%ri_rpa_im_time%eps_compress)

         ibounds(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

         CALL dbt_copy(t_3c_overl_int(1, 1), t_3c_overl_int_ao_ao, move_data=.TRUE.)

         CALL dbt_contract(1.0_dp, mo_coeff_gw_t, t_3c_overl_int_ao_ao, 1.0_dp, &
                           t_3c_overl_int_ao_mo, contract_1=[1], notcontract_1=[2], &
                           contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                           bounds_2=ibounds, move_data=.FALSE., unit_nr=unit_nr_prv)

      END DO

      CALL cp_fm_release(fm_mat_mo_coeff_gw)

      IF (do_ic_model) THEN
         pdims = 0
         CALL dbt_pgrid_create(para_env%group, pdims, pgrid_ic, &
                               tensor_dims=[SIZE(sizes_RI_split), nmo, nmo])

         CALL create_3c_tensor(t_3c_overl_int_mo_ao, dist1, dist2, dist3, pgrid_ic, &
                               sizes_RI_split, sizes_MO_1, sizes_AO_split, [1, 2], [3], name="(RI MO | AO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL create_3c_tensor(t_3c_overl_int_mo_mo, dist1, dist2, dist3, pgrid_ic, &
                               sizes_RI_split, sizes_MO_1, sizes_MO_1, [1, 2], [3], name="(RI MO | MO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL dbt_create(t_3c_overl_int_mo_mo, t_3c_overl_nnP_ic)
         CALL create_3c_tensor(t_3c_overl_nnP_ic_reflected, dist1, dist2, dist3, pgrid_ic, &
                               sizes_RI_split, sizes_MO_1, sizes_MO_1, [1], [2, 3], name="(RI | MO MO)")
         DEALLOCATE (dist1, dist2, dist3)

         CALL dbt_pgrid_destroy(pgrid_ic)

         CALL dbt_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_mo_ao, order=[1, 3, 2])
         CALL dbt_contract(1.0_dp, mo_coeff_gw_t, t_3c_overl_int_mo_ao, 0.0_dp, &
                           t_3c_overl_int_mo_mo, contract_1=[1], notcontract_1=[2], &
                           contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                           bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)
         CALL dbt_copy(t_3c_overl_int_mo_mo, t_3c_overl_nnP_ic)

         NULLIFY (mat_mo_coeff_gw_reflected)
         CALL dbcsr_init_p(mat_mo_coeff_gw_reflected)

         CALL cp_dbcsr_m_by_n_from_row_template(mat_mo_coeff_gw_reflected, template=matrix_s(1)%matrix, n=nmo, &
                                                sym=dbcsr_type_no_symmetry)

         CALL reflect_mat_row(mat_mo_coeff_gw_reflected, mat_mo_coeff_gw, para_env, qs_env, unit_nr, do_beta=my_do_beta)

         ! normalize reflected MOs (they are not properly normalized since high angular momentum basis functions
         ! of the image molecule are not exactly reflected at the image plane (sign problem in p_z function)
         CALL dbcsr_create(matrix=mat_work, template=mat_mo_coeff_gw_reflected, matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_get_info(mat_work, distribution=dist_templ, nblkcols_total=size_MO, col_blk_size=sizes_MO_blocked)

         CALL dbcsr_distribution_get(dist_templ, nprows=nprows, npcols=npcols)

         ALLOCATE (distp_1(size_MO), distp_2(size_MO))
         CALL dbt_default_distvec(size_MO, nprows, sizes_MO_blocked, distp_1)
         CALL dbt_default_distvec(size_MO, npcols, sizes_MO_blocked, distp_2)
         CALL dbcsr_distribution_new(dist, template=dist_templ, row_dist=distp_1, col_dist=distp_2, reuse_arrays=.TRUE.)

         ALLOCATE (sizes_MO_p1(size_MO))
         ALLOCATE (sizes_MO_p2(size_MO))
         sizes_MO_p1(:) = sizes_MO_blocked
         sizes_MO_p2(:) = sizes_MO_blocked
         CALL dbcsr_create(mat_norm, "mo norm", dist, dbcsr_type_no_symmetry, sizes_MO_p1, sizes_MO_p2, &
                           reuse_arrays=.TRUE.)
         CALL dbcsr_distribution_release(dist)

         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, mat_mo_coeff_gw_reflected, 0.0_dp, mat_work)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mat_mo_coeff_gw_reflected, mat_work, 0.0_dp, mat_norm)

         CALL dbcsr_release(mat_work)

         ALLOCATE (norm(nmo))
         norm = 0.0_dp

         CALL dbcsr_get_diag(mat_norm, norm)
         CALL mp_sum(norm, para_env%group)

         DO imo = bounds(1, 1), bounds(2, 1)
            norm(imo) = 1.0_dp/SQRT(norm(imo))
         END DO

         CALL dbcsr_create(mat_norm_diag, template=mat_norm)
         CALL dbcsr_release(mat_norm)

         CALL dbcsr_add_on_diag(mat_norm_diag, 1.0_dp)

         CALL dbcsr_set_diag(mat_norm_diag, norm)

         CALL dbcsr_create(mat_mo_coeff_gw_reflected_norm, template=mat_mo_coeff_gw_reflected)
         CALL dbcsr_multiply("N", "N", 1.0_dp, mat_mo_coeff_gw_reflected, mat_norm_diag, 0.0_dp, mat_mo_coeff_gw_reflected_norm)
         CALL dbcsr_release(mat_norm_diag)

         CALL dbcsr_filter(mat_mo_coeff_gw_reflected_norm, 1.0E-20_dp)

         CALL dbt_copy_matrix_to_tensor(mat_mo_coeff_gw_reflected_norm, mo_coeff_gw_t_tmp)
         CALL dbcsr_release(mat_mo_coeff_gw_reflected_norm)
         CALL dbt_copy(mo_coeff_gw_t_tmp, mo_coeff_gw_t)

         CALL dbt_contract(1.0_dp, mo_coeff_gw_t, t_3c_overl_int_ao_ao, 0.0_dp, &
                           t_3c_overl_int_ao_mo, contract_1=[1], notcontract_1=[2], &
                           contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                           bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)

         CALL dbt_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_mo_ao, order=[1, 3, 2])
         CALL dbt_contract(1.0_dp, mo_coeff_gw_t, t_3c_overl_int_mo_ao, 0.0_dp, &
                           t_3c_overl_int_mo_mo, contract_1=[1], notcontract_1=[2], &
                           contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                           bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)
         CALL dbt_copy(t_3c_overl_int_mo_mo, t_3c_overl_nnP_ic_reflected)
         CALL dbt_destroy(t_3c_overl_int_mo_ao)
         CALL dbt_destroy(t_3c_overl_int_mo_mo)

         CALL dbcsr_release_p(mat_mo_coeff_gw_reflected)

      END IF

      CALL alloc_containers(t_3c_O_mo_compressed, 1)
      CALL get_tensor_occupancy(t_3c_overl_int_ao_mo, nze, occ)
      memory_3c = 0.0_dp

      CALL compress_tensor(t_3c_overl_int_ao_mo, t_3c_O_mo_ind, t_3c_O_mo_compressed, &
                           qs_env%mp2_env%ri_rpa_im_time%eps_compress, memory_3c)

      CALL mp_sum(memory_3c, para_env%group)
      compression_factor = REAL(nze, dp)*1.0E-06*8.0_dp/memory_3c

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="((T3,A,T66,F11.2,A4))") &
            "MEMORY_INFO| Memory of MO-contracted tensor (compressed):", memory_3c, ' MiB'

         WRITE (UNIT=unit_nr, FMT="((T3,A,T60,F21.2))") &
            "MEMORY_INFO| Compression factor:                  ", compression_factor
      END IF

      CALL dbcsr_release_p(mat_mo_coeff_gw)

      CALL dbt_destroy(t_3c_overl_int_ao_ao)
      CALL dbt_destroy(mo_coeff_gw_t)
      CALL dbt_destroy(mo_coeff_gw_t_tmp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief reflect from V = (A,B|B,A) to V_reflected = (B,A|A,B) where A belongs to the block of the molecule
!>        and B to the off diagonal block between molecule and image of the molecule
!> \param mat_reflected ...
!> \param mat_orig ...
!> \param para_env ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE reflect_mat_row(mat_reflected, mat_orig, para_env, qs_env, unit_nr, do_beta)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_reflected
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_orig
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: do_beta

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reflect_mat_row'

      INTEGER :: block, block_size, col, col_rec, col_size, handle, i_atom, i_block, imepos, &
         j_atom, natom, nblkcols_total, nblkrows_total, offset, row, row_rec, row_reflected, &
         row_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, entry_counter, image_atom, &
         num_blocks_rec, num_blocks_send, num_entries_rec, num_entries_send, sizes_rec, sizes_send
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL                                            :: found_image_atom
      REAL(KIND=dp)                                      :: avg_z_dist, delta, eps_dist2, &
                                                            min_z_dist, ra(3), rb(3), sum_z, &
                                                            z_reflection
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL dbcsr_reserve_all_blocks(mat_reflected)

      CALL get_qs_env(qs_env, cell=cell, &
                      particle_set=particle_set)

      ! first check, whether we have an image molecule
      CALL dbcsr_get_info(mat_orig, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          row_blk_size=row_blk_sizes, &
                          col_blk_size=col_blk_sizes)

      natom = SIZE(particle_set)
      CPASSERT(natom == nblkrows_total)

      eps_dist2 = qs_env%mp2_env%ri_g0w0%eps_dist
      eps_dist2 = eps_dist2*eps_dist2

      sum_z = 0.0_dp

      DO i_atom = 1, natom

         ra(:) = pbc(particle_set(i_atom)%r, cell)

         sum_z = sum_z + ra(3)

      END DO

      z_reflection = sum_z/REAL(natom, KIND=dp)

      sum_z = 0.0_dp

      DO i_atom = 1, natom

         ra(:) = pbc(particle_set(i_atom)%r, cell)

         sum_z = sum_z + ABS(ra(3) - z_reflection)

      END DO

      avg_z_dist = sum_z/REAL(natom, KIND=dp)

      min_z_dist = avg_z_dist

      DO i_atom = 1, natom

         ra(:) = pbc(particle_set(i_atom)%r, cell)

         IF (ABS(ra(3) - z_reflection) < min_z_dist) THEN
            min_z_dist = ABS(ra(3) - z_reflection)
         END IF

      END DO

      IF (unit_nr > 0 .AND. .NOT. do_beta) THEN
         WRITE (unit_nr, '(T3,A,T70,F9.2,A2)') 'IC_MODEL| Average distance of the molecule to the image plane:', &
            avg_z_dist*0.529_dp, ' A'
         WRITE (unit_nr, '(T3,A,T70,F9.2,A2)') 'IC_MODEL| Minimum distance of the molecule to the image plane:', &
            min_z_dist*0.529_dp, ' A'
      END IF

      ALLOCATE (image_atom(nblkrows_total))
      image_atom = 0

      DO i_atom = 1, natom

         found_image_atom = .FALSE.

         ra(:) = pbc(particle_set(i_atom)%r, cell)

         DO j_atom = 1, natom

            rb(:) = pbc(particle_set(j_atom)%r, cell)

            delta = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) + rb(3) - 2.0_dp*z_reflection)**2

            ! SQRT(delta) < eps_dist
            IF (delta < eps_dist2) THEN
               ! this CPASSERT ensures that there is at most one image atom for each atom
               CPASSERT(.NOT. found_image_atom)
               image_atom(i_atom) = j_atom
               found_image_atom = .TRUE.
               ! check whether we have the same basis at the image atom
               ! if this is wrong, check whether you have the same basis sets for the molecule and the image
               CPASSERT(row_blk_sizes(i_atom) == row_blk_sizes(j_atom))
            END IF

         END DO

         ! this CPASSERT ensures that there is at least one image atom for each atom
         CPASSERT(found_image_atom)

      END DO

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      ALLOCATE (num_blocks_send(0:para_env%num_pe - 1))
      num_entries_rec = 0
      num_blocks_rec = 0
      num_entries_send = 0
      num_blocks_send = 0

      CALL dbcsr_iterator_start(iter, mat_orig)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         row_reflected = image_atom(row)

         CALL dbcsr_get_stored_coordinates(mat_reflected, row_reflected, col, imepos)

         num_entries_send(imepos) = num_entries_send(imepos) + row_size*col_size
         num_blocks_send(imepos) = num_blocks_send(imepos) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe - 1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe - 1))

         DO imepos = 0, para_env%num_pe - 1

            sizes_send(2*imepos) = num_entries_send(imepos)
            sizes_send(2*imepos + 1) = num_blocks_send(imepos)

         END DO

         CALL mp_alltoall(sizes_send, sizes_rec, 2, para_env%group)

         DO imepos = 0, para_env%num_pe - 1
            num_entries_rec(imepos) = sizes_rec(2*imepos)
            num_blocks_rec(imepos) = sizes_rec(2*imepos + 1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send(0)
         num_blocks_rec(0) = num_blocks_send(0)

      END IF

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send(imepos), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = 0

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      CALL dbcsr_iterator_start(iter, mat_orig)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         row_reflected = image_atom(row)

         CALL dbcsr_get_stored_coordinates(mat_reflected, row_reflected, col, imepos)

         block_size = row_size*col_size

         offset = entry_counter(imepos)

         buffer_send(imepos)%msg(offset + 1:offset + block_size) = &
            RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

         block = block_counter(imepos) + 1

         buffer_send(imepos)%indx(block, 1) = row_reflected
         buffer_send(imepos)%indx(block, 2) = col
         buffer_send(imepos)%indx(block, 3) = offset

         entry_counter(imepos) = entry_counter(imepos) + block_size

         block_counter(imepos) = block_counter(imepos) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array)

      ! fill the reflected matrix
      DO imepos = 0, para_env%num_pe - 1

         DO i_block = 1, num_blocks_rec(imepos)

            row_rec = buffer_rec(imepos)%indx(i_block, 1)
            col_rec = buffer_rec(imepos)%indx(i_block, 2)

            CALL dbcsr_iterator_start(iter, mat_reflected)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                              row_size=row_size, col_size=col_size)

               IF (row_rec == row .AND. col_rec == col) THEN

                  offset = buffer_rec(imepos)%indx(i_block, 3)

                  data_block(:, :) = RESHAPE(buffer_rec(imepos)%msg(offset + 1:offset + row_size*col_size), &
                                             (/row_size, col_size/))

               END IF

            END DO

            CALL dbcsr_iterator_stop(iter)

         END DO

      END DO

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (block_counter, entry_counter)
      DEALLOCATE (num_entries_rec)
      DEALLOCATE (num_blocks_rec)
      DEALLOCATE (num_entries_send)
      DEALLOCATE (num_blocks_send)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param atom_from_basis_index ...
!> \param basis_size ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE get_atom_index_from_basis_function_index(qs_env, atom_from_basis_index, basis_size, basis_type)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_basis_index
      INTEGER                                            :: basis_size
      CHARACTER(LEN=*)                                   :: basis_type

      INTEGER                                            :: iatom, LLL, natom, nkind
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_end, row_blk_start
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, particle_set)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, natom=natom, nkind=nkind, &
                      particle_set=particle_set)

      ALLOCATE (row_blk_start(natom))
      ALLOCATE (row_blk_end(natom))
      ALLOCATE (basis_set(nkind))
      CALL basis_set_list_setup(basis_set, basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=row_blk_start, last_sgf=row_blk_end, &
                            basis=basis_set)
      DO LLL = 1, basis_size
         DO iatom = 1, natom
            IF (LLL >= row_blk_start(iatom) .AND. LLL <= row_blk_end(iatom)) THEN
               atom_from_basis_index(LLL) = iatom
            END IF
         END DO
      END DO

      DEALLOCATE (basis_set)
      DEALLOCATE (row_blk_start)
      DEALLOCATE (row_blk_end)

   END SUBROUTINE get_atom_index_from_basis_function_index

! **************************************************************************************************
!> \brief ...
!> \param weight_re ...
!> \param weight_im ...
!> \param num_cells ...
!> \param iatom ...
!> \param jatom ...
!> \param xkp ...
!> \param wkp_W ...
!> \param cell ...
!> \param index_to_cell ...
!> \param hmat ...
!> \param particle_set ...
!> \param abs_cutoffs_chi_W ...
! **************************************************************************************************
   SUBROUTINE compute_weight_re_im(weight_re, weight_im, &
                                   num_cells, iatom, jatom, xkp, wkp_W, &
                                   cell, index_to_cell, hmat, particle_set, abs_cutoffs_chi_W)

      REAL(KIND=dp)                                      :: weight_re, weight_im
      INTEGER                                            :: num_cells, iatom, jatom
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      REAL(KIND=dp)                                      :: wkp_W
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(2)                        :: abs_cutoffs_chi_W

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_weight_re_im'

      INTEGER                                            :: handle, icell, xcell, ycell, zcell
      REAL(KIND=dp)                                      :: abs_rab_cell, arg, coskl, &
                                                            dampening_function, first_cutoff, &
                                                            second_cutoff, sinkl
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, rab_cell_i

      CALL timeset(routineN, handle)

      first_cutoff = abs_cutoffs_chi_W(1)
      second_cutoff = abs_cutoffs_chi_W(2)

      weight_re = 0.0_dp
      weight_im = 0.0_dp

      DO icell = 1, num_cells

         xcell = index_to_cell(1, icell)
         ycell = index_to_cell(2, icell)
         zcell = index_to_cell(3, icell)

         arg = REAL(xcell, dp)*xkp(1) + REAL(ycell, dp)*xkp(2) + REAL(zcell, dp)*xkp(3)

         coskl = wkp_W*COS(twopi*arg)
         sinkl = wkp_W*SIN(twopi*arg)

         cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, icell), dp))

         rab_cell_i(1:3) = pbc(particle_set(iatom)%r(1:3), cell) - &
                           (pbc(particle_set(jatom)%r(1:3), cell) + cell_vector(1:3))

         abs_rab_cell = SQRT(rab_cell_i(1)**2 + rab_cell_i(2)**2 + rab_cell_i(3)**2)

         ! smooth function that is only one if atoms are very close and becomes smoothly zero
         ! for larger distances
         IF (abs_rab_cell < first_cutoff) THEN
            weight_re = 1.0_dp*coskl
            weight_im = 1.0_dp*sinkl
         ELSE IF (abs_rab_cell < second_cutoff) THEN
            dampening_function = 0.5_dp*(COS(pi*(abs_rab_cell - first_cutoff)/ &
                                             (second_cutoff - first_cutoff)) + 1.0_dp)
            weight_re = dampening_function*coskl
            weight_im = dampening_function*sinkl
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_weight_re_im

! **************************************************************************************************
!> \brief ...
!> \param shortest_lattice_vector ...
!> \param hmat ...
!> \param index_to_cell ...
!> \param num_cells ...
! **************************************************************************************************
   SUBROUTINE compute_shortest_lattice_vector(shortest_lattice_vector, hmat, index_to_cell, num_cells)

      REAL(KIND=dp)                                      :: shortest_lattice_vector
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER                                            :: num_cells

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_shortest_lattice_vector'

      INTEGER                                            :: handle, icell, xcell, ycell, zcell
      REAL(KIND=dp)                                      :: abs_cell_vector
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector

      CALL timeset(routineN, handle)

      shortest_lattice_vector = 1.0E9_dp

      DO icell = 1, num_cells

         xcell = index_to_cell(1, icell)
         ycell = index_to_cell(2, icell)
         zcell = index_to_cell(3, icell)

         IF (xcell == 0 .AND. ycell == 0 .AND. zcell == 0) CYCLE

         cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, icell), dp))

         abs_cell_vector = SQRT(cell_vector(1)**2 + cell_vector(2)**2 + cell_vector(3)**2)

         IF (abs_cell_vector < shortest_lattice_vector) shortest_lattice_vector = abs_cell_vector

      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_shortest_lattice_vector

END MODULE rpa_gw_im_time_util

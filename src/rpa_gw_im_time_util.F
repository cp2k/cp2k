!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for GW with imaginary time
!> \par History
!>      06.2019 split from rpa_im_time.F [Frederick Stein]
! **************************************************************************************************
MODULE rpa_gw_im_time_util

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
        dbcsr_scalar, dbcsr_set_diag, dbcsr_type, dbcsr_type_no_symmetry
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, dbcsr_t_create, &
        dbcsr_t_destroy, dbcsr_t_get_info, dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, &
        dbcsr_t_pgrid_type, dbcsr_t_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_sum
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              cyclic_tensor_dist,&
                                              split_block_sizes
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_im_time_util'

   PUBLIC :: get_tensor_3c_overl_int_gw
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param t_3c_overl_int ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
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
   SUBROUTINE get_tensor_3c_overl_int_gw(t_3c_overl_int, &
                                         t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                         mo_coeff, matrix_s, &
                                         gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, &
                                         para_env, &
                                         do_ic_model, &
                                         t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                         qs_env, unit_nr, do_beta)

      TYPE(dbcsr_t_type), DIMENSION(:, :)                :: t_3c_overl_int
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                                            nmo
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: do_ic_model
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_tensor_3c_overl_int_gw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icol_global, imo, irow_global, &
                                                            npcols, nprows, size_MO, &
                                                            sqrt_max_bsize, unit_nr_prv
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2, dist3, sizes_AO, &
                                                            sizes_AO_split, sizes_MO, sizes_RI, &
                                                            sizes_RI_split
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(2, 1)                           :: bounds
      INTEGER, DIMENSION(3)                              :: pdims
      INTEGER, DIMENSION(:), POINTER                     :: distp_1, distp_2, sizes_MO_blocked, &
                                                            sizes_MO_p1, sizes_MO_p2
      LOGICAL                                            :: memory_info, my_do_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: norm
      TYPE(cp_fm_type), POINTER                          :: fm_mat_mo_coeff_gw
      TYPE(dbcsr_distribution_type)                      :: dist, dist_templ
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid, pgrid_2d, pgrid_ic
      TYPE(dbcsr_t_type) :: mo_coeff_gw_t, mo_coeff_gw_t_tmp, t_3c_overl_int_ao_ao, &
         t_3c_overl_int_ao_mo, t_3c_overl_int_mo_ao, t_3c_overl_int_mo_mo
      TYPE(dbcsr_type)                                   :: mat_mo_coeff_gw_reflected_norm, &
                                                            mat_norm, mat_norm_diag, mat_work
      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_gw, &
                                                            mat_mo_coeff_gw_reflected

      memory_info = qs_env%mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_prv = unit_nr
      ELSE
         unit_nr_prv = 0
      ENDIF

      my_do_beta = .FALSE.
      IF (PRESENT(do_beta)) my_do_beta = do_beta

      CALL timeset(routineN, handle)

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

      sqrt_max_bsize = qs_env%mp2_env%ri_rpa_im_time%max_bsize_sqrt

      CALL dbcsr_t_get_info(t_3c_overl_int(1, 1), blk_size_1=sizes_RI, blk_size_2=sizes_AO)
      CALL split_block_sizes(sizes_AO, sizes_AO_split, sqrt_max_bsize)
      CALL split_block_sizes(sizes_RI, sizes_RI_split, sqrt_max_bsize)
      pdims = 0

      ALLOCATE (sizes_MO(nmo))
      sizes_MO = 1

      CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid, &
                                tensor_dims=[SIZE(sizes_RI_split), SIZE(sizes_AO_split), SIZE(sizes_MO)])

      pdims_2d = 0
      CALL dbcsr_t_pgrid_create(para_env%group, pdims_2d, pgrid_2d, &
                                tensor_dims=[SIZE(sizes_AO), SIZE(sizes_MO)])

      CALL create_3c_tensor(t_3c_overl_int_ao_ao, dist1, dist2, dist3, pgrid, &
                            sizes_RI_split, sizes_AO_split, sizes_AO, [1, 2], [3], "(RI AO | AO)")
      DEALLOCATE (dist1, dist2, dist3)
      CALL dbcsr_t_copy(t_3c_overl_int(1, 1), t_3c_overl_int_ao_ao)
      CALL create_3c_tensor(t_3c_overl_int_ao_mo, dist1, dist2, dist3, pgrid, &
                            sizes_RI_split, sizes_AO_split, sizes_MO, [1, 2], [3], "(RI AO | MO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(t_3c_overl_int_gw_RI, dist1, dist2, dist3, pgrid, &
                            sizes_RI, sizes_AO_split, sizes_MO, [1], [2, 3], "(RI | AO MO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(t_3c_overl_int_gw_AO, dist1, dist2, dist3, pgrid, &
                            sizes_AO, sizes_RI_split, sizes_MO, [1], [2, 3], "(AO | RI MO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL dbcsr_t_pgrid_destroy(pgrid)

      CALL create_2c_tensor(mo_coeff_gw_t, dist1, dist2, pgrid_2d, sizes_AO, sizes_MO, name="(AO|MO)")
      DEALLOCATE (dist1, dist2)
      CALL dbcsr_t_pgrid_destroy(pgrid_2d)

      CALL dbcsr_t_create(mat_mo_coeff_gw, mo_coeff_gw_t_tmp, name="MO coeffs")
      CALL dbcsr_t_copy_matrix_to_tensor(mat_mo_coeff_gw, mo_coeff_gw_t_tmp)

      CALL dbcsr_t_copy(mo_coeff_gw_t_tmp, mo_coeff_gw_t)

      bounds(1, 1) = homo - gw_corr_lev_occ + 1
      bounds(2, 1) = homo + gw_corr_lev_virt

      CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_gw_t, t_3c_overl_int_ao_ao, dbcsr_scalar(0.0_dp), &
                            t_3c_overl_int_ao_mo, contract_1=[1], notcontract_1=[2], &
                            contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                            bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)

      CALL dbcsr_t_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_gw_RI)
      CALL dbcsr_t_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_gw_AO, order=[2, 1, 3])

      CALL cp_fm_release(fm_mat_mo_coeff_gw)

      IF (do_ic_model) THEN
         pdims = 0
         CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid_ic, &
                                   tensor_dims=[SIZE(sizes_RI_split), SIZE(sizes_MO), SIZE(sizes_MO)])

         CALL create_3c_tensor(t_3c_overl_int_mo_ao, dist1, dist2, dist3, pgrid_ic, &
                               sizes_RI_split, sizes_MO, sizes_AO, [1, 2], [3], "(RI MO | AO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL create_3c_tensor(t_3c_overl_int_mo_mo, dist1, dist2, dist3, pgrid_ic, &
                               sizes_RI_split, sizes_MO, sizes_MO, [1, 2], [3], "(RI MO | MO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL dbcsr_t_create(t_3c_overl_int_mo_mo, t_3c_overl_nnP_ic)
         CALL create_3c_tensor(t_3c_overl_nnP_ic_reflected, dist1, dist2, dist3, pgrid_ic, &
                               sizes_RI_split, sizes_MO, sizes_MO, [1], [2, 3], "(RI | MO MO)")
         DEALLOCATE (dist1, dist2, dist3)

         CALL dbcsr_t_pgrid_destroy(pgrid_ic)

         CALL dbcsr_t_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_mo_ao, order=[1, 3, 2])
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_gw_t, t_3c_overl_int_mo_ao, dbcsr_scalar(0.0_dp), &
                               t_3c_overl_int_mo_mo, contract_1=[1], notcontract_1=[2], &
                               contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                               bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)
         CALL dbcsr_t_copy(t_3c_overl_int_mo_mo, t_3c_overl_nnP_ic)

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
         CALL cyclic_tensor_dist(size_MO, nprows, sizes_MO_blocked, distp_1)
         CALL cyclic_tensor_dist(size_MO, npcols, sizes_MO_blocked, distp_2)
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
         ENDDO

         CALL dbcsr_create(mat_norm_diag, template=mat_norm)
         CALL dbcsr_release(mat_norm)

         CALL dbcsr_add_on_diag(mat_norm_diag, 1.0_dp)

         CALL dbcsr_set_diag(mat_norm_diag, norm)

         CALL dbcsr_create(mat_mo_coeff_gw_reflected_norm, template=mat_mo_coeff_gw_reflected)
         CALL dbcsr_multiply("N", "N", 1.0_dp, mat_mo_coeff_gw_reflected, mat_norm_diag, 0.0_dp, mat_mo_coeff_gw_reflected_norm)
         CALL dbcsr_release(mat_norm_diag)

         CALL dbcsr_filter(mat_mo_coeff_gw_reflected_norm, 1.0E-20_dp)

         CALL dbcsr_t_copy_matrix_to_tensor(mat_mo_coeff_gw_reflected_norm, mo_coeff_gw_t_tmp)
         CALL dbcsr_release(mat_mo_coeff_gw_reflected_norm)
         CALL dbcsr_t_copy(mo_coeff_gw_t_tmp, mo_coeff_gw_t)

         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_gw_t, t_3c_overl_int_ao_ao, dbcsr_scalar(0.0_dp), &
                               t_3c_overl_int_ao_mo, contract_1=[1], notcontract_1=[2], &
                               contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                               bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)

         CALL dbcsr_t_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_mo_ao, order=[1, 3, 2])
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_gw_t, t_3c_overl_int_mo_ao, dbcsr_scalar(0.0_dp), &
                               t_3c_overl_int_mo_mo, contract_1=[1], notcontract_1=[2], &
                               contract_2=[3], notcontract_2=[1, 2], map_1=[3], map_2=[1, 2], &
                               bounds_2=bounds, move_data=.FALSE., unit_nr=unit_nr_prv)
         CALL dbcsr_t_copy(t_3c_overl_int_mo_mo, t_3c_overl_nnP_ic_reflected)
         CALL dbcsr_t_destroy(t_3c_overl_int_mo_ao)
         CALL dbcsr_t_destroy(t_3c_overl_int_mo_mo)

         CALL dbcsr_release_p(mat_mo_coeff_gw_reflected)

      ENDIF

      CALL dbcsr_release_p(mat_mo_coeff_gw)

      CALL dbcsr_t_destroy(t_3c_overl_int_ao_ao)
      CALL dbcsr_t_destroy(mo_coeff_gw_t)
      CALL dbcsr_t_destroy(mo_coeff_gw_t_tmp)
      CALL dbcsr_t_destroy(t_3c_overl_int_ao_mo)

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reflect_mat_row', &
         routineP = moduleN//':'//routineN

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

END MODULE rpa_gw_im_time_util

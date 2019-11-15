!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_finalize, dbcsr_get_diag, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_init_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release_p, dbcsr_reserve_all_blocks, &
        dbcsr_reserve_blocks, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_sum
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_im_time_util'

   PUBLIC :: get_mat_3c_overl_int_gw, &
             replicate_mat_to_subgroup_simple, &
             fill_mat_3c_overl_int_gw
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int ...
!> \param mat_3c_overl_int_gw ...
!> \param mo_coeff ...
!> \param matrix_s ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param mat_dm_virt_local ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param cut_RI ...
!> \param row_from_LLL ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param do_ic_model ...
!> \param do_ic_opt_homo_lumo ...
!> \param mat_3c_overl_nnP_ic ...
!> \param mat_3c_overl_nnP_ic_reflected ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE get_mat_3c_overl_int_gw(mat_3c_overl_int, mat_3c_overl_int_gw, mo_coeff, matrix_s, &
                                      gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, mat_dm_virt_local, &
                                      para_env, para_env_sub, cut_RI, row_from_LLL, &
                                      my_group_L_starts_im_time, my_group_L_sizes_im_time, do_ic_model, &
                                      do_ic_opt_homo_lumo, mat_3c_overl_nnP_ic, mat_3c_overl_nnP_ic_reflected, &
                                      qs_env, unit_nr, do_beta)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_3c_overl_int_gw
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                                            nmo
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_dm_virt_local
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: cut_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_from_LLL, my_group_L_starts_im_time, &
                                                            my_group_L_sizes_im_time
      LOGICAL, INTENT(IN)                                :: do_ic_model, do_ic_opt_homo_lumo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_3c_overl_nnP_ic, &
                                                            mat_3c_overl_nnP_ic_reflected
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_3c_overl_int_gw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cut_RI, icol_global, &
                                                            irow_global, n_level_gw, n_level_gw_ref
      LOGICAL                                            :: my_do_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: norm
      TYPE(cp_fm_type), POINTER                          :: fm_mat_mo_coeff_gw
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_3c_overl_int_gw_dummy, &
                                                            mat_3c_overl_int_gw_for_mult
      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_gw, mat_mo_coeff_gw_local, &
                                                            mat_mo_coeff_gw_reflected, mat_norm, &
                                                            mat_work

      CALL timeset(routineN, handle)

      my_do_beta = .FALSE.
      IF (PRESENT(do_beta)) my_do_beta = do_beta

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
      CALL dbcsr_create(matrix=mat_mo_coeff_gw, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mat_mo_coeff_gw, &
                            mat_mo_coeff_gw, &
                            keep_sparsity=.FALSE.)

      ! just remove the blocks which have been set to zero
      CALL dbcsr_filter(mat_mo_coeff_gw, 1.0E-20_dp)

      NULLIFY (mat_mo_coeff_gw_local)
      CALL dbcsr_init_p(mat_mo_coeff_gw_local)
      CALL dbcsr_create(matrix=mat_mo_coeff_gw_local, &
                        template=mat_dm_virt_local%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL replicate_mat_to_subgroup_simple(para_env, para_env_sub, mat_mo_coeff_gw, nmo, &
                                            mat_mo_coeff_gw_local)

      NULLIFY (mat_3c_overl_int_gw_for_mult)
      CALL dbcsr_allocate_matrix_set(mat_3c_overl_int_gw_for_mult, cut_RI)

      DO i_cut_RI = 1, cut_RI

         ALLOCATE (mat_3c_overl_int_gw_for_mult(i_cut_RI)%matrix)
         CALL dbcsr_create(matrix=mat_3c_overl_int_gw_for_mult(i_cut_RI)%matrix, &
                           template=mat_3c_overl_int(i_cut_RI, 1, 1)%matrix)

      END DO

      DO i_cut_RI = 1, cut_RI

         CALL dbcsr_multiply("T", "N", 1.0_dp, mat_mo_coeff_gw_local, &
                             mat_3c_overl_int(i_cut_RI, 1, 1)%matrix, &
                             0.0_dp, mat_3c_overl_int_gw_for_mult(i_cut_RI)%matrix)

      END DO

      CALL fill_mat_3c_overl_int_gw(mat_3c_overl_int_gw(:), mat_3c_overl_int_gw_for_mult, row_from_LLL, &
                                    my_group_L_starts_im_time, my_group_L_sizes_im_time, cut_RI, &
                                    para_env, gw_corr_lev_occ, gw_corr_lev_virt, homo)

      IF (do_ic_model) THEN

         CALL fill_mat_3c_overl_nnP_ic(mat_3c_overl_nnP_ic, mat_3c_overl_int_gw(:), &
                                       mat_mo_coeff_gw, gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                       do_ic_opt_homo_lumo)

         NULLIFY (mat_mo_coeff_gw_reflected)
         CALL dbcsr_init_p(mat_mo_coeff_gw_reflected)
         CALL dbcsr_create(matrix=mat_mo_coeff_gw_reflected, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL reflect_mat_row(mat_mo_coeff_gw_reflected, mat_mo_coeff_gw, para_env, qs_env, unit_nr, do_beta=my_do_beta)

         CALL dbcsr_filter(mat_mo_coeff_gw_reflected, 1.0E-20_dp)

         CALL replicate_mat_to_subgroup_simple(para_env, para_env_sub, mat_mo_coeff_gw_reflected, nmo, &
                                               mat_mo_coeff_gw_local)

         DO i_cut_RI = 1, cut_RI

            CALL dbcsr_multiply("T", "N", 1.0_dp, mat_mo_coeff_gw_local, &
                                mat_3c_overl_int(i_cut_RI, 1, 1)%matrix, &
                                0.0_dp, mat_3c_overl_int_gw_for_mult(i_cut_RI)%matrix)

         END DO

         NULLIFY (mat_3c_overl_int_gw_dummy)
         CALL dbcsr_allocate_matrix_set(mat_3c_overl_int_gw_dummy, gw_corr_lev_occ + gw_corr_lev_virt)

         DO n_level_gw = 1, gw_corr_lev_occ + gw_corr_lev_virt

            ALLOCATE (mat_3c_overl_int_gw_dummy(n_level_gw)%matrix)
            CALL dbcsr_create(matrix=mat_3c_overl_int_gw_dummy(n_level_gw)%matrix, &
                              template=mat_3c_overl_int_gw(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
         END DO

         CALL fill_mat_3c_overl_int_gw(mat_3c_overl_int_gw_dummy, mat_3c_overl_int_gw_for_mult, row_from_LLL, &
                                       my_group_L_starts_im_time, my_group_L_sizes_im_time, cut_RI, &
                                       para_env, gw_corr_lev_occ, gw_corr_lev_virt, homo)

         CALL fill_mat_3c_overl_nnP_ic(mat_3c_overl_nnP_ic_reflected, mat_3c_overl_int_gw_dummy, &
                                       mat_mo_coeff_gw_reflected, gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                       do_ic_opt_homo_lumo)

         ! normalize reflected MOs (they are not properly normalized since high angular momentum basis functions
         ! of the image molecule are not exactly reflected at the image plane (sign problem in p_z function)
         NULLIFY (mat_work)
         CALL dbcsr_init_p(mat_work)
         CALL dbcsr_create(matrix=mat_work, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_norm)
         CALL dbcsr_init_p(mat_norm)
         CALL dbcsr_create(matrix=mat_norm, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_multiply("T", "N", 1.0_dp, mat_mo_coeff_gw_reflected, matrix_s(1)%matrix, 0.0_dp, mat_work)

         CALL dbcsr_multiply("N", "N", 1.0_dp, mat_work, mat_mo_coeff_gw_reflected, 0.0_dp, mat_norm)

         ALLOCATE (norm(nmo))
         norm = 0.0_dp

         CALL dbcsr_get_diag(mat_norm, norm)

         CALL mp_sum(norm, para_env%group)

         DO n_level_gw = 1, gw_corr_lev_occ + gw_corr_lev_virt

            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

            CALL dbcsr_scale(mat_3c_overl_nnP_ic_reflected(n_level_gw)%matrix, 1.0_dp/norm(n_level_gw_ref))

         END DO

         CALL dbcsr_release_p(mat_work)
         CALL dbcsr_release_p(mat_norm)
         CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_gw_dummy)
         CALL dbcsr_release_p(mat_mo_coeff_gw_reflected)

      END IF ! IC Model

      CALL dbcsr_release_p(mat_mo_coeff_gw)
      CALL dbcsr_release_p(mat_mo_coeff_gw_local)
      CALL cp_fm_release(fm_mat_mo_coeff_gw)
      CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_gw_for_mult)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_nnP_ic ...
!> \param mat_3c_overl_int_gw ...
!> \param mat_mo_coeff_gw ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param do_ic_opt_homo_lumo ...
! **************************************************************************************************
   SUBROUTINE fill_mat_3c_overl_nnP_ic(mat_3c_overl_nnP_ic, mat_3c_overl_int_gw, mat_mo_coeff_gw, &
                                       gw_corr_lev_occ, gw_corr_lev_virt, homo, do_ic_opt_homo_lumo)

      TYPE(dbcsr_p_type), DIMENSION(:)                   :: mat_3c_overl_nnP_ic, mat_3c_overl_int_gw
      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_gw
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo
      LOGICAL, INTENT(IN)                                :: do_ic_opt_homo_lumo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_mat_3c_overl_nnP_ic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, col_global, col_offset, col_size, &
                                                            handle, i_col, n_level_gw, &
                                                            n_level_gw_ref, row, row_offset, &
                                                            row_size
      LOGICAL                                            :: is_occ
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_gw_copy

      CALL timeset(routineN, handle)

      NULLIFY (mat_mo_coeff_gw_copy)
      CALL dbcsr_init_p(mat_mo_coeff_gw_copy)
      CALL dbcsr_create(matrix=mat_mo_coeff_gw_copy, &
                        template=mat_mo_coeff_gw, &
                        matrix_type=dbcsr_type_no_symmetry)

      DO n_level_gw = 1, gw_corr_lev_occ + gw_corr_lev_virt

         n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

         CALL dbcsr_copy(mat_mo_coeff_gw_copy, mat_mo_coeff_gw)

         ! set all MO coeff to zero except C_nmu where n=n_level_gw
         CALL dbcsr_iterator_start(iter, mat_mo_coeff_gw_copy)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            ! new image charge model with orbital localization
            IF (do_ic_opt_homo_lumo) THEN

               is_occ = (n_level_gw_ref <= homo)

               DO i_col = 1, col_size

                  col_global = i_col + col_offset - 1

                  IF (is_occ) THEN

                     IF (col_global > homo) THEN

                        data_block(:, i_col) = 0.0_dp

                     END IF

                  ELSE

                     IF (col_global <= homo) THEN

                        data_block(:, i_col) = 0.0_dp

                     END IF

                  END IF

               END DO

               ! Neaton model
            ELSE

               DO i_col = 1, col_size

                  col_global = i_col + col_offset - 1

                  IF (col_global .NE. n_level_gw_ref) THEN

                     data_block(:, i_col) = 0.0_dp

                  END IF

               END DO

            END IF

         END DO

         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_multiply("N", "N", 1.0_dp, mat_3c_overl_int_gw(n_level_gw)%matrix, mat_mo_coeff_gw_copy, &
                             0.0_dp, mat_3c_overl_nnP_ic(n_level_gw)%matrix)

         CALL dbcsr_filter(mat_3c_overl_nnP_ic(n_level_gw)%matrix, 1.0E-30_dp)

      END DO

      CALL dbcsr_release_p(mat_mo_coeff_gw_copy)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int_gw ...
!> \param mat_3c_overl_int_gw_for_mult ...
!> \param row_from_LLL ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param cut_RI ...
!> \param para_env ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
! **************************************************************************************************
   SUBROUTINE fill_mat_3c_overl_int_gw(mat_3c_overl_int_gw, mat_3c_overl_int_gw_for_mult, row_from_LLL, &
                                       my_group_L_starts_im_time, my_group_L_sizes_im_time, cut_RI, &
                                       para_env, gw_corr_lev_occ, gw_corr_lev_virt, homo)

      TYPE(dbcsr_p_type), DIMENSION(:)                   :: mat_3c_overl_int_gw, &
                                                            mat_3c_overl_int_gw_for_mult
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_from_LLL, my_group_L_starts_im_time, &
                                                            my_group_L_sizes_im_time
      INTEGER, INTENT(IN)                                :: cut_RI
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_mat_3c_overl_int_gw', &
         routineP = moduleN//':'//routineN

      INTEGER :: block, col, col_offset, col_offset_data_block, col_prim, col_size, col_size_orig, &
         handle, handle1, i_block, i_col_orig, i_cut_RI, i_row_rec, imepos, imepos_dest, &
         integer_block_counter, LLL, my_group_L_size, my_group_L_start, n_level_gw, offset, &
         offset_buffer, old_block, RI_index_data_block, row, row_offset, row_offset_rec, row_RI, &
         row_size, row_size_rec
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, cols_tmp, cols_to_allocate, &
         entry_counter, num_blocks_rec, num_blocks_send, num_entries_rec, num_entries_send, &
         rows_tmp, rows_to_allocate, sizes_rec, sizes_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL                                            :: first_cycle, is_new_block
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset("GW3c_allocate_stuff", handle1)

      NULLIFY (data_block)

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      num_entries_send(:) = 0

      ALLOCATE (num_blocks_send(0:para_env%num_pe - 1))
      num_blocks_send(:) = 0

      CALL timestop(handle1)

      CALL timeset("GW3c_get_sizes", handle1)

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         CALL dbcsr_iterator_start(iter, mat_3c_overl_int_gw_for_mult(i_cut_RI)%matrix)

         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO LLL = my_group_L_start, my_group_L_start + my_group_L_size - 1

               row_RI = row_from_LLL(LLL)

               CALL dbcsr_get_stored_coordinates(mat_3c_overl_int_gw(1)%matrix, row_RI, col, imepos_dest)

               num_entries_send(imepos_dest) = num_entries_send(imepos_dest) + row_size*col_size/my_group_L_size

               num_blocks_send(imepos_dest) = num_blocks_send(imepos_dest) + col_size/my_group_L_size

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO ! i_cut_RI

      CALL timestop(handle1)

      CALL timeset("GW3c_send_sizes", handle1)

      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe - 1))

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

      CALL timestop(handle1)

      CALL timeset("GW3c_fill_buffer_send", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 10))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send(imepos), 10))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = 0

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         CALL dbcsr_iterator_start(iter, mat_3c_overl_int_gw_for_mult(i_cut_RI)%matrix)

         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO LLL = my_group_L_start, my_group_L_start + my_group_L_size - 1

               row_RI = row_from_LLL(LLL)

               ! the distribution of every mat_3c_overl_int_gw is the same
               CALL dbcsr_get_stored_coordinates(mat_3c_overl_int_gw(1)%matrix, row_RI, col, imepos_dest)

               col_size_orig = col_size/my_group_L_size

               col_offset_data_block = (LLL - my_group_L_start)*col_size_orig + 1

               DO i_col_orig = 1, col_size_orig

                  block = block_counter(imepos_dest) + 1

                  CALL dbcsr_get_stored_coordinates(mat_3c_overl_int_gw(1)%matrix, row_RI, col, imepos_dest)

                  offset = entry_counter(imepos_dest)

                  buffer_send(imepos_dest)%msg(offset + 1:offset + row_size) = &
                     data_block(1:row_size, col_offset_data_block + i_col_orig - 1)

                  buffer_send(imepos_dest)%indx(block, 1) = LLL
                  buffer_send(imepos_dest)%indx(block, 2) = col
                  buffer_send(imepos_dest)%indx(block, 3) = col_offset
                  buffer_send(imepos_dest)%indx(block, 4) = col_size_orig
                  buffer_send(imepos_dest)%indx(block, 5) = row
                  buffer_send(imepos_dest)%indx(block, 6) = row_offset
                  buffer_send(imepos_dest)%indx(block, 7) = row_size
                  buffer_send(imepos_dest)%indx(block, 8) = offset
                  buffer_send(imepos_dest)%indx(block, 9) = i_col_orig

                  entry_counter(imepos_dest) = entry_counter(imepos_dest) + row_size

                  block_counter(imepos_dest) = block_counter(imepos_dest) + 1

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      DEALLOCATE (entry_counter, block_counter)

      CALL timestop(handle1)

      CALL timeset("GW3c_comm_data", handle1)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      CALL timestop(handle1)

      CALL timeset("GW3c_alloc_blocks", handle1)

      ALLOCATE (rows_to_allocate(1))
      rows_to_allocate = 0

      ALLOCATE (cols_to_allocate(1))
      cols_to_allocate = 0

      integer_block_counter = 1

      first_cycle = .TRUE.

      DO imepos = 0, para_env%num_pe - 1

         DO i_block = 1, num_blocks_rec(imepos)

            row_RI = row_from_LLL(buffer_rec(imepos)%indx(i_block, 1))
            col_prim = buffer_rec(imepos)%indx(i_block, 2)

            is_new_block = .TRUE.

            ! check whether block is already there
            DO old_block = 1, integer_block_counter

               IF (row_RI == rows_to_allocate(old_block) .AND. col_prim == cols_to_allocate(old_block)) THEN

                  is_new_block = .FALSE.

               END IF

            END DO

            IF (is_new_block) THEN

               IF (first_cycle) THEN

                  rows_to_allocate(1) = row_RI
                  cols_to_allocate(1) = col_prim

                  first_cycle = .FALSE.

               ELSE

                  ALLOCATE (rows_tmp(integer_block_counter))
                  ALLOCATE (cols_tmp(integer_block_counter))

                  rows_tmp(1:integer_block_counter) = rows_to_allocate(1:integer_block_counter)
                  cols_tmp(1:integer_block_counter) = cols_to_allocate(1:integer_block_counter)

                  DEALLOCATE (rows_to_allocate)
                  DEALLOCATE (cols_to_allocate)

                  ALLOCATE (rows_to_allocate(integer_block_counter + 1))
                  ALLOCATE (cols_to_allocate(integer_block_counter + 1))

                  rows_to_allocate(1:integer_block_counter) = rows_tmp(1:integer_block_counter)
                  cols_to_allocate(1:integer_block_counter) = cols_tmp(1:integer_block_counter)

                  DEALLOCATE (rows_tmp)
                  DEALLOCATE (cols_tmp)

                  rows_to_allocate(integer_block_counter + 1) = row_RI
                  cols_to_allocate(integer_block_counter + 1) = col_prim

                  integer_block_counter = integer_block_counter + 1

               END IF

            END IF

         END DO

      END DO

      CALL timestop(handle1)

      CALL timeset("GW3c_fill_mat", handle1)

      DO n_level_gw = 1, gw_corr_lev_occ + gw_corr_lev_virt

         IF (rows_to_allocate(1) .NE. 0) THEN

            CALL dbcsr_reserve_blocks(mat_3c_overl_int_gw(n_level_gw)%matrix, &
                                      rows=rows_to_allocate, &
                                      cols=cols_to_allocate)

         END IF

         CALL dbcsr_finalize(mat_3c_overl_int_gw(n_level_gw)%matrix)

         DO imepos = 0, para_env%num_pe - 1

            DO i_block = 1, num_blocks_rec(imepos)

               row_RI = row_from_LLL(buffer_rec(imepos)%indx(i_block, 1))
               col_prim = buffer_rec(imepos)%indx(i_block, 2)

               CALL dbcsr_iterator_start(iter, mat_3c_overl_int_gw(n_level_gw)%matrix)
               DO WHILE (dbcsr_iterator_blocks_left(iter))

                  CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                                 row_size=row_size, col_size=col_size, &
                                                 row_offset=row_offset, col_offset=col_offset)

                  IF (row_RI == row .AND. col_prim == col) THEN

                     ! now check whether we have the n_level_gw in the buffer_rec
                     LLL = buffer_rec(imepos)%indx(i_block, 1)
                     row_offset_rec = buffer_rec(imepos)%indx(i_block, 6)
                     row_size_rec = buffer_rec(imepos)%indx(i_block, 7)
                     offset_buffer = buffer_rec(imepos)%indx(i_block, 8)
                     i_col_orig = buffer_rec(imepos)%indx(i_block, 9)

                     DO i_row_rec = row_offset_rec, row_offset_rec + row_size_rec - 1

                        IF (i_row_rec == n_level_gw + homo - gw_corr_lev_occ) THEN

                           RI_index_data_block = LLL - row_offset + 1

                           data_block(RI_index_data_block, i_col_orig) = &
                              buffer_rec(imepos)%msg(offset_buffer + 1 + i_row_rec - row_offset_rec)

                        END IF

                     END DO

                  END IF

               END DO

               CALL dbcsr_iterator_stop(iter)

            END DO

         END DO

      END DO

      DEALLOCATE (rows_to_allocate)
      DEALLOCATE (cols_to_allocate)

      CALL timestop(handle1)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
      END DO

      DEALLOCATE (buffer_send, buffer_rec)

      DEALLOCATE (num_entries_send, num_blocks_send, req_array, num_entries_rec, num_blocks_rec)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param mat_dm_global ...
!> \param nmo ...
!> \param mat_local ...
! **************************************************************************************************
   SUBROUTINE replicate_mat_to_subgroup_simple(para_env, para_env_sub, mat_dm_global, nmo, &
                                               mat_local)

      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(dbcsr_type), POINTER                          :: mat_dm_global
      INTEGER, INTENT(IN)                                :: nmo
      TYPE(dbcsr_type), POINTER                          :: mat_local

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_mat_to_subgroup_simple', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_counter, block_offset, block_size, col, col_from_buffer, col_offset, &
         col_size, handle, handle1, i_block, i_entry, i_mepos, igroup, imepos, imepos_sub, &
         msg_offset, nblkrows_total, ngroup, num_blocks, num_pe_sub, offset, row, row_from_buffer, &
         row_offset, row_size, total_num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_counter, cols_to_allocate_all, entry_counter, &
         num_entries_blocks_rec, num_entries_blocks_send, row_block_from_index, &
         rows_to_allocate_all, sizes_rec, sizes_send
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_offset, row_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset("get_sizes_D", handle1)

      NULLIFY (data_block)

      ALLOCATE (num_entries_blocks_send(0:2*para_env%num_pe - 1))
      num_entries_blocks_send(:) = 0

      ALLOCATE (num_entries_blocks_rec(0:2*para_env%num_pe - 1))
      num_entries_blocks_rec(:) = 0

      ngroup = para_env%num_pe/para_env_sub%num_pe

      num_pe_sub = para_env_sub%num_pe

      CALL dbcsr_iterator_start(iter, mat_dm_global)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imepos_sub)

         DO igroup = 0, ngroup - 1

            imepos = imepos_sub + igroup*num_pe_sub

            num_entries_blocks_send(2*imepos) = num_entries_blocks_send(2*imepos) + row_size*col_size
            num_entries_blocks_send(2*imepos + 1) = num_entries_blocks_send(2*imepos + 1) + 1

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("send_sizes_D_1", handle1)

      total_num_entries = SUM(num_entries_blocks_send)
      CALL mp_sum(total_num_entries, para_env%group)

      CALL timestop(handle1)

      CALL timeset("send_sizes_D_2", handle1)

      IF (para_env%num_pe > 1) THEN
         CALL mp_alltoall(num_entries_blocks_send, num_entries_blocks_rec, 2, para_env%group)
      ELSE
         num_entries_blocks_rec(0:1) = num_entries_blocks_send(0:1)
      END IF

      CALL timestop(handle1)

      CALL timeset("get_data_D", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_blocks_rec(2*imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_blocks_send(2*imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_blocks_rec(2*imepos + 1), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_blocks_send(2*imepos + 1), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ALLOCATE (blk_counter(0:para_env%num_pe - 1))
      blk_counter = 0

      CALL dbcsr_iterator_start(iter, mat_dm_global)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imepos_sub)

         DO igroup = 0, ngroup - 1

            imepos = imepos_sub + igroup*num_pe_sub

            msg_offset = entry_counter(imepos)

            block_size = row_size*col_size

            buffer_send(imepos)%msg(msg_offset + 1:msg_offset + block_size) = &
               RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

            entry_counter(imepos) = entry_counter(imepos) + block_size

            blk_counter(imepos) = blk_counter(imepos) + 1

            block_offset = blk_counter(imepos)

            buffer_send(imepos)%indx(block_offset, 1) = row
            buffer_send(imepos)%indx(block_offset, 2) = col
            buffer_send(imepos)%indx(block_offset, 3) = msg_offset

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("send_data_D", handle1)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (sizes_rec(0:para_env%num_pe - 1))
      ALLOCATE (sizes_send(0:para_env%num_pe - 1))

      DO imepos = 0, para_env%num_pe - 1

         sizes_send(imepos) = num_entries_blocks_send(2*imepos)
         sizes_rec(imepos) = num_entries_blocks_rec(2*imepos)

      END DO

      CALL communicate_buffer(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array, sizes_rec, sizes_send)

      CALL timestop(handle1)

      CALL timeset("row_block_from_index", handle1)

      CALL dbcsr_get_info(mat_local, &
                          nblkrows_total=nblkrows_total, &
                          row_blk_offset=row_blk_offset, &
                          row_blk_size=row_blk_size)

      ALLOCATE (row_block_from_index(nmo))
      row_block_from_index = 0

      DO i_entry = 1, nmo
         DO i_block = 1, nblkrows_total

            IF (i_entry >= row_blk_offset(i_block) .AND. i_entry <= row_blk_offset(i_block) + row_blk_size(i_block) - 1) THEN

               row_block_from_index(i_entry) = i_block

            END IF

         END DO
      END DO

      CALL timestop(handle1)

      CALL timeset("reserve_blocks_D", handle1)

      num_blocks = 0

      ! get the number of blocks, which have to be allocated
      DO imepos = 0, para_env%num_pe - 1
         num_blocks = num_blocks + num_entries_blocks_rec(2*imepos + 1)
      END DO

      ALLOCATE (rows_to_allocate_all(num_blocks))
      rows_to_allocate_all = 0

      ALLOCATE (cols_to_allocate_all(num_blocks))
      cols_to_allocate_all = 0

      block_counter = 0

      DO i_mepos = 0, para_env%num_pe - 1

         DO i_block = 1, num_entries_blocks_rec(2*i_mepos + 1)

            block_counter = block_counter + 1

            rows_to_allocate_all(block_counter) = buffer_rec(i_mepos)%indx(i_block, 1)
            cols_to_allocate_all(block_counter) = buffer_rec(i_mepos)%indx(i_block, 2)

         END DO

      END DO

      CALL dbcsr_reserve_blocks(mat_local, &
                                rows=rows_to_allocate_all, &
                                cols=cols_to_allocate_all)
      CALL dbcsr_finalize(mat_local)
      CALL dbcsr_set(mat_local, 0.0_dp)

      CALL timestop(handle1)

      CALL timeset("fill_mat_local", handle1)

      ! fill the dbcsr matrix
      CALL dbcsr_iterator_start(iter, mat_local)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         DO imepos = 0, para_env%num_pe - 1

            DO i_block = 1, num_entries_blocks_rec(2*imepos + 1)

               row_from_buffer = buffer_rec(imepos)%indx(i_block, 1)
               col_from_buffer = buffer_rec(imepos)%indx(i_block, 2)
               offset = buffer_rec(imepos)%indx(i_block, 3)

               IF (row == row_from_buffer .AND. col == col_from_buffer) THEN

                  data_block(1:row_size, 1:col_size) = &
                     RESHAPE(buffer_rec(imepos)%msg(offset + 1:offset + row_size*col_size), &
                             (/row_size, col_size/))

               END IF

            END DO

         END DO

      END DO ! blocks

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (entry_counter)
      DEALLOCATE (blk_counter)
      DEALLOCATE (row_block_from_index)
      DEALLOCATE (num_entries_blocks_send, num_entries_blocks_rec)

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

      ! matrix has to be quadratic of course
      CPASSERT(nblkrows_total == nblkcols_total)

      natom = SIZE(particle_set)
      CPASSERT(natom == nblkcols_total)

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

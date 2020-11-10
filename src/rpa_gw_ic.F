!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate image charge corrections
!> \par History
!>      06.2019 Moved from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE rpa_gw_ic
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_get_info, dbcsr_get_stored_coordinates, dbcsr_init_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_release_p, &
        dbcsr_reserve_all_blocks, dbcsr_scalar, dbcsr_set, dbcsr_transposed, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, dbcsr_t_create, &
        dbcsr_t_destroy, dbcsr_t_get_block, dbcsr_t_get_info, dbcsr_t_iterator_blocks_left, &
        dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, dbcsr_t_iterator_stop, &
        dbcsr_t_iterator_type, dbcsr_t_nblks_total, dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, &
        dbcsr_t_pgrid_type, dbcsr_t_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_dims_create,&
                                              mp_sum
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE physcon,                         ONLY: evolt
   USE qs_tensors_types,                ONLY: create_2c_tensor
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_ic'

   PUBLIC :: calculate_ic_correction

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param mat_SinvVSinv ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param gw_corr_lev_tot ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param unit_nr ...
!> \param print_ic_values ...
!> \param para_env ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_ic_correction(Eigenval, mat_SinvVSinv, &
                                      t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, gw_corr_lev_tot, &
                                      gw_corr_lev_occ, gw_corr_lev_virt, homo, unit_nr, &
                                      print_ic_values, para_env, &
                                      do_alpha, do_beta)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_SinvVSinv
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      INTEGER, INTENT(IN)                                :: gw_corr_lev_tot, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, homo, unit_nr
      LOGICAL, INTENT(IN)                                :: print_ic_values
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_ic_correction'

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, mo_end, mo_start, n_level_gw, &
                                                            n_level_gw_ref
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2, sizes_RI_split
      INTEGER, DIMENSION(2)                              :: pdims
      LOGICAL                                            :: do_closed_shell, my_do_alpha, my_do_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Delta_Sigma_Neaton
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid_2d
      TYPE(dbcsr_t_type) :: t_3c_overl_nnP_ic_reflected_ctr, t_SinvVSinv, t_SinvVSinv_tmp

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

      IF (PRESENT(do_beta)) THEN
         my_do_beta = do_beta
      ELSE
         my_do_beta = .FALSE.
      END IF

      do_closed_shell = .NOT. (my_do_alpha .OR. my_do_beta)

      ALLOCATE (Delta_Sigma_Neaton(gw_corr_lev_tot))
      Delta_Sigma_Neaton = 0.0_dp

      mo_start = homo - gw_corr_lev_occ + 1
      mo_end = homo + gw_corr_lev_virt
      CPASSERT(mo_end - mo_start + 1 == gw_corr_lev_tot)

      ALLOCATE (sizes_RI_split(dbcsr_t_nblks_total(t_3c_overl_nnP_ic_reflected, 1)))
      CALL dbcsr_t_get_info(t_3c_overl_nnP_ic_reflected, blk_size_1=sizes_RI_split)

      CALL dbcsr_t_create(mat_SinvVSinv%matrix, t_SinvVSinv_tmp)
      CALL dbcsr_t_copy_matrix_to_tensor(mat_SinvVSinv%matrix, t_SinvVSinv_tmp)
      pdims = 0
      CALL mp_dims_create(para_env%num_pe, pdims)
      CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid_2d)
      CALL create_2c_tensor(t_SinvVSinv, dist_1, dist_2, pgrid_2d, sizes_RI_split, sizes_RI_split, &
                            name="(RI|RI)")
      DEALLOCATE (dist_1, dist_2)
      CALL dbcsr_t_pgrid_destroy(pgrid_2d)

      CALL dbcsr_t_copy(t_SinvVSinv_tmp, t_SinvVSinv)
      CALL dbcsr_t_destroy(t_SinvVSinv_tmp)
      CALL dbcsr_t_create(t_3c_overl_nnP_ic_reflected, t_3c_overl_nnP_ic_reflected_ctr)
      CALL dbcsr_t_contract(dbcsr_scalar(0.5_dp), t_SinvVSinv, t_3c_overl_nnP_ic_reflected, &
                            dbcsr_scalar(0.0_dp), t_3c_overl_nnP_ic_reflected_ctr, &
                            contract_1=[2], notcontract_1=[1], &
                            contract_2=[1], notcontract_2=[2, 3], &
                            map_1=[1], map_2=[2, 3])

      CALL trace_ic_gw(t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected_ctr, Delta_Sigma_Neaton, [mo_start, mo_end], para_env)

      Delta_Sigma_Neaton(gw_corr_lev_occ + 1:) = -Delta_Sigma_Neaton(gw_corr_lev_occ + 1:)

      CALL dbcsr_t_destroy(t_SinvVSinv)
      CALL dbcsr_t_destroy(t_3c_overl_nnP_ic_reflected_ctr)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') 'Single-electron energies with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '----------------------------------------------------------'
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') 'Single-electron energies of alpha spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '-------------------------------------------------------------------------'
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') 'Single-electron energies of beta spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '------------------------------------------------------------------------'
         END IF

         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T3,A)') 'Reference for the ic: Neaton et al., PRL 97, 216405 (2006)'
         WRITE (unit_nr, *) ' '

         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T14,2A)') 'MO     E_n before ic corr           Delta E_ic', &
            '    E_n after ic corr'

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
               n_level_gw_ref, ' ( ', occ_virt, ')  ', &
               Eigenval(n_level_gw_ref)*evolt, &
               Delta_Sigma_Neaton(n_level_gw)*evolt, &
               (Eigenval(n_level_gw_ref) + Delta_Sigma_Neaton(n_level_gw))*evolt

         END DO

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F57.2)') 'IC HOMO-LUMO gap (eV)', (Eigenval(homo + 1) + &
                                                                      Delta_Sigma_Neaton(gw_corr_lev_occ + 1) - &
                                                                      Eigenval(homo) - &
                                                                      Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F51.2)') 'Alpha IC HOMO-LUMO gap (eV)', (Eigenval(homo + 1) + &
                                                                            Delta_Sigma_Neaton(gw_corr_lev_occ + 1) - &
                                                                            Eigenval(homo) - &
                                                                            Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F52.2)') 'Beta IC HOMO-LUMO gap (eV)', (Eigenval(homo + 1) + &
                                                                           Delta_Sigma_Neaton(gw_corr_lev_occ + 1) - &
                                                                           Eigenval(homo) - &
                                                                           Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
         END IF

         IF (print_ic_values) THEN

            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'Horizontal list for copying the image charge corrections for use as input:'
            WRITE (unit_nr, '(*(F7.3))') (Delta_Sigma_Neaton(n_level_gw)*evolt, &
                                          n_level_gw=1, gw_corr_lev_tot)

         END IF

      END IF

      Eigenval(homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt) = Eigenval(homo - gw_corr_lev_occ + 1: &
                                                                              homo + gw_corr_lev_virt) &
                                                                     + Delta_Sigma_Neaton(1:gw_corr_lev_tot)

      CALL timestop(handle)

   END SUBROUTINE calculate_ic_correction

! **************************************************************************************************
!> \brief ...
!> \param t3c_1 ...
!> \param t3c_2 ...
!> \param Delta_Sigma_Neaton ...
!> \param mo_bounds ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE trace_ic_gw(t3c_1, t3c_2, Delta_Sigma_Neaton, mo_bounds, para_env)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t3c_1, t3c_2
      REAL(dp), DIMENSION(:), INTENT(INOUT)              :: Delta_Sigma_Neaton
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mo_bounds
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trace_ic_gw'

      INTEGER                                            :: blk, handle, n, n_end, n_end_block, &
                                                            n_start, n_start_block
      INTEGER, DIMENSION(3)                              :: boff, bsize, ind
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_1, block_2
      TYPE(dbcsr_t_iterator_type)                        :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_t_iterator_start(iter, t3c_1)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_size=bsize, blk_offset=boff)
         IF (ind(2) /= ind(3)) CYCLE
         CALL dbcsr_t_get_block(t3c_1, ind, block_1, found)
         CPASSERT(found)
         CALL dbcsr_t_get_block(t3c_2, ind, block_2, found)
         IF (.NOT. found) CYCLE

         IF (boff(3) < mo_bounds(1)) THEN
            n_start_block = mo_bounds(1) - boff(3) + 1
            n_start = 1
         ELSE
            n_start_block = 1
            n_start = boff(3) - mo_bounds(1) + 1
         END IF

         IF (boff(3) + bsize(3) - 1 > mo_bounds(2)) THEN
            n_end_block = mo_bounds(2) - boff(3) + 1
            n_end = mo_bounds(2) - mo_bounds(1) + 1
         ELSE
            n_end_block = bsize(3)
            n_end = boff(3) + bsize(3) - mo_bounds(1)
         END IF

         Delta_Sigma_Neaton(n_start:n_end) = &
            Delta_Sigma_Neaton(n_start:n_end) + &
            (/(DOT_PRODUCT(block_1(:, n, n), &
                           block_2(:, n, n)), &
               n=n_start_block, n_end_block)/)
         DEALLOCATE (block_1, block_2)
      END DO
      CALL dbcsr_t_iterator_stop(iter)

      CALL mp_sum(Delta_Sigma_Neaton, para_env%group)

      CALL timestop(handle)

   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param fm_mat_M_occ ...
!> \param unit_nr ...
!> \param gw_corr_lev_occ ...
!> \param para_env ...
!> \param do_homo ...
!> \param do_lumo ...
! **************************************************************************************************
   SUBROUTINE print_Neaton_value(fm_mat_M_occ, unit_nr, gw_corr_lev_occ, para_env, do_homo, do_lumo)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_occ
      INTEGER, INTENT(IN)                                :: unit_nr, gw_corr_lev_occ
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_homo, do_lumo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_Neaton_value'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_homo, my_do_lumo
      REAL(KIND=dp)                                      :: Neaton_value

      CALL timeset(routineN, handle)

      my_do_homo = .FALSE.
      IF (PRESENT(do_homo)) my_do_homo = do_homo

      my_do_lumo = .FALSE.
      IF (PRESENT(do_lumo)) my_do_lumo = do_lumo

      CPASSERT(my_do_homo .NEQV. my_do_lumo)

      CALL cp_fm_get_info(matrix=fm_mat_M_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      Neaton_value = 0.0_dp

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == gw_corr_lev_occ .AND. i_global == gw_corr_lev_occ) THEN
               Neaton_value = fm_mat_M_occ%local_data(iiB, jjB)
            END IF
         END DO
      END DO

      CALL mp_sum(Neaton_value, para_env%group)

      IF (unit_nr > 0 .AND. my_do_homo) WRITE (unit_nr, '(T3,A,F47.2)') 'Neaton value of the HOMO (eV): ', &
         Neaton_value*evolt
      IF (unit_nr > 0 .AND. my_do_lumo) WRITE (unit_nr, '(T3,A,F47.2)') 'Neaton value of the LUMO (eV): ', &
         Neaton_value*evolt

      CALL timestop(handle)

   END SUBROUTINE print_Neaton_value

! **************************************************************************************************
!> \brief ...
!> \param coeff_homo ...
!> \param fm_mat_U_occ ...
!> \param para_env ...
!> \param homo ...
!> \param gw_corr_lev_occ ...
!> \param do_lumo ...
! **************************************************************************************************
   SUBROUTINE update_coeff_homo(coeff_homo, fm_mat_U_occ, para_env, homo, gw_corr_lev_occ, do_lumo)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: coeff_homo
      TYPE(cp_fm_type), POINTER                          :: fm_mat_U_occ
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: homo, gw_corr_lev_occ
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_lumo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_coeff_homo'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_homo
      REAL(KIND=dp)                                      :: norm_coeff_homo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff_homo_update, &
                                                            coeff_homo_update_orthog

      CALL timeset(routineN, handle)

      my_do_homo = .TRUE.
      IF (PRESENT(do_lumo)) my_do_homo = .NOT. do_lumo

      CALL cp_fm_get_info(matrix=fm_mat_U_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ALLOCATE (coeff_homo_update(homo))
      coeff_homo_update = 0.0_dp

      ALLOCATE (coeff_homo_update_orthog(homo))
      coeff_homo_update_orthog = 0.0_dp

      ! take the eigenvector belongning to the largest eigenvalue
      DO iiB = 1, nrow_local

         i_global = row_indices(iiB)

         DO jjB = 1, ncol_local

            j_global = col_indices(jjB)

            IF (j_global .NE. 1) CYCLE

            IF (my_do_homo) THEN

               coeff_homo_update(i_global + homo - gw_corr_lev_occ) = fm_mat_U_occ%local_data(iiB, jjB)

            ELSE

               coeff_homo_update(i_global) = fm_mat_U_occ%local_data(iiB, jjB)

            END IF

         END DO

      END DO

      CALL mp_sum(coeff_homo_update, para_env%group)

      norm_coeff_homo = NORM2(coeff_homo_update)

      coeff_homo(:) = coeff_homo_update(:)/norm_coeff_homo

      CALL timestop(handle)

   END SUBROUTINE update_coeff_homo

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_M_virt ...
!> \param mat_N_virt_dbcsr ...
!> \param matrix_s ...
!> \param Eigenval ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fill_fm_mat_M_virt(fm_mat_M_virt, mat_N_virt_dbcsr, matrix_s, Eigenval, gw_corr_lev_virt, homo, para_env)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_virt
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_N_virt_dbcsr, matrix_s
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: gw_corr_lev_virt, homo
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_mat_M_virt'

      INTEGER :: col, col_offset, col_size, handle, i_col, i_global, i_index, i_row, iiB, imepos, &
         j_global, jjB, m_level_gw, n_level_gw, nblkrows_total, ncol_local, nfullrows_total, &
         nrow_local, offset, row, row_block, row_index, row_offset, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_from_indx, entry_counter, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_blk_offset, &
                                                            row_blk_sizes, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_tmp, matrix_tmp_2
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      NULLIFY (matrix_tmp)
      CALL dbcsr_init_p(matrix_tmp)
      CALL cp_dbcsr_m_by_n_from_template(matrix_tmp, template=matrix_s(1)%matrix, m=gw_corr_lev_virt, n=gw_corr_lev_virt, &
                                         sym=dbcsr_type_no_symmetry)

      NULLIFY (matrix_tmp_2)
      CALL dbcsr_init_p(matrix_tmp_2)
      CALL dbcsr_create(matrix=matrix_tmp_2, &
                        template=matrix_tmp, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_reserve_all_blocks(matrix_tmp)
      CALL dbcsr_set(matrix_tmp, 0.0_dp)

      CALL dbcsr_get_info(matrix_tmp, &
                          nblkrows_total=nblkrows_total, &
                          nfullrows_total=nfullrows_total, &
                          row_blk_offset=row_blk_offset, &
                          row_blk_size=row_blk_sizes)

      ALLOCATE (blk_from_indx(nfullrows_total))

      DO row_index = 1, nfullrows_total

         DO row_block = 1, nblkrows_total

            IF (row_index >= row_blk_offset(row_block) .AND. &
                row_index <= row_blk_offset(row_block) + row_blk_sizes(row_block) - 1) THEN

               blk_from_indx(row_index) = row_block

            END IF

         END DO

      END DO

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      num_entries_send = 0
      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      num_entries_rec = 0

      DO n_level_gw = 1, gw_corr_lev_virt

         CALL dbcsr_iterator_start(iter, mat_N_virt_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            DO i_row = 1, row_size

               m_level_gw = row_offset - 1 + i_row - homo

               IF (m_level_gw < 1) CYCLE

               IF (m_level_gw > gw_corr_lev_virt) CYCLE

               CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                 blk_from_indx(n_level_gw), imepos)

               num_entries_send(imepos) = num_entries_send(imepos) + 1

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 1

      DO n_level_gw = 1, gw_corr_lev_virt

         CALL dbcsr_iterator_start(iter, mat_N_virt_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            DO i_row = 1, row_size

               m_level_gw = row_offset - 1 + i_row - homo

               IF (m_level_gw < 1) CYCLE

               IF (m_level_gw > gw_corr_lev_virt) CYCLE

               CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                 blk_from_indx(n_level_gw), imepos)

               offset = entry_counter(imepos)

               buffer_send(imepos)%msg(offset) = data_block(i_row, 1)
               buffer_send(imepos)%indx(offset, 1) = m_level_gw
               buffer_send(imepos)%indx(offset, 2) = n_level_gw

               entry_counter(imepos) = entry_counter(imepos) + 1

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array)

      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size, &
                                        col_offset=col_offset, col_size=col_size)

         DO imepos = 0, para_env%num_pe - 1

            DO i_index = 1, num_entries_rec(imepos)

               IF (buffer_rec(imepos)%indx(i_index, 1) >= row_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 1) <= row_offset + row_size - 1 .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) >= col_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) <= col_offset + col_size - 1) THEN

                  i_row = buffer_rec(imepos)%indx(i_index, 1) - row_offset + 1
                  i_col = buffer_rec(imepos)%indx(i_index, 2) - col_offset + 1

                  data_block(i_row, i_col) = buffer_rec(imepos)%msg(i_index)

               END IF

            END DO

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      ! symmetrize the result
      CALL dbcsr_transposed(matrix_tmp_2, matrix_tmp)
      CALL dbcsr_add(matrix_tmp, matrix_tmp_2, 0.5_dp, 0.5_dp)

      CALL copy_dbcsr_to_fm(matrix_tmp, fm_mat_M_virt)

      ! add the eigenvalue on the diag of M
      CALL cp_fm_get_info(matrix=fm_mat_M_virt, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global) THEN
               fm_mat_M_virt%local_data(iiB, jjB) = fm_mat_M_virt%local_data(iiB, jjB) + &
                                                    Eigenval(i_global + homo)
            END IF
         END DO
      END DO

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      CALL dbcsr_release_p(matrix_tmp)
      CALL dbcsr_release_p(matrix_tmp_2)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_mat_M_virt

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_M_occ ...
!> \param mat_N_occ_dbcsr ...
!> \param matrix_s ...
!> \param Eigenval ...
!> \param gw_corr_lev_occ ...
!> \param homo ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fill_fm_mat_M_occ(fm_mat_M_occ, mat_N_occ_dbcsr, matrix_s, Eigenval, gw_corr_lev_occ, homo, para_env)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_N_occ_dbcsr, matrix_s
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, homo
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fill_fm_mat_M_occ'

      INTEGER :: col, col_offset, col_size, handle, i_col, i_global, i_index, i_row, iiB, imepos, &
         j_global, jjB, m_level_gw, n_level_gw, nblkrows_total, ncol_local, nfullrows_total, &
         nrow_local, offset, row, row_block, row_index, row_offset, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_from_indx, entry_counter, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_blk_offset, &
                                                            row_blk_sizes, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_tmp
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      NULLIFY (matrix_tmp)
      CALL dbcsr_init_p(matrix_tmp)
      CALL cp_dbcsr_m_by_n_from_template(matrix_tmp, template=matrix_s(1)%matrix, m=gw_corr_lev_occ, n=gw_corr_lev_occ, &
                                         sym=dbcsr_type_no_symmetry)

      CALL dbcsr_reserve_all_blocks(matrix_tmp)
      CALL dbcsr_set(matrix_tmp, 0.0_dp)

      CALL dbcsr_get_info(matrix_tmp, &
                          nblkrows_total=nblkrows_total, &
                          nfullrows_total=nfullrows_total, &
                          row_blk_offset=row_blk_offset, &
                          row_blk_size=row_blk_sizes)

      ALLOCATE (blk_from_indx(nfullrows_total))

      DO row_index = 1, nfullrows_total

         DO row_block = 1, nblkrows_total

            IF (row_index >= row_blk_offset(row_block) .AND. &
                row_index <= row_blk_offset(row_block) + row_blk_sizes(row_block) - 1) THEN

               blk_from_indx(row_index) = row_block

            END IF

         END DO

      END DO

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      num_entries_send = 0
      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      num_entries_rec = 0

      DO n_level_gw = 1, gw_corr_lev_occ

         CALL dbcsr_iterator_start(iter, mat_N_occ_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            IF (row_offset + row_size - 1 <= homo) THEN

               DO i_row = 1, row_size

                  m_level_gw = row_offset - 1 + i_row - (homo - gw_corr_lev_occ)

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  num_entries_send(imepos) = num_entries_send(imepos) + 1

               END DO

            ELSE IF (row_offset <= homo) THEN

               DO m_level_gw = row_offset - (homo - gw_corr_lev_occ), gw_corr_lev_occ

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  num_entries_send(imepos) = num_entries_send(imepos) + 1

               END DO

            END IF

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 1

      DO n_level_gw = 1, gw_corr_lev_occ

         CALL dbcsr_iterator_start(iter, mat_N_occ_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            IF (row_offset + row_size - 1 <= homo) THEN

               DO i_row = 1, row_size

                  m_level_gw = row_offset - 1 + i_row - (homo - gw_corr_lev_occ)

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  offset = entry_counter(imepos)

                  buffer_send(imepos)%msg(offset) = data_block(i_row, 1)
                  buffer_send(imepos)%indx(offset, 1) = m_level_gw
                  buffer_send(imepos)%indx(offset, 2) = n_level_gw

                  entry_counter(imepos) = entry_counter(imepos) + 1

               END DO

            ELSE IF (row_offset <= homo) THEN

               DO m_level_gw = row_offset - (homo - gw_corr_lev_occ), gw_corr_lev_occ

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  offset = entry_counter(imepos)

                  i_row = m_level_gw + (homo - gw_corr_lev_occ) - row_offset + 1

                  buffer_send(imepos)%msg(offset) = data_block(i_row, 1)
                  buffer_send(imepos)%indx(offset, 1) = m_level_gw
                  buffer_send(imepos)%indx(offset, 2) = n_level_gw

                  entry_counter(imepos) = entry_counter(imepos) + 1

               END DO

            END IF

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array)

      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size, &
                                        col_offset=col_offset, col_size=col_size)

         DO imepos = 0, para_env%num_pe - 1

            DO i_index = 1, num_entries_rec(imepos)

               IF (buffer_rec(imepos)%indx(i_index, 1) >= row_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 1) <= row_offset + row_size - 1 .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) >= col_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) <= col_offset + col_size - 1) THEN

                  i_row = buffer_rec(imepos)%indx(i_index, 1) - row_offset + 1
                  i_col = buffer_rec(imepos)%indx(i_index, 2) - col_offset + 1

                  data_block(i_row, i_col) = buffer_rec(imepos)%msg(i_index)

               END IF

            END DO

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL copy_dbcsr_to_fm(matrix_tmp, fm_mat_M_occ)

      ! add the eigenvalue on the diag of M
      CALL cp_fm_get_info(matrix=fm_mat_M_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global) THEN
               fm_mat_M_occ%local_data(iiB, jjB) = fm_mat_M_occ%local_data(iiB, jjB) + &
                                                   Eigenval(i_global + homo - gw_corr_lev_occ)
            END IF
         END DO
      END DO

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      CALL dbcsr_release_p(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_mat_M_occ

! **************************************************************************************************
!> \brief ...
!> \param coeff_dbcsr ...
!> \param coeff ...
!> \param homo ...
! **************************************************************************************************
   SUBROUTINE fill_coeff_dbcsr_occ(coeff_dbcsr, coeff, homo)
      TYPE(dbcsr_type), POINTER                          :: coeff_dbcsr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coeff
      INTEGER, INTENT(IN)                                :: homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_coeff_dbcsr_occ'

      INTEGER                                            :: col, end_data_block, handle, row, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, coeff_dbcsr)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size)

         IF (row_offset + row_size - 1 <= homo) THEN

            data_block(1:row_size, 1) = coeff(row_offset:row_offset + row_size - 1)

         ELSE IF (row_offset <= homo) THEN

            end_data_block = homo - row_offset + 1

            data_block(1:end_data_block, 1) = coeff(row_offset:row_offset + end_data_block - 1)

         END IF

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE fill_coeff_dbcsr_occ

! **************************************************************************************************
!> \brief ...
!> \param coeff_dbcsr ...
!> \param coeff ...
!> \param homo ...
! **************************************************************************************************
   SUBROUTINE fill_coeff_dbcsr_virt(coeff_dbcsr, coeff, homo)
      TYPE(dbcsr_type), POINTER                          :: coeff_dbcsr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coeff
      INTEGER, INTENT(IN)                                :: homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_coeff_dbcsr_virt'

      INTEGER                                            :: col, handle, row, row_offset, row_size, &
                                                            start_data_block
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, coeff_dbcsr)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size)

         IF (row_offset > homo) THEN

            data_block(1:row_size, 1) = coeff(row_offset - homo:row_offset + row_size - homo - 1)

         ELSE IF (row_offset + row_size - 1 > homo) THEN

            start_data_block = homo - row_offset + 1

            data_block(start_data_block:row_size, 1) = coeff(1:row_offset + row_size - homo)

         END IF

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE fill_coeff_dbcsr_virt

END MODULE rpa_gw_ic

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 08.2023
! **************************************************************************************************
MODULE gw_utils_communication
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_filter,&
                                              dbcsr_get_stored_coordinates,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_dbcsr_redistribution,         ONLY: cp_dbcsr_to_subgroups
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbt_api,                         ONLY: dbt_clear,&
                                              dbt_copy,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE kinds,                           ONLY: dp
   USE matrix_communication,            ONLY: buffer_type => matrix_buffer_type,&
                                              exchange_matrix_buffers
   USE message_passing,                 ONLY: mp_para_env_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils_communication'

   PUBLIC :: local_dbt_to_global_mat, fm_to_local_tensor, local_matrix_to_global_matrix, local_dbt_to_global_fm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_global ...
!> \param mat_global ...
!> \param mat_local ...
!> \param tensor ...
!> \param bs_env Band-structure environment containing GW parameters.
!> \param atom_ranges ...
! **************************************************************************************************
   SUBROUTINE fm_to_local_tensor(fm_global, mat_global, mat_local, tensor, bs_env, atom_ranges)

      TYPE(cp_fm_type)                                   :: fm_global
      TYPE(dbcsr_type)                                   :: mat_global, mat_local
      TYPE(dbt_type)                                     :: tensor
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:, :), OPTIONAL                 :: atom_ranges

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_to_local_tensor'

      INTEGER                                            :: handle
      TYPE(dbt_type)                                     :: tensor_tmp

      CALL timeset(routineN, handle)

      CALL dbt_clear(tensor)
      CALL copy_fm_to_dbcsr(fm_global, mat_global, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(mat_global, bs_env%eps_filter)
      IF (PRESENT(atom_ranges)) THEN
         CALL cp_dbcsr_to_subgroups(mat_global, mat_local, bs_env%para_env, &
                                            bs_env%para_env_tensor%num_pe, atom_ranges)
      ELSE
         CALL cp_dbcsr_to_subgroups(mat_global, mat_local, bs_env%para_env, &
                                            bs_env%para_env_tensor%num_pe)
      END IF
      CALL dbt_create(mat_local, tensor_tmp)
      CALL dbt_copy_matrix_to_tensor(mat_local, tensor_tmp)
      CALL dbt_copy(tensor_tmp, tensor, move_data=.TRUE.)
      CALL dbt_destroy(tensor_tmp)
      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)

      CALL timestop(handle)

   END SUBROUTINE fm_to_local_tensor

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param mat_tensor ...
!> \param mat_global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE local_dbt_to_global_mat(tensor, mat_tensor, mat_global, para_env)

      TYPE(dbt_type)                                     :: tensor
      TYPE(dbcsr_type)                                   :: mat_tensor, mat_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_dbt_to_global_mat'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_copy_tensor_to_matrix(tensor, mat_tensor)
      CALL dbt_clear(tensor)
      ! the next para_env%sync is not mandatory, but it makes the timing output
      ! of local_matrix_to_global_matrix correct
      CALL para_env%sync()
      CALL local_matrix_to_global_matrix(mat_tensor, mat_global, para_env)

      CALL timestop(handle)

   END SUBROUTINE local_dbt_to_global_mat


! **************************************************************************************************
!> \brief ...
!> \param mat_local ...
!> \param mat_global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE local_matrix_to_global_matrix(mat_local, mat_global, para_env)

      TYPE(dbcsr_type)                                   :: mat_local, mat_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_matrix_to_global_matrix'

      INTEGER                                            :: block_size, c, col, col_size, handle, &
                                                            handle1, i_block, imep, o, offset, r, &
                                                            row, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, entry_counter, num_blocks_rec, &
         num_blocks_send, num_entries_rec, num_entries_send, sizes_rec, sizes_send
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(buffer_type), ALLOCATABLE, DIMENSION(:)       :: buffer_rec, buffer_send
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: mat_global_copy

      CALL timeset(routineN, handle)

      CALL timeset("get_coord", handle1)

      CALL dbcsr_create(mat_global_copy, template=mat_global)
      CALL dbcsr_reserve_all_blocks(mat_global_copy)

      CALL dbcsr_set(mat_global, 0.0_dp)
      CALL dbcsr_set(mat_global_copy, 0.0_dp)

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

      CALL dbcsr_iterator_start(iter, mat_local)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         CALL dbcsr_get_stored_coordinates(mat_global, row, col, imep)

         num_entries_send(imep) = num_entries_send(imep) + row_size*col_size
         num_blocks_send(imep) = num_blocks_send(imep) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("comm_size", handle1)

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe - 1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe - 1))

         DO imep = 0, para_env%num_pe - 1

            sizes_send(2*imep) = num_entries_send(imep)
            sizes_send(2*imep + 1) = num_blocks_send(imep)

         END DO

         CALL para_env%alltoall(sizes_send, sizes_rec, 2)

         DO imep = 0, para_env%num_pe - 1
            num_entries_rec(imep) = sizes_rec(2*imep)
            num_blocks_rec(imep) = sizes_rec(2*imep + 1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send(0)
         num_blocks_rec(0) = num_blocks_send(0)

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_buffer", handle1)

      ! allocate data message and corresponding indices
      DO imep = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imep)%msg(num_entries_rec(imep)))
         buffer_rec(imep)%msg = 0.0_dp

         ALLOCATE (buffer_send(imep)%msg(num_entries_send(imep)))
         buffer_send(imep)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imep)%indx(num_blocks_rec(imep), 5))
         buffer_rec(imep)%indx = 0

         ALLOCATE (buffer_send(imep)%indx(num_blocks_send(imep), 5))
         buffer_send(imep)%indx = 0

      END DO

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = 0

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ! fill buffer_send
      CALL dbcsr_iterator_start(iter, mat_local)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         CALL dbcsr_get_stored_coordinates(mat_global, row, col, imep)

         block_size = row_size*col_size

         offset = entry_counter(imep)

         buffer_send(imep)%msg(offset + 1:offset + block_size) = &
            RESHAPE(data_block(1:row_size, 1:col_size), [block_size])

         i_block = block_counter(imep) + 1

         buffer_send(imep)%indx(i_block, 1) = row
         buffer_send(imep)%indx(i_block, 2) = col
         buffer_send(imep)%indx(i_block, 3) = offset

         entry_counter(imep) = entry_counter(imep) + block_size

         block_counter(imep) = block_counter(imep) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("comm_data", handle1)

      CALL exchange_matrix_buffers(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send)

      CALL timestop(handle1)

      CALL timeset("set_blocks", handle1)

      ! fill mat_global_copy
      CALL dbcsr_iterator_start(iter, mat_global_copy)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         DO imep = 0, para_env%num_pe - 1

            DO i_block = 1, num_blocks_rec(imep)

               IF (row == buffer_rec(imep)%indx(i_block, 1) .AND. &
                   col == buffer_rec(imep)%indx(i_block, 2)) THEN

                  offset = buffer_rec(imep)%indx(i_block, 3)

                  r = row_size
                  c = col_size
                  o = offset

                  data_block(1:r, 1:c) = data_block(1:r, 1:c) + &
                                         RESHAPE(buffer_rec(imep)%msg(o + 1:o + r*c), [r, c])

               END IF

            END DO

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_copy(mat_global, mat_global_copy)

      CALL dbcsr_release(mat_global_copy)

      ! remove the blocks which are exactly zero from mat_global
      CALL dbcsr_filter(mat_global, 1.0E-30_dp)

      DO imep = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imep)%msg)
         DEALLOCATE (buffer_send(imep)%msg)
         DEALLOCATE (buffer_rec(imep)%indx)
         DEALLOCATE (buffer_send(imep)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)

      DEALLOCATE (block_counter, entry_counter)

      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE local_matrix_to_global_matrix



! **************************************************************************************************
!> \brief ...
!> \param t_R ...
!> \param fm_R ...
!> \param mat_global ...
!> \param mat_local ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE local_dbt_to_global_fm(t_R, fm_R, mat_global, mat_local, bs_env)
      TYPE(dbt_type), DIMENSION(:)                       :: t_R
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_R
      TYPE(dbcsr_p_type)                                 :: mat_global, mat_local
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_dbt_to_global_fm'

      INTEGER                                            :: handle, i_cell, n_images

      CALL timeset(routineN, handle)

      n_images = SIZE(t_R)

      CPASSERT(n_images == SIZE(fm_R))

      DO i_cell = 1, n_images
         CALL dbcsr_set(mat_global%matrix, 0.0_dp)
         CALL dbcsr_set(mat_local%matrix, 0.0_dp)
         CALL local_dbt_to_global_mat(t_R(i_cell), mat_local%matrix, mat_global%matrix, &
                                      bs_env%para_env)
         CALL copy_dbcsr_to_fm(mat_global%matrix, fm_R(i_cell))
      END DO

      CALL timestop(handle)

   END SUBROUTINE local_dbt_to_global_fm

END MODULE gw_utils_communication

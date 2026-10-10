!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Redistribution of DBCSR blocks and distributed matrix slices.
! **************************************************************************************************
MODULE cp_dbcsr_redistribution
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_filter, dbcsr_finalize, dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_reserve_blocks, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE matrix_communication,            ONLY: exchange_matrix_buffers,&
                                              matrix_buffer_type
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_redistribution'
   PUBLIC :: cp_dbcsr_to_subgroups

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mat_global ...
!> \param mat_local ...
!> \param para_env ...
!> \param num_pe_sub ...
!> \param block_row_ranges ...
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_to_subgroups(mat_global, mat_local, para_env, num_pe_sub, block_row_ranges)
      TYPE(dbcsr_type)                                   :: mat_global, mat_local
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER                                            :: num_pe_sub
      INTEGER, DIMENSION(:, :), OPTIONAL                 :: block_row_ranges

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_to_subgroups'

      INTEGER :: block_counter, block_offset, block_size, col, col_offset, col_size, handle, &
         handle1, i_block, i_mepos, igroup, imep, imep_sub, msg_offset, ngroup, num_blocks, &
         offset, row, row_offset, row_size, total_num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_counter, cols_to_alloc, entry_counter, &
         num_entries_blocks_rec, num_entries_blocks_send, rows_to_alloc, sizes_rec, sizes_send
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(matrix_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset("get_sizes", handle1)

      NULLIFY (data_block)

      ALLOCATE (num_entries_blocks_send(0:2*para_env%num_pe - 1))
      num_entries_blocks_send(:) = 0

      ALLOCATE (num_entries_blocks_rec(0:2*para_env%num_pe - 1))
      num_entries_blocks_rec(:) = 0

      ngroup = para_env%num_pe/num_pe_sub

      CALL dbcsr_iterator_start(iter, mat_global)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imep_sub)

         DO igroup = 0, ngroup - 1

            IF (PRESENT(block_row_ranges)) THEN
               IF (row < block_row_ranges(1, igroup + 1) .OR. row > block_row_ranges(2, igroup + 1)) CYCLE
            END IF
            imep = imep_sub + igroup*num_pe_sub

            num_entries_blocks_send(2*imep) = num_entries_blocks_send(2*imep) + row_size*col_size
            num_entries_blocks_send(2*imep + 1) = num_entries_blocks_send(2*imep + 1) + 1

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("send_sizes_1", handle1)

      total_num_entries = SUM(num_entries_blocks_send)
      CALL para_env%sum(total_num_entries)

      CALL timestop(handle1)

      CALL timeset("send_sizes_2", handle1)

      IF (para_env%num_pe > 1) THEN
         CALL para_env%alltoall(num_entries_blocks_send, num_entries_blocks_rec, 2)
      ELSE
         num_entries_blocks_rec(0:1) = num_entries_blocks_send(0:1)
      END IF

      CALL timestop(handle1)

      CALL timeset("get_data", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO imep = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imep)%msg(num_entries_blocks_rec(2*imep)))
         buffer_rec(imep)%msg = 0.0_dp

         ALLOCATE (buffer_send(imep)%msg(num_entries_blocks_send(2*imep)))
         buffer_send(imep)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imep)%indx(num_entries_blocks_rec(2*imep + 1), 3))
         buffer_rec(imep)%indx = 0

         ALLOCATE (buffer_send(imep)%indx(num_entries_blocks_send(2*imep + 1), 3))
         buffer_send(imep)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ALLOCATE (blk_counter(0:para_env%num_pe - 1))
      blk_counter = 0

      CALL dbcsr_iterator_start(iter, mat_global)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imep_sub)

         DO igroup = 0, ngroup - 1

            IF (PRESENT(block_row_ranges)) THEN
               IF (row < block_row_ranges(1, igroup + 1) .OR. row > block_row_ranges(2, igroup + 1)) CYCLE
            END IF

            imep = imep_sub + igroup*num_pe_sub

            msg_offset = entry_counter(imep)

            block_size = row_size*col_size

            buffer_send(imep)%msg(msg_offset + 1:msg_offset + block_size) = &
               RESHAPE(data_block(1:row_size, 1:col_size), [block_size])

            entry_counter(imep) = entry_counter(imep) + block_size

            blk_counter(imep) = blk_counter(imep) + 1

            block_offset = blk_counter(imep)

            buffer_send(imep)%indx(block_offset, 1) = row
            buffer_send(imep)%indx(block_offset, 2) = col
            buffer_send(imep)%indx(block_offset, 3) = msg_offset

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("send_data", handle1)

      ALLOCATE (sizes_rec(0:para_env%num_pe - 1))
      ALLOCATE (sizes_send(0:para_env%num_pe - 1))

      DO imep = 0, para_env%num_pe - 1
         sizes_send(imep) = num_entries_blocks_send(2*imep)
         sizes_rec(imep) = num_entries_blocks_rec(2*imep)
      END DO

      CALL exchange_matrix_buffers(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send)

      CALL timestop(handle1)

      CALL timeset("reserve_blocks", handle1)

      num_blocks = 0

      ! get the number of blocks, which have to be allocated
      DO imep = 0, para_env%num_pe - 1
         num_blocks = num_blocks + num_entries_blocks_rec(2*imep + 1)
      END DO

      ALLOCATE (rows_to_alloc(num_blocks))
      rows_to_alloc = 0

      ALLOCATE (cols_to_alloc(num_blocks))
      cols_to_alloc = 0

      block_counter = 0

      DO i_mepos = 0, para_env%num_pe - 1

         DO i_block = 1, num_entries_blocks_rec(2*i_mepos + 1)

            block_counter = block_counter + 1

            rows_to_alloc(block_counter) = buffer_rec(i_mepos)%indx(i_block, 1)
            cols_to_alloc(block_counter) = buffer_rec(i_mepos)%indx(i_block, 2)

         END DO

      END DO

      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)
      CALL dbcsr_reserve_blocks(mat_local, rows=rows_to_alloc(:), cols=cols_to_alloc(:))
      CALL dbcsr_finalize(mat_local)
      CALL dbcsr_set(mat_local, 0.0_dp)

      CALL timestop(handle1)

      CALL timeset("fill_mat_local", handle1)

      DO imep = 0, para_env%num_pe - 1
         DO i_block = 1, num_entries_blocks_rec(2*imep + 1)
            row = buffer_rec(imep)%indx(i_block, 1)
            col = buffer_rec(imep)%indx(i_block, 2)
            offset = buffer_rec(imep)%indx(i_block, 3)
            CALL dbcsr_get_block_p(mat_local, row, col, data_block, found=found, &
                                   row_size=row_size, col_size=col_size)
            CPASSERT(found)
            data_block(:, :) = RESHAPE(buffer_rec(imep)%msg(offset + 1:offset + row_size*col_size), &
                                       [row_size, col_size])
         END DO
      END DO

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_to_subgroups

END MODULE cp_dbcsr_redistribution

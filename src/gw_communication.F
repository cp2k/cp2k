!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 08.2023
! **************************************************************************************************
MODULE gw_communication
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_finalize, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_release, &
        dbcsr_reserve_all_blocks, dbcsr_reserve_blocks, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type,&
                                              mp_waitall
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_communication'

   PUBLIC :: global_matrix_to_local_matrix, local_matrix_to_global_matrix

   TYPE buffer_type
      REAL(KIND=dp), DIMENSION(:), POINTER  :: msg
      INTEGER, DIMENSION(:), POINTER  :: sizes
      INTEGER, DIMENSION(:, :), POINTER  :: indx
      INTEGER :: proc
      INTEGER :: msg_req
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mat_global ...
!> \param mat_local ...
!> \param para_env ...
!> \param num_pe_sub ...
! **************************************************************************************************
   SUBROUTINE global_matrix_to_local_matrix(mat_global, mat_local, para_env, num_pe_sub)
      TYPE(dbcsr_type)                                   :: mat_global, mat_local
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: num_pe_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'global_matrix_to_local_matrix'

      INTEGER :: block_counter, block_offset, block_size, col, col_from_buffer, col_offset, &
         col_size, handle, handle1, i_block, i_entry, i_mepos, igroup, imep, imep_sub, msg_offset, &
         nblkrows_total, ngroup, nmo, num_blocks, offset, row, row_from_buffer, row_offset, &
         row_size, total_num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_counter, cols_to_alloc, entry_counter, &
         num_entries_blocks_rec, num_entries_blocks_send, row_block_from_index, rows_to_alloc, &
         sizes_rec, sizes_send
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_offset, row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(buffer_type), ALLOCATABLE, DIMENSION(:)       :: buffer_rec, buffer_send
      TYPE(dbcsr_iterator_type)                          :: iter

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

         ! JW here additional check is needed for iatom and jatom range
         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imep_sub)

         DO igroup = 0, ngroup - 1

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

         ! JW check here whether block is in iatom/jatom range
         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imep_sub)

         DO igroup = 0, ngroup - 1

            imep = imep_sub + igroup*num_pe_sub

            msg_offset = entry_counter(imep)

            block_size = row_size*col_size

            buffer_send(imep)%msg(msg_offset + 1:msg_offset + block_size) = &
               RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

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

      CALL communicate_buffer(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send)

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

            IF (i_entry >= row_blk_offset(i_block) .AND. &
                i_entry <= row_blk_offset(i_block) + row_blk_size(i_block) - 1) THEN

               row_block_from_index(i_entry) = i_block

            END IF

         END DO
      END DO

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

      CALL dbcsr_iterator_start(iter, mat_local)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         DO imep = 0, para_env%num_pe - 1

            DO i_block = 1, num_entries_blocks_rec(2*imep + 1)

               row_from_buffer = buffer_rec(imep)%indx(i_block, 1)
               col_from_buffer = buffer_rec(imep)%indx(i_block, 2)
               offset = buffer_rec(imep)%indx(i_block, 3)

               IF (row == row_from_buffer .AND. col == col_from_buffer) THEN

                  data_block(1:row_size, 1:col_size) = &
                     RESHAPE(buffer_rec(imep)%msg(offset + 1:offset + row_size*col_size), &
                             (/row_size, col_size/))

               END IF

            END DO

         END DO

      END DO ! blocks

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      DO imep = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imep)%msg)
         DEALLOCATE (buffer_rec(imep)%indx)
         DEALLOCATE (buffer_send(imep)%msg)
         DEALLOCATE (buffer_send(imep)%indx)
      END DO

      CALL timestop(handle)

   END SUBROUTINE global_matrix_to_local_matrix

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param num_entries_rec ...
!> \param num_entries_send ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param do_indx ...
!> \param do_msg ...
! **************************************************************************************************
   SUBROUTINE communicate_buffer(para_env, num_entries_rec, num_entries_send, &
                                 buffer_rec, buffer_send, do_indx, do_msg)

      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_entries_rec, num_entries_send
      TYPE(buffer_type), ALLOCATABLE, DIMENSION(:)       :: buffer_rec, buffer_send
      LOGICAL, OPTIONAL                                  :: do_indx, do_msg

      CHARACTER(LEN=*), PARAMETER :: routineN = 'communicate_buffer'

      INTEGER                                            :: handle, imep, rec_counter, send_counter
      LOGICAL                                            :: my_do_indx, my_do_msg
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req

      CALL timeset(routineN, handle)

      NULLIFY (req)
      ALLOCATE (req(1:para_env%num_pe, 4))

      my_do_indx = .TRUE.
      IF (PRESENT(do_indx)) my_do_indx = do_indx
      my_do_msg = .TRUE.
      IF (PRESENT(do_msg)) my_do_msg = do_msg

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imep = 0, para_env%num_pe - 1
            IF (num_entries_rec(imep) > 0) THEN
               rec_counter = rec_counter + 1
               IF (my_do_indx) THEN
                  CALL para_env%irecv(buffer_rec(imep)%indx, imep, req(rec_counter, 3), tag=4)
               END IF
               IF (my_do_msg) THEN
                  CALL para_env%irecv(buffer_rec(imep)%msg, imep, req(rec_counter, 4), tag=7)
               END IF
            END IF
         END DO

         DO imep = 0, para_env%num_pe - 1
            IF (num_entries_send(imep) > 0) THEN
               send_counter = send_counter + 1
               IF (my_do_indx) THEN
                  CALL para_env%isend(buffer_send(imep)%indx, imep, req(send_counter, 1), tag=4)
               END IF
               IF (my_do_msg) THEN
                  CALL para_env%isend(buffer_send(imep)%msg, imep, req(send_counter, 2), tag=7)
               END IF
            END IF
         END DO

         IF (my_do_indx) THEN
            CALL mp_waitall(req(1:send_counter, 1))
            CALL mp_waitall(req(1:rec_counter, 3))
         END IF

         IF (my_do_msg) THEN
            CALL mp_waitall(req(1:send_counter, 2))
            CALL mp_waitall(req(1:rec_counter, 4))
         END IF

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      DEALLOCATE (req)

      CALL timestop(handle)

   END SUBROUTINE communicate_buffer

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
                                                            rec_counter, row, row_size, &
                                                            send_counter
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, entry_counter, num_blocks_rec, &
         num_blocks_send, num_entries_rec, num_entries_send, sizes_rec, sizes_send
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(buffer_type), ALLOCATABLE, DIMENSION(:)       :: buffer_rec, buffer_send
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: mat_global_copy
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req

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
            RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

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

      NULLIFY (req)
      ALLOCATE (req(1:para_env%num_pe, 4))

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imep = 0, para_env%num_pe - 1
            IF (num_entries_rec(imep) > 0) THEN
               rec_counter = rec_counter + 1
               CALL para_env%irecv(buffer_rec(imep)%indx, imep, req(rec_counter, 3), tag=4)
            END IF
            IF (num_entries_rec(imep) > 0) THEN
               CALL para_env%irecv(buffer_rec(imep)%msg, imep, req(rec_counter, 4), tag=7)
            END IF
         END DO

         DO imep = 0, para_env%num_pe - 1
            IF (num_entries_send(imep) > 0) THEN
               send_counter = send_counter + 1
               CALL para_env%isend(buffer_send(imep)%indx, imep, req(send_counter, 1), tag=4)
            END IF
            IF (num_entries_send(imep) > 0) THEN
               CALL para_env%isend(buffer_send(imep)%msg, imep, req(send_counter, 2), tag=7)
            END IF
         END DO

         CALL mp_waitall(req(1:send_counter, 1:2))
         CALL mp_waitall(req(1:rec_counter, 3:4))

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

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
                                         RESHAPE(buffer_rec(imep)%msg(o + 1:o + r*c), (/r, c/))

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

      DEALLOCATE (req)

      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE local_matrix_to_global_matrix

END MODULE gw_communication

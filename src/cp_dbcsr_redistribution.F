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
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_clear, dbcsr_filter, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_readonly_start, dbcsr_iterator_start, dbcsr_iterator_stop, &
        dbcsr_iterator_type, dbcsr_p_type, dbcsr_reserve_blocks, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE matrix_communication,            ONLY: allocate_matrix_buffers,&
                                              exchange_matrix_buffers,&
                                              matrix_buffer_type
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_redistribution'
   PUBLIC :: cp_dbcsr_to_subgroups, cp_dbcsr_from_subgroups, cp_array_to_dbcsr_slices

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

      INTEGER :: block_offset, block_size, col, col_size, handle, handle2, i_block, igroup, imep, &
         imep_sub, msg_offset, ngroup, offset, row, row_size
      INTEGER, ALLOCATABLE                               :: counts_rec(:, :), counts_send(:, :)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_counter, entry_counter
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(matrix_buffer_type)                           :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset(routineN//"_count", handle2)

      NULLIFY (data_block)

      ALLOCATE (counts_send(2, 0:para_env%num_pe - 1), counts_rec(2, 0:para_env%num_pe - 1))
      counts_send = 0

      ngroup = para_env%num_pe/num_pe_sub

      CALL dbcsr_iterator_start(iter, mat_global)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imep_sub)

         DO igroup = 0, ngroup - 1

            IF (PRESENT(block_row_ranges)) THEN
               IF (row < block_row_ranges(1, igroup + 1) .OR. row > block_row_ranges(2, igroup + 1)) CYCLE
            END IF
            imep = imep_sub + igroup*num_pe_sub

            counts_send(1, imep) = counts_send(1, imep) + row_size*col_size
            counts_send(2, imep) = counts_send(2, imep) + 1

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle2)

      CALL timeset(routineN//"_exchange_counts", handle2)

      IF (para_env%num_pe > 1) THEN
         CALL para_env%alltoall(counts_send, counts_rec, 2)
      ELSE
         counts_rec(:, :) = counts_send
      END IF

      CALL timestop(handle2)

      CALL timeset(routineN//"_pack", handle2)

      CALL allocate_matrix_buffers(counts_send(1, :), counts_rec(1, :), 3, buffer_send, buffer_rec, &
                                   counts_send(2, :), counts_rec(2, :))

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = buffer_send%msg_offset(:para_env%num_pe - 1)

      ALLOCATE (blk_counter(0:para_env%num_pe - 1))
      blk_counter(:) = buffer_send%indx_offset(:para_env%num_pe - 1)

      CALL dbcsr_iterator_start(iter, mat_global)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         CALL dbcsr_get_stored_coordinates(mat_local, row, col, imep_sub)

         DO igroup = 0, ngroup - 1

            IF (PRESENT(block_row_ranges)) THEN
               IF (row < block_row_ranges(1, igroup + 1) .OR. row > block_row_ranges(2, igroup + 1)) CYCLE
            END IF

            imep = imep_sub + igroup*num_pe_sub

            msg_offset = entry_counter(imep)

            block_size = row_size*col_size

            buffer_send%msg(msg_offset + 1:msg_offset + block_size) = &
               RESHAPE(data_block(1:row_size, 1:col_size), [block_size])

            entry_counter(imep) = entry_counter(imep) + block_size

            blk_counter(imep) = blk_counter(imep) + 1

            block_offset = blk_counter(imep)

            ! The payload offset in each record is relative to its peer's segment.
            buffer_send%indx(:, block_offset) = [row, col, msg_offset - buffer_send%msg_offset(imep)]

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle2)

      CALL timeset(routineN//"_exchange_data", handle2)

      CALL exchange_matrix_buffers(para_env, buffer_rec, buffer_send)

      CALL timestop(handle2)

      CALL timeset(routineN//"_reserve", handle2)

      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)
      CALL dbcsr_reserve_blocks(mat_local, rows=buffer_rec%indx(1, :), cols=buffer_rec%indx(2, :))
      CALL dbcsr_finalize(mat_local)
      CALL dbcsr_set(mat_local, 0.0_dp)

      CALL timestop(handle2)

      CALL timeset(routineN//"_unpack", handle2)

      DO imep = 0, para_env%num_pe - 1
         DO i_block = buffer_rec%indx_offset(imep) + 1, buffer_rec%indx_offset(imep + 1)
            row = buffer_rec%indx(1, i_block)
            col = buffer_rec%indx(2, i_block)
            offset = buffer_rec%msg_offset(imep) + buffer_rec%indx(3, i_block)
            CALL dbcsr_get_block_p(mat_local, row, col, data_block, found=found, &
                                   row_size=row_size, col_size=col_size)
            CPASSERT(found)
            data_block(:, :) = RESHAPE(buffer_rec%msg(offset + 1:offset + row_size*col_size), &
                                       [row_size, col_size])
         END DO
      END DO

      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_to_subgroups

! **************************************************************************************************
!> \brief Sum subgroup blocks into the global layout, or copy mapped blocks.
!> \param mat_local ...
!> \param mat_global ...
!> \param para_env communicator covering all source and destination ranks
!> \param row_mapping optional source-to-destination block-row permutation
!> \param copy_only copy without summation, filtering, or consuming the source (default: false)
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_from_subgroups(mat_local, mat_global, para_env, row_mapping, copy_only)

      TYPE(dbcsr_type)                                   :: mat_local, mat_global
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: row_mapping
      LOGICAL, INTENT(IN), OPTIONAL                      :: copy_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_from_subgroups'

      INTEGER                                            :: block_size, col, col_size, handle, &
                                                            handle2, i_block, imep, mapped_row, &
                                                            offset, row, row_size
      INTEGER, ALLOCATABLE                               :: counts_rec(:, :), counts_send(:, :)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_counter, entry_counter
      LOGICAL                                            :: found, my_copy_only
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(matrix_buffer_type)                           :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)
      my_copy_only = .FALSE.
      IF (PRESENT(copy_only)) my_copy_only = copy_only

      CALL timeset(routineN//"_prepare", handle2)

      IF (.NOT. my_copy_only) THEN
         CALL dbcsr_clear(mat_global)
         CALL dbcsr_reserve_all_blocks(mat_global)
      END IF

      ALLOCATE (counts_send(2, 0:para_env%num_pe - 1), counts_rec(2, 0:para_env%num_pe - 1))
      counts_send = 0

      CALL dbcsr_iterator_readonly_start(iter, mat_local)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         mapped_row = row
         IF (PRESENT(row_mapping)) mapped_row = row_mapping(row)
         CALL dbcsr_get_stored_coordinates(mat_global, mapped_row, col, imep)

         counts_send(1, imep) = counts_send(1, imep) + row_size*col_size
         counts_send(2, imep) = counts_send(2, imep) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle2)

      CALL timeset(routineN//"_exchange_counts", handle2)

      IF (para_env%num_pe > 1) THEN
         CALL para_env%alltoall(counts_send, counts_rec, 2)
      ELSE
         counts_rec(:, :) = counts_send
      END IF

      CALL timestop(handle2)

      CALL timeset(routineN//"_pack", handle2)

      CALL allocate_matrix_buffers(counts_send(1, :), counts_rec(1, :), 3, buffer_send, buffer_rec, &
                                   counts_send(2, :), counts_rec(2, :))

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = buffer_send%indx_offset(:para_env%num_pe - 1)

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = buffer_send%msg_offset(:para_env%num_pe - 1)

      ! fill buffer_send
      CALL dbcsr_iterator_readonly_start(iter, mat_local)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         mapped_row = row
         IF (PRESENT(row_mapping)) mapped_row = row_mapping(row)
         CALL dbcsr_get_stored_coordinates(mat_global, mapped_row, col, imep)

         block_size = row_size*col_size

         offset = entry_counter(imep)

         buffer_send%msg(offset + 1:offset + block_size) = &
            RESHAPE(data_block(1:row_size, 1:col_size), [block_size])

         i_block = block_counter(imep) + 1

         ! The receiver adds its own segment start to this peer-relative offset.
         buffer_send%indx(:, i_block) = [mapped_row, col, offset - buffer_send%msg_offset(imep)]

         entry_counter(imep) = entry_counter(imep) + block_size

         block_counter(imep) = block_counter(imep) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle2)

      CALL timeset(routineN//"_exchange_data", handle2)

      CALL exchange_matrix_buffers(para_env, buffer_rec, buffer_send)

      CALL timestop(handle2)

      CALL timeset(routineN//"_unpack", handle2)

      ! Keep the rank and receive-buffer order of the block-wise accumulation.
      DO imep = 0, para_env%num_pe - 1
         DO i_block = buffer_rec%indx_offset(imep) + 1, buffer_rec%indx_offset(imep + 1)
            row = buffer_rec%indx(1, i_block)
            col = buffer_rec%indx(2, i_block)
            offset = buffer_rec%msg_offset(imep) + buffer_rec%indx(3, i_block)
            CALL dbcsr_get_block_p(mat_global, row, col, data_block, found=found, &
                                   row_size=row_size, col_size=col_size)
            CPASSERT(found)
            IF (my_copy_only) THEN
               data_block(:, :) = RESHAPE(buffer_rec%msg(offset + 1:offset + row_size*col_size), &
                                          [row_size, col_size])
            ELSE
               data_block(:, :) = data_block + &
                                  RESHAPE(buffer_rec%msg(offset + 1:offset + row_size*col_size), &
                                          [row_size, col_size])
            END IF
         END DO
      END DO

      IF (.NOT. my_copy_only) THEN
         ! Remove numerically zero blocks from the summed matrix.
         CALL dbcsr_filter(mat_global, 1.0E-30_dp)
      END IF

      IF (.NOT. my_copy_only) THEN
         CALL dbcsr_set(mat_local, 0.0_dp)
         CALL dbcsr_filter(mat_local, 1.0_dp)
      END IF

      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_from_subgroups

! **************************************************************************************************
!> \brief Consume a distributed flattened array and create a DBCSR matrix per slice.
!> \param array local flattened rows; all ranks have the same number of slices; deallocated on return
!> \param nrows number of rows in each output matrix
!> \param ncols number of columns in each output matrix
!> \param para_env_sub communicator for the output matrices
!> \param first_row first flattened row held locally
!> \param last_row last flattened row held locally
!> \param template DBCSR template defining output distribution and block sizes
!> \param matrices preallocated array of output matrices
!> \param matrix_ptrs preallocated array of output pointers; pointees are allocated here
!> \note Supply exactly one output representation. Flattened rows use (i-1)*ncols+j.
! **************************************************************************************************
   SUBROUTINE cp_array_to_dbcsr_slices(array, nrows, ncols, para_env_sub, &
                                       first_row, last_row, template, matrices, matrix_ptrs)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: array
      INTEGER, INTENT(IN)                                :: nrows, ncols
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      INTEGER, INTENT(IN)                                :: first_row, last_row
      TYPE(dbcsr_type), INTENT(INOUT)                    :: template
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: matrices
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: matrix_ptrs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_array_to_dbcsr_slices'

      INTEGER                                            :: flat_row, handle, i, ientry, islice, j, &
                                                            nprocs, pass, peer
      INTEGER, ALLOCATABLE                               :: counter(:), recv_count(:), &
                                                            send_count(:), source_rows(:)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_ia
      TYPE(dbcsr_type), POINTER                          :: output_matrix
      TYPE(matrix_buffer_type)                           :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)
      CPASSERT(PRESENT(matrices) .NEQV. PRESENT(matrix_ptrs))
      NULLIFY (blacs_env, fm_struct)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_sub)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrows, &
                               ncol_global=ncols, para_env=para_env_sub)
      CALL cp_fm_create(fm_ia, fm_struct, name="fm_ia")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_blacs_env_release(blacs_env)

      nprocs = para_env_sub%num_pe
      ALLOCATE (send_count(0:nprocs - 1), recv_count(0:nprocs - 1), counter(0:nprocs - 1))
      send_count = 0

      ! Build the coordinate and packing maps once for all slices.
      DO pass = 1, 2
         IF (pass == 2) THEN
            CALL para_env_sub%alltoall(send_count, recv_count, 1)
            CALL allocate_matrix_buffers(send_count, recv_count, 2, buffer_send, buffer_rec)
            ALLOCATE (source_rows(SIZE(buffer_send%msg)))
            counter(:) = buffer_send%msg_offset(:nprocs - 1)
         END IF
         DO flat_row = first_row, last_row
            i = (flat_row - 1)/ncols + 1
            j = MOD(flat_row - 1, ncols) + 1
            peer = fm_ia%matrix_struct%context%blacs2mpi( &
                   fm_ia%matrix_struct%g2p_row(i), fm_ia%matrix_struct%g2p_col(j))
            IF (pass == 1) THEN
               send_count(peer) = send_count(peer) + 1
            ELSE
               counter(peer) = counter(peer) + 1
               ientry = counter(peer)
               buffer_send%indx(:, ientry) = [i, j]
               source_rows(ientry) = flat_row - first_row + 1
            END IF
         END DO
      END DO

      CALL exchange_matrix_buffers(para_env_sub, buffer_rec, buffer_send, do_msg=.FALSE.)
      DO ientry = 1, SIZE(buffer_rec%msg)
         buffer_rec%indx(1, ientry) = fm_ia%matrix_struct%g2l_row(buffer_rec%indx(1, ientry))
         buffer_rec%indx(2, ientry) = fm_ia%matrix_struct%g2l_col(buffer_rec%indx(2, ientry))
      END DO

      DO islice = 1, SIZE(array, 2)
         DO ientry = 1, SIZE(buffer_send%msg)
            buffer_send%msg(ientry) = array(source_rows(ientry), islice)
         END DO
         CALL exchange_matrix_buffers(para_env_sub, buffer_rec, buffer_send, do_indx=.FALSE.)
         CALL cp_fm_set_all(fm_ia, 0.0_dp)
         DO ientry = 1, SIZE(buffer_rec%msg)
            i = buffer_rec%indx(1, ientry)
            j = buffer_rec%indx(2, ientry)
            fm_ia%local_data(i, j) = buffer_rec%msg(ientry)
         END DO

         IF (PRESENT(matrices)) THEN
            output_matrix => matrices(islice)
         ELSE
            ALLOCATE (matrix_ptrs(islice)%matrix)
            output_matrix => matrix_ptrs(islice)%matrix
         END IF
         CALL cp_dbcsr_m_by_n_from_template(output_matrix, template=template, &
                                            m=nrows, n=ncols, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_ia, output_matrix, keep_sparsity=.FALSE.)
      END DO

      DEALLOCATE (array)
      CALL cp_fm_release(fm_ia)
      CALL timestop(handle)
   END SUBROUTINE cp_array_to_dbcsr_slices

END MODULE cp_dbcsr_redistribution

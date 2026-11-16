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
        dbcsr_clear, dbcsr_filter, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_get_readonly_block_p, dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_readonly_start, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_reserve_blocks, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_m_by_n_from_template
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE matrix_communication,            ONLY: allocate_matrix_buffers,&
                                              exchange_matrix_buffers,&
                                              matrix_buffer_type
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_redistribution'
   PUBLIC :: cp_dbcsr_to_subgroups, cp_dbcsr_from_subgroups, cp_array_to_dbcsr_slices, cp_array_to_dbcsr_blocks

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
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT), TARGET  :: array(:, :)
      INTEGER, INTENT(IN)                                :: nrows, ncols
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      INTEGER, INTENT(IN)                                :: first_row, last_row
      TYPE(dbcsr_type), INTENT(INOUT)                    :: template
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL, TARGET  :: matrices(:)
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: matrix_ptrs(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_array_to_dbcsr_slices'

      INTEGER                                            :: handle, i, islice, j, nblock, nslice
      INTEGER, ALLOCATABLE                               :: blocks(:, :)
      INTEGER, POINTER                                   :: local_cols(:), local_rows(:)
      REAL(KIND=dp), POINTER                             :: values(:, :, :)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: outputs(:)

      CALL timeset(routineN, handle)
      CPASSERT(PRESENT(matrices) .NEQV. PRESENT(matrix_ptrs))
      nslice = SIZE(array, 2)
      ALLOCATE (outputs(nslice))
      DO islice = 1, nslice
         IF (PRESENT(matrices)) THEN
            outputs(islice)%matrix => matrices(islice)
         ELSE
            ALLOCATE (matrix_ptrs(islice)%matrix)
            outputs(islice)%matrix => matrix_ptrs(islice)%matrix
         END IF
         CALL cp_dbcsr_m_by_n_from_template(outputs(islice)%matrix, template=template, &
                                            m=nrows, n=ncols, sym=dbcsr_type_no_symmetry)
         CALL dbcsr_reserve_all_blocks(outputs(islice)%matrix)
         CALL dbcsr_finalize(outputs(islice)%matrix)
      END DO
      IF (nslice > 0) THEN
         CALL dbcsr_get_info(outputs(1)%matrix, local_rows=local_rows, local_cols=local_cols)
         ALLOCATE (blocks(3, SIZE(local_rows)*SIZE(local_cols)*nslice))
         nblock = 0
         DO islice = 1, nslice
            DO j = 1, SIZE(local_cols)
               DO i = 1, SIZE(local_rows)
                  nblock = nblock + 1
                  blocks(:, nblock) = [local_rows(i), local_cols(j), islice]
               END DO
            END DO
         END DO
         ! A rank-remapped view preserves the row-major flattened ia convention without copying.
         values(1:SIZE(array, 1), 1:1, 1:nslice) => array
         CALL cp_array_to_dbcsr_blocks(values, first_row, last_row, para_env_sub, outputs, blocks, &
                                       row_major=.TRUE.)
      END IF
      DEALLOCATE (array)
      CALL timestop(handle)
   END SUBROUTINE cp_array_to_dbcsr_slices

! **************************************************************************************************
!> \brief Fill existing DBCSR blocks directly from distributed flattened array slices.
!>        Requests split contiguous block fragments at source-rank boundaries. Sources combine
!>        components before replying, so each target element transfers only one real value.
!>        Slice tiles bound request/reply storage; existing sparsity is preserved.
!> \param array local flattened entries, components, slices; borrowed, never deallocated
!> \param first_entry first flattened entry held locally
!> \param last_entry last flattened entry held locally; may precede first_entry on empty ranks
!> \param comm communicator covering sources and output matrices, in any source-range order
!> \param matrices existing output slices with identical dimensions and block sizes
!> \param blocks requested (block row, block column, slice) columns; absent local blocks are skipped
!> \param factors optional integer component multiplicities per requested block; default one component
!> \param row_major flatten with column index fastest; default Fortran column-major ordering
!> \note Source ranges must partition the whole flattened matrix. Requested blocks must be unique.
!>       Unrequested output entries are zeroed. The 64 MiB workspace target includes the borrowed
!>       array and descriptors; an irreducible single-slice exchange may exceed this target.
! **************************************************************************************************
   SUBROUTINE cp_array_to_dbcsr_blocks(array, first_entry, last_entry, comm, matrices, blocks, factors, row_major)
      REAL(KIND=dp), INTENT(IN)                          :: array(:, :, :)
      INTEGER, INTENT(IN)                                :: first_entry, last_entry, blocks(:, :)

      CLASS(mp_comm_type), INTENT(IN)                    :: comm
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: matrices(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: factors(:, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: row_major

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_array_to_dbcsr_blocks'

      INTEGER(KIND=int_8), PARAMETER                     :: max_bytes = 64_int_8*1024**2

      INTEGER :: a, b, col, component, element, first, handle, i, image, last, length, lo, &
                 lower, middle, nc, nfast, nprocs, pass, peer, pos, rank, record, row, upper, width
      INTEGER(KIND=int_8)                                :: fixed_bytes, nel
      INTEGER, ALLOCATABLE                               :: counter(:), counts_rec(:, :), &
                                                            counts_send(:, :), ends(:), order(:), &
                                                            ranges(:, :), targets(:, :)
      INTEGER, POINTER                                   :: fast_offset(:), fast_size(:), &
                                                            slow_offset(:), slow_size(:)
      LOGICAL                                            :: found, row_order
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(matrix_buffer_type)                           :: requests_rec, requests_send

      CALL timeset(routineN, handle)
      IF (SIZE(matrices) == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF
      row_order = .FALSE.
      IF (PRESENT(row_major)) row_order = row_major
      CALL dbcsr_get_info(matrices(1)%matrix, row_blk_size=fast_size, col_blk_size=slow_size, &
                          row_blk_offset=fast_offset, col_blk_offset=slow_offset)
      IF (row_order) THEN
         CALL dbcsr_get_info(matrices(1)%matrix, row_blk_size=slow_size, col_blk_size=fast_size, &
                             row_blk_offset=slow_offset, col_blk_offset=fast_offset)
      END IF
      nfast = SUM(fast_size)
      nel = INT(nfast, int_8)*SUM(slow_size)
      nc = SIZE(array, 2)
      CPASSERT(nel > 0 .AND. nel < HUGE(element))
      CPASSERT(nc >= 1)
      CPASSERT(SIZE(array, 1) == MAX(0, last_entry - first_entry + 1))
      CPASSERT(SIZE(array, 3) == SIZE(matrices) .AND. SIZE(blocks, 1) == 3)
      IF (PRESENT(factors)) THEN
         CPASSERT(SIZE(factors, 1) == nc .AND. SIZE(factors, 2) == SIZE(blocks, 2))
      ELSE
         CPASSERT(nc == 1)
      END IF
      nprocs = comm%num_pe
      ALLOCATE (ends(nprocs), order(nprocs), ranges(2, nprocs), counts_send(2, 0:nprocs - 1), &
                counts_rec(2, 0:nprocs - 1), counter(0:nprocs - 1))
      CALL comm%allgather([first_entry, last_entry], ranges)
      ends(:) = ranges(2, :)
      WHERE (ranges(2, :) < ranges(1, :)) ends = HUGE(element)
      CALL sort(ends, nprocs, order)
      ! Validate coverage, independently of communicator rank ordering and empty ranges.
      element = 0
      DO rank = 1, nprocs
         IF (ends(rank) == HUGE(element)) EXIT
         CPASSERT(ranges(1, order(rank)) == element + 1)
         element = ends(rank)
      END DO
      CPASSERT(INT(element, int_8) == nel)
      fixed_bytes = 8_int_8*SIZE(array, KIND=int_8) + 4_int_8*SIZE(blocks, KIND=int_8) + 64_int_8*nprocs
      IF (PRESENT(factors)) fixed_bytes = fixed_bytes + 4_int_8*SIZE(factors, KIND=int_8)
      width = INT(MIN(INT(SIZE(matrices), int_8), &
                      MAX(1_int_8, (max_bytes - fixed_bytes)/((80_int_8 + 8*nc)*nel))))
      CALL comm%min(width)
      DO image = 1, SIZE(matrices)
         CALL dbcsr_set(matrices(image)%matrix, 0.0_dp)
      END DO
      DO first = 1, SIZE(matrices), width
         last = MIN(SIZE(matrices), first + width - 1)
         counts_send(:, :) = 0
         DO pass = 1, 2
            IF (pass == 2) THEN
               CALL comm%alltoall(counts_send, counts_rec, 2)
               CALL allocate_matrix_buffers(counts_send(1, :), counts_rec(1, :), 3 + nc, requests_send, requests_rec, &
                                            send_blocks=counts_send(2, :), recv_blocks=counts_rec(2, :))
               ALLOCATE (targets(3, SIZE(requests_send%indx, 2)))
               counter(:) = requests_send%indx_offset(:nprocs - 1)
            END IF
            DO i = 1, SIZE(blocks, 2)
               image = blocks(3, i)
               IF (image < first .OR. image > last) CYCLE
               CALL dbcsr_get_readonly_block_p(matrices(image)%matrix, blocks(1, i), blocks(2, i), block, found=found)
               IF (.NOT. found) CYCLE
               a = blocks(1, i)
               b = blocks(2, i)
               IF (row_order) THEN
                  a = blocks(2, i)
                  b = blocks(1, i)
               END IF
               DO col = 1, slow_size(b)
                  row = 1
                  DO WHILE (row <= fast_size(a))
                     element = (slow_offset(b) + col - 2)*nfast + fast_offset(a) + row - 1
                     ! Find the first nonempty source range whose upper bound covers this entry.
                     lower = 1
                     upper = nprocs
                     DO WHILE (lower < upper)
                        middle = (lower + upper)/2
                        IF (ends(middle) < element) THEN
                           lower = middle + 1
                        ELSE
                           upper = middle
                        END IF
                     END DO
                     peer = order(lower) - 1
                     length = MIN(fast_size(a) - row + 1, ends(lower) - element + 1)
                     IF (pass == 1) THEN
                        counts_send(:, peer) = counts_send(:, peer) + [length, 1]
                     ELSE
                        counter(peer) = counter(peer) + 1
                        record = counter(peer)
                        requests_send%indx(1:3, record) = [element, length, image]
                        requests_send%indx(4:, record) = 1
                        IF (PRESENT(factors)) requests_send%indx(4:, record) = factors(:, i)
                        targets(:, record) = [i, row, col]
                     END IF
                     row = row + length
                  END DO
               END DO
            END DO
         END DO
         CALL exchange_matrix_buffers(comm, requests_rec, requests_send, do_msg=.FALSE.)
         DO peer = 0, nprocs - 1
            pos = requests_rec%msg_offset(peer)
            DO record = requests_rec%indx_offset(peer) + 1, requests_rec%indx_offset(peer + 1)
               lo = requests_rec%indx(1, record) - first_entry + 1
               length = requests_rec%indx(2, record)
               image = requests_rec%indx(3, record)
               CPASSERT(lo >= 1 .AND. lo + length - 1 <= SIZE(array, 1))
               requests_rec%msg(pos + 1:pos + length) = &
                  REAL(requests_rec%indx(4, record), dp)*array(lo:lo + length - 1, 1, image)
               DO component = 2, nc
                  requests_rec%msg(pos + 1:pos + length) = requests_rec%msg(pos + 1:pos + length) + &
                                      REAL(requests_rec%indx(3 + component, record), dp)*array(lo:lo + length - 1, component, image)
               END DO
               pos = pos + length
            END DO
         END DO
         CALL exchange_matrix_buffers(comm, requests_send, requests_rec, do_indx=.FALSE.)
         DO peer = 0, nprocs - 1
            pos = requests_send%msg_offset(peer)
            DO record = requests_send%indx_offset(peer) + 1, requests_send%indx_offset(peer + 1)
               i = targets(1, record)
               row = targets(2, record)
               col = targets(3, record)
               length = requests_send%indx(2, record)
               CALL dbcsr_get_block_p(matrices(blocks(3, i))%matrix, blocks(1, i), blocks(2, i), block, found=found)
               CPASSERT(found)
               IF (row_order) THEN
                  block(col, row:row + length - 1) = requests_send%msg(pos + 1:pos + length)
               ELSE
                  block(row:row + length - 1, col) = requests_send%msg(pos + 1:pos + length)
               END IF
               pos = pos + length
            END DO
         END DO
         DEALLOCATE (targets)
      END DO
      CALL timestop(handle)
   END SUBROUTINE cp_array_to_dbcsr_blocks

END MODULE cp_dbcsr_redistribution

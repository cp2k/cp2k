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
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE group_dist_types,                ONLY: get_group_dist,&
                                              group_dist_d1_type
   USE kinds,                           ONLY: dp
   USE matrix_communication,            ONLY: exchange_matrix_buffers,&
                                              matrix_buffer_type
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_null,&
                                              mp_request_type,&
                                              mp_waitall
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_redistribution'
   PUBLIC :: cp_dbcsr_to_subgroups, cp_dbcsr_from_subgroups, cp_array_to_dbcsr_slices

   TYPE index_map
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: map
   END TYPE index_map

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
                                                            handle1, i_block, imep, mapped_row, &
                                                            offset, row, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, entry_counter, num_blocks_rec, &
         num_blocks_send, num_entries_rec, num_entries_send, sizes_rec, sizes_send
      LOGICAL                                            :: found, my_copy_only
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(matrix_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)
      my_copy_only = .FALSE.
      IF (PRESENT(copy_only)) my_copy_only = copy_only

      CALL timeset("get_coord", handle1)

      IF (.NOT. my_copy_only) THEN
         CALL dbcsr_clear(mat_global)
         CALL dbcsr_reserve_all_blocks(mat_global)
      END IF

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

      CALL dbcsr_iterator_readonly_start(iter, mat_local)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         mapped_row = row
         IF (PRESENT(row_mapping)) mapped_row = row_mapping(row)
         CALL dbcsr_get_stored_coordinates(mat_global, mapped_row, col, imep)

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

         ALLOCATE (buffer_rec(imep)%indx(num_blocks_rec(imep), 3))
         buffer_rec(imep)%indx = 0

         ALLOCATE (buffer_send(imep)%indx(num_blocks_send(imep), 3))
         buffer_send(imep)%indx = 0

      END DO

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = 0

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

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

         buffer_send(imep)%msg(offset + 1:offset + block_size) = &
            RESHAPE(data_block(1:row_size, 1:col_size), [block_size])

         i_block = block_counter(imep) + 1

         buffer_send(imep)%indx(i_block, 1) = mapped_row
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

      ! Keep the rank and receive-buffer order of the block-wise accumulation.
      DO imep = 0, para_env%num_pe - 1
         DO i_block = 1, num_blocks_rec(imep)
            row = buffer_rec(imep)%indx(i_block, 1)
            col = buffer_rec(imep)%indx(i_block, 2)
            offset = buffer_rec(imep)%indx(i_block, 3)
            CALL dbcsr_get_block_p(mat_global, row, col, data_block, found=found, &
                                   row_size=row_size, col_size=col_size)
            CPASSERT(found)
            IF (my_copy_only) THEN
               data_block(:, :) = RESHAPE(buffer_rec(imep)%msg(offset + 1:offset + row_size*col_size), &
                                          [row_size, col_size])
            ELSE
               data_block(:, :) = data_block + &
                                  RESHAPE(buffer_rec(imep)%msg(offset + 1:offset + row_size*col_size), &
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

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_from_subgroups

! **************************************************************************************************
!> \brief Consume a distributed flattened array and create a DBCSR matrix per slice.
!> \param array local rows of the flattened matrices; deallocated on return
!> \param nrows number of rows in each output matrix
!> \param ncols number of columns in each output matrix
!> \param para_env_sub communicator for the output matrices
!> \param first_row first flattened row held locally
!> \param last_row last flattened row held locally
!> \param nslices number of local matrix slices
!> \param row_distribution distribution of flattened rows over the communicator
!> \param template DBCSR template defining output distribution and block sizes
!> \param matrices preallocated array of output matrices
!> \param matrix_ptrs preallocated array of output pointers; pointees are allocated here
!> \note Supply exactly one output representation. Flattened rows use (i-1)*ncols+j.
! **************************************************************************************************
   SUBROUTINE cp_array_to_dbcsr_slices(array, nrows, ncols, para_env_sub, &
                                       first_row, last_row, nslices, row_distribution, template, &
                                       matrices, matrix_ptrs)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: array
      INTEGER, INTENT(IN)                                :: nrows, ncols
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      INTEGER, INTENT(IN)                                :: first_row, last_row, nslices
      TYPE(group_dist_d1_type), INTENT(IN)               :: row_distribution
      TYPE(dbcsr_type), INTENT(INOUT)                    :: template
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: matrices
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: matrix_ptrs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_array_to_dbcsr_slices'

      INTEGER :: dimen_ia, dummy_proc, handle, i_global, i_local, iaia, iiB, iii, itmp(2), &
         j_global, j_local, jjB, jjj, kkB, mypcol, myprow, ncol_local, npcol, nprow, nrow_local, &
         number_of_rec, number_of_send, proc_receive, proc_send, proc_shift, rec_counter, &
         rec_iaia_end, rec_iaia_size, rec_iaia_start, rec_pcol, rec_prow, ref_send_pcol, &
         ref_send_prow, send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii_vet, map_rec_size, map_send_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: grid_2_mepos, grid_ref_2_send_pos, &
                                                            indices_map_my, mepos_2_grid
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: part_ia
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_ia
      TYPE(dbcsr_type), POINTER                          :: output_matrix
      TYPE(index_map), ALLOCATABLE, DIMENSION(:)         :: indices_rec
      TYPE(matrix_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: req_send

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(matrices) .NEQV. PRESENT(matrix_ptrs))
      dimen_ia = nrows*ncols

      ! create sub blacs env
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_sub)

      ! create the fm_ia buffer matrix
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrows, &
                               ncol_global=ncols, para_env=para_env_sub)
      CALL cp_fm_create(fm_ia, fm_struct, name="fm_ia")
      ! release structure
      CALL cp_fm_struct_release(fm_struct)
      ! release blacs_env
      CALL cp_blacs_env_release(blacs_env)

      ! get array information
      CALL cp_fm_get_info(matrix=fm_ia, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      myprow = fm_ia%matrix_struct%context%mepos(1)
      mypcol = fm_ia%matrix_struct%context%mepos(2)
      nprow = fm_ia%matrix_struct%context%num_pe(1)
      npcol = fm_ia%matrix_struct%context%num_pe(2)

      ! 0) create array containing the processes position and supporting infos
      ALLOCATE (grid_2_mepos(0:nprow - 1, 0:npcol - 1))
      grid_2_mepos = 0
      ALLOCATE (mepos_2_grid(2, 0:para_env_sub%num_pe - 1))
      ! fill the info array
      grid_2_mepos(myprow, mypcol) = para_env_sub%mepos
      ! sum infos
      CALL para_env_sub%sum(grid_2_mepos)
      CALL para_env_sub%allgather([myprow, mypcol], mepos_2_grid)

      ! loop over local index range and define the sending map
      ALLOCATE (map_send_size(0:para_env_sub%num_pe - 1))
      map_send_size = 0
      dummy_proc = 0
      DO iaia = first_row, last_row
         i_global = (iaia - 1)/ncols + 1
         j_global = MOD(iaia - 1, ncols) + 1
         send_prow = fm_ia%matrix_struct%g2p_row(i_global)
         send_pcol = fm_ia%matrix_struct%g2p_col(j_global)
         proc_send = grid_2_mepos(send_prow, send_pcol)
         map_send_size(proc_send) = map_send_size(proc_send) + 1
      END DO

      ! loop over local data of fm_ia and define the receiving map
      ALLOCATE (map_rec_size(0:para_env_sub%num_pe - 1))
      map_rec_size = 0
      part_ia = REAL(dimen_ia, KIND=dp)/REAL(para_env_sub%num_pe, KIND=dp)

      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            iaia = (i_global - 1)*ncols + j_global
            proc_receive = INT(REAL(iaia - 1, KIND=dp)/part_ia)
            proc_receive = MAX(0, proc_receive)
            proc_receive = MIN(proc_receive, para_env_sub%num_pe - 1)
            DO
               itmp = get_limit(dimen_ia, para_env_sub%num_pe, proc_receive)
               IF (iaia >= itmp(1) .AND. iaia <= itmp(2)) EXIT
               IF (iaia < itmp(1)) proc_receive = proc_receive - 1
               IF (iaia > itmp(2)) proc_receive = proc_receive + 1
            END DO
            map_rec_size(proc_receive) = map_rec_size(proc_receive) + 1
         END DO
      END DO

      ! allocate the buffer for sending data
      number_of_send = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be sent
      ALLOCATE (buffer_send(number_of_send))
      ! and the map from the grid of processess to the message position
      ALLOCATE (grid_ref_2_send_pos(0:nprow - 1, 0:npcol - 1))
      grid_ref_2_send_pos = 0
      ! finally allocate each message
      send_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
         size_send_buffer = map_send_size(proc_send)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            ! allocate the sending buffer (msg)
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%proc = proc_send
            ! get the pointer to prow, pcol of the process that has
            ! to receive this message
            ref_send_prow = mepos_2_grid(1, proc_send)
            ref_send_pcol = mepos_2_grid(2, proc_send)
            ! save the rank of the process that has to receive this message
            grid_ref_2_send_pos(ref_send_prow, ref_send_pcol) = send_counter
         END IF
      END DO

      ! allocate the buffer for receiving data
      number_of_rec = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be received
      ! and relative indices
      ALLOCATE (buffer_rec(number_of_rec))
      ALLOCATE (indices_rec(number_of_rec))
      ! finally allocate each message and fill the array of indices
      rec_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%proc = proc_receive
            ! create the indices array
            ALLOCATE (indices_rec(rec_counter)%map(2, size_rec_buffer))
            indices_rec(rec_counter)%map = 0
            CALL get_group_dist(row_distribution, proc_receive, rec_iaia_start, rec_iaia_end, rec_iaia_size)
            iii = 0
            DO iaia = rec_iaia_start, rec_iaia_end
               i_global = (iaia - 1)/ncols + 1
               j_global = MOD(iaia - 1, ncols) + 1
               rec_prow = fm_ia%matrix_struct%g2p_row(i_global)
               rec_pcol = fm_ia%matrix_struct%g2p_col(j_global)
               IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
               iii = iii + 1
               i_local = fm_ia%matrix_struct%g2l_row(i_global)
               j_local = fm_ia%matrix_struct%g2l_col(j_global)
               indices_rec(rec_counter)%map(1, iii) = i_local
               indices_rec(rec_counter)%map(2, iii) = j_local
            END DO
         END IF
      END DO
      ! and create the index map for my local data
      IF (map_rec_size(para_env_sub%mepos) > 0) THEN
         size_rec_buffer = map_rec_size(para_env_sub%mepos)
         ALLOCATE (indices_map_my(2, size_rec_buffer))
         indices_map_my = 0
         iii = 0
         DO iaia = first_row, last_row
            i_global = (iaia - 1)/ncols + 1
            j_global = MOD(iaia - 1, ncols) + 1
            rec_prow = fm_ia%matrix_struct%g2p_row(i_global)
            rec_pcol = fm_ia%matrix_struct%g2p_col(j_global)
            IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
            iii = iii + 1
            i_local = fm_ia%matrix_struct%g2l_row(i_global)
            j_local = fm_ia%matrix_struct%g2l_col(j_global)
            indices_map_my(1, iii) = i_local
            indices_map_my(2, iii) = j_local
         END DO
      END IF

      ! auxiliary vector of indices for the send buffer
      ALLOCATE (iii_vet(number_of_send))
      ! vector for the send requests
      ALLOCATE (req_send(number_of_send))
      ! loop over auxiliary basis function and redistribute into a fm
      ! and then compy the fm into a dbcsr matrix
      DO kkB = 1, nslices
         ! zero the matries of the buffers and post the messages to be received
         CALL cp_fm_set_all(matrix=fm_ia, alpha=0.0_dp)
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter + 1
               buffer_rec(rec_counter)%msg = 0.0_dp
               CALL para_env_sub%irecv(buffer_rec(rec_counter)%msg, proc_receive, &
                                       buffer_rec(rec_counter)%msg_req)
            END IF
         END DO
         ! fill the sending buffer and send the messages
         DO send_counter = 1, number_of_send
            buffer_send(send_counter)%msg = 0.0_dp
         END DO
         iii_vet = 0
         jjj = 0
         DO iaia = first_row, last_row
            i_global = (iaia - 1)/ncols + 1
            j_global = MOD(iaia - 1, ncols) + 1
            send_prow = fm_ia%matrix_struct%g2p_row(i_global)
            send_pcol = fm_ia%matrix_struct%g2p_col(j_global)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            ! we don't need to send to ourselves
            IF (grid_2_mepos(send_prow, send_pcol) == para_env_sub%mepos) THEN
               ! filling fm_ia with local data
               jjj = jjj + 1
               i_local = indices_map_my(1, jjj)
               j_local = indices_map_my(2, jjj)
               fm_ia%local_data(i_local, j_local) = array(iaia - first_row + 1, kkB)
            ELSE
               send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
               iii_vet(send_counter) = iii_vet(send_counter) + 1
               iii = iii_vet(send_counter)
               buffer_send(send_counter)%msg(iii) = array(iaia - first_row + 1, kkB)
            END IF
         END DO
         req_send = mp_request_null
         send_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
            IF (map_send_size(proc_send) > 0) THEN
               send_counter = send_counter + 1
               CALL para_env_sub%isend(buffer_send(send_counter)%msg, proc_send, &
                                       buffer_send(send_counter)%msg_req)
               req_send(send_counter) = buffer_send(send_counter)%msg_req
            END IF
         END DO

         ! receive the massages and fill the fm_ia
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
            size_rec_buffer = map_rec_size(proc_receive)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter + 1
               ! wait for the message
               CALL buffer_rec(rec_counter)%msg_req%wait()
               DO iii = 1, size_rec_buffer
                  i_local = indices_rec(rec_counter)%map(1, iii)
                  j_local = indices_rec(rec_counter)%map(2, iii)
                  fm_ia%local_data(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
               END DO
            END IF
         END DO

         ! wait all
         CALL mp_waitall(req_send(:))

         ! now create the DBCSR matrix and copy fm_ia into it
         IF (PRESENT(matrices)) THEN
            output_matrix => matrices(kkB)
         ELSE
            ALLOCATE (matrix_ptrs(kkB)%matrix)
            output_matrix => matrix_ptrs(kkB)%matrix
         END IF
         CALL cp_dbcsr_m_by_n_from_template(output_matrix, template=template, &
                                            m=nrows, n=ncols, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_ia, output_matrix, keep_sparsity=.FALSE.)

      END DO

      DEALLOCATE (array)

      ! release buffer matrix
      CALL cp_fm_release(fm_ia)

      CALL timestop(handle)

   END SUBROUTINE cp_array_to_dbcsr_slices

END MODULE cp_dbcsr_redistribution

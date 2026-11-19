!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Binary file I/O for DBCSR/DBM matrices built on cp_dbcsr_api (DBCSRv_1.0 format)
! **************************************************************************************************
MODULE cp_dbcsr_io
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_hold, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_info, dbcsr_get_stored_coordinates, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_readonly_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_put_block, dbcsr_release, dbcsr_type, &
        dbcsr_valid_index, dbcsr_work_create
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8,&
                                              real_4,&
                                              real_4_size,&
                                              real_8_size
   USE message_passing,                 ONLY: &
        file_amode_create, file_amode_rdonly, file_amode_wronly, file_offset, mp_comm_type, &
        mp_file_delete, mp_file_type, mpi_character_size, mpi_integer_size
   USE util,                            ONLY: sort
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_io'

   INTEGER, PARAMETER, PRIVATE :: version_len = 10
   INTEGER, PARAMETER, PRIVATE :: name_v_1_0_len = 80
   INTEGER, PARAMETER, PRIVATE :: char_count = version_len + name_v_1_0_len + 1
   CHARACTER(LEN=version_len), PARAMETER, PRIVATE :: version_v_1_0 = "DBCSRv_1.0"

   INTEGER, PARAMETER, PRIVATE :: dbcsr_type_real_4_id = 1
   INTEGER, PARAMETER, PRIVATE :: dbcsr_type_real_8_id = 3

   PUBLIC :: dbcsr_binary_write
   PUBLIC :: dbcsr_binary_read

CONTAINS

! **************************************************************************************************
!> \brief Write a DBCSR matrix to a binary file in DBCSRv_1.0 format
!> \param matrix ...
!> \param filepath ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_write(matrix, filepath)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: filepath

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_binary_write'
      INTEGER(KIND=file_offset), PARAMETER               :: bof = 0_file_offset

      CHARACTER                                          :: matrix_type
      CHARACTER(LEN=default_string_length)               :: matrix_name
      CHARACTER(LEN=name_v_1_0_len)                      :: matrix_name_v_1_0
      INTEGER :: bcol, brow, cnt, csize, cur_off, data_area_size, handle, k, linfo_size, &
         mp_handle, my_off, mynode, nblkcols_total, nblkrows_total, nblks, numnodes, nze, &
         orig_slot, r, rsize, s_end, s_start, slot
      INTEGER(KIND=file_offset) :: bdata_my_offset, bdata_start, prev_da, prev_linfo, &
         subh1_chars_bytes, subh2_my_offset, subh2_start, subh3_my_offset, subh3_start, tot_linfo
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_p, col_i, col_unsorted, da_sizes, &
                                                            hdr1_ints, linfo_sizes, linfo_vec, &
                                                            perm, raw_off, row_cnt, row_p, row_pos
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: col_blk_size, row_blk_size
      INTEGER, DIMENSION(2)                              :: subh2_pair
      LOGICAL                                            :: tr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: r_dp, raw_data
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_ptr
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: mp_group
      TYPE(mp_file_type)                                 :: thefile

      CALL timeset(routineN, handle)

      IF (.NOT. dbcsr_valid_index(matrix)) THEN
         CPABORT("Invalid DBCSR matrix in dbcsr_binary_write")
      END IF

      NULLIFY (row_blk_size, col_blk_size, blk_ptr)
      CALL dbcsr_get_info(matrix, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          matrix_type=matrix_type, &
                          name=matrix_name, &
                          distribution=distribution)

      CALL dbcsr_distribution_get(distribution, group=mp_handle, mynode=mynode, numnodes=numnodes)
      CALL mp_group%set_handle(mp_handle)

      ! Pass 1: count local blocks per row and total local element count
      ALLOCATE (row_cnt(nblkrows_total))
      row_cnt(:) = 0
      nblks = 0
      data_area_size = 0

!$OMP PARALLEL DEFAULT(NONE) SHARED(matrix, row_cnt) &
!$OMP PRIVATE(iter, brow, bcol, blk_ptr, rsize, csize, tr) &
!$OMP REDUCTION(+:nblks, data_area_size)
      NULLIFY (blk_ptr)
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, brow, bcol, blk_ptr, &
                                        row_size=rsize, col_size=csize, transposed=tr)
         nblks = nblks + 1
!$OMP ATOMIC UPDATE
         row_cnt(brow) = row_cnt(brow) + 1
         data_area_size = data_area_size + rsize*csize
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      ALLOCATE (row_p(nblkrows_total + 1))
      row_p(1) = 0
      DO r = 1, nblkrows_total
         row_p(r + 1) = row_p(r) + row_cnt(r)
      END DO

      ! Pass 2: collect block columns and untranposed block data, then sort columns per row
      ALLOCATE (col_unsorted(nblks), raw_off(nblks), col_i(nblks), blk_p(nblks))
      ALLOCATE (raw_data(MAX(1, data_area_size)), r_dp(MAX(1, data_area_size)))
      raw_data(:) = 0.0_dp
      r_dp(:) = 0.0_dp
      ALLOCATE (row_pos(nblkrows_total))
      IF (nblkrows_total > 0) row_pos(:) = row_p(1:nblkrows_total)

      cur_off = 1
!$OMP PARALLEL DEFAULT(NONE) SHARED(matrix, row_pos, col_unsorted, raw_off, raw_data, cur_off) &
!$OMP PRIVATE(iter, brow, bcol, blk_ptr, rsize, csize, tr, nze, slot, my_off)
      NULLIFY (blk_ptr)
      CALL dbcsr_iterator_readonly_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, brow, bcol, blk_ptr, &
                                        row_size=rsize, col_size=csize, transposed=tr)
         nze = rsize*csize
!$OMP ATOMIC CAPTURE
         row_pos(brow) = row_pos(brow) + 1
         slot = row_pos(brow)
!$OMP END ATOMIC
!$OMP ATOMIC CAPTURE
         my_off = cur_off
         cur_off = cur_off + nze
!$OMP END ATOMIC
         col_unsorted(slot) = bcol
         raw_off(slot) = my_off
         IF (nze > 0) THEN
            IF (tr) THEN
               raw_data(my_off:my_off + nze - 1) = RESHAPE(TRANSPOSE(blk_ptr), [nze])
            ELSE
               raw_data(my_off:my_off + nze - 1) = RESHAPE(blk_ptr, [nze])
            END IF
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      IF (nblkrows_total > 0) THEN
         ALLOCATE (perm(MAX(1, MAXVAL(row_cnt))))
         cur_off = 1
         DO r = 1, nblkrows_total
            cnt = row_cnt(r)
            IF (cnt > 0) THEN
               s_start = row_p(r) + 1
               s_end = row_p(r + 1)
               col_i(s_start:s_end) = col_unsorted(s_start:s_end)
               IF (cnt > 1) THEN
                  CALL sort(col_i(s_start:s_end), cnt, perm(1:cnt))
               ELSE
                  perm(1) = 1
               END IF
               rsize = row_blk_size(r)
               DO k = 1, cnt
                  orig_slot = s_start + perm(k) - 1
                  bcol = col_i(s_start + k - 1)
                  nze = rsize*col_blk_size(bcol)
                  blk_p(s_start + k - 1) = cur_off
                  IF (nze > 0) THEN
                     r_dp(cur_off:cur_off + nze - 1) = raw_data(raw_off(orig_slot):raw_off(orig_slot) + nze - 1)
                  END IF
                  cur_off = cur_off + nze
               END DO
            END IF
         END DO
         DEALLOCATE (perm)
      END IF
      DEALLOCATE (row_cnt, row_pos, col_unsorted, raw_off, raw_data)

      linfo_size = 1 + nblkrows_total + 2*nblks
      ALLOCATE (linfo_vec(linfo_size))
      linfo_vec(1:nblkrows_total + 1) = row_p(1:nblkrows_total + 1)
      IF (nblks > 0) THEN
         linfo_vec(nblkrows_total + 2:nblkrows_total + 1 + nblks) = col_i(1:nblks)
         linfo_vec(nblkrows_total + 2 + nblks:linfo_size) = blk_p(1:nblks)
      END IF
      DEALLOCATE (row_p, col_i, blk_p)

      ALLOCATE (linfo_sizes(numnodes), da_sizes(numnodes))
      CALL mp_group%allgather(linfo_size, linfo_sizes)
      CALL mp_group%allgather(data_area_size, da_sizes)

      subh1_chars_bytes = INT(char_count, KIND=file_offset)*INT(mpi_character_size, KIND=file_offset)
      subh2_start = subh1_chars_bytes + &
                    INT(4 + nblkrows_total + nblkcols_total, KIND=file_offset)* &
                    INT(mpi_integer_size, KIND=file_offset)
      subh2_my_offset = subh2_start + INT(2*mynode, KIND=file_offset)*INT(mpi_integer_size, KIND=file_offset)

      subh3_start = subh2_start + INT(2*numnodes, KIND=file_offset)*INT(mpi_integer_size, KIND=file_offset)
      prev_linfo = SUM(INT(linfo_sizes(1:mynode), KIND=file_offset))
      tot_linfo = SUM(INT(linfo_sizes(1:numnodes), KIND=file_offset))
      subh3_my_offset = subh3_start + prev_linfo*INT(mpi_integer_size, KIND=file_offset)

      bdata_start = subh3_start + tot_linfo*INT(mpi_integer_size, KIND=file_offset)
      prev_da = SUM(INT(da_sizes(1:mynode), KIND=file_offset))
      bdata_my_offset = bdata_start + prev_da*INT(real_8_size, KIND=file_offset)

      DEALLOCATE (linfo_sizes, da_sizes)

      IF (mynode == 0) CALL mp_file_delete(filepath)
      CALL mp_group%sync()

      CALL thefile%open(mp_group, filepath, file_amode_create + file_amode_wronly)

      IF (mynode == 0) THEN
         matrix_name_v_1_0 = matrix_name
         CALL thefile%write_at(bof, version_v_1_0)
         CALL thefile%write_at(bof + INT(version_len*mpi_character_size, KIND=file_offset), matrix_name_v_1_0)
         CALL thefile%write_at(bof + INT((version_len + name_v_1_0_len)*mpi_character_size, KIND=file_offset), &
                               matrix_type)

         ALLOCATE (hdr1_ints(4 + nblkrows_total + nblkcols_total))
         hdr1_ints(1:4) = [numnodes, dbcsr_type_real_8_id, nblkrows_total, nblkcols_total]
         IF (nblkrows_total > 0) THEN
            hdr1_ints(5:4 + nblkrows_total) = row_blk_size(1:nblkrows_total)
         END IF
         IF (nblkcols_total > 0) THEN
            hdr1_ints(5 + nblkrows_total:4 + nblkrows_total + nblkcols_total) = col_blk_size(1:nblkcols_total)
         END IF
         CALL thefile%write_at(bof + subh1_chars_bytes, hdr1_ints)
         DEALLOCATE (hdr1_ints)
      END IF

      subh2_pair = [nblks, data_area_size]
      CALL thefile%write_at_all(subh2_my_offset, subh2_pair)
      CALL thefile%write_at_all(subh3_my_offset, linfo_vec)
      CALL thefile%write_at_all(bdata_my_offset, r_dp, msglen=data_area_size)

      CALL thefile%close()

      DEALLOCATE (linfo_vec, r_dp)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_binary_write

! **************************************************************************************************
!> \brief Read a DBCSR matrix from a binary DBCSRv_1.0 file into the given distribution
!> \param filepath ...
!> \param distribution ...
!> \param matrix_new ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_read(filepath, distribution, matrix_new)
      CHARACTER(len=*), INTENT(IN)                       :: filepath
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: distribution
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_new

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbcsr_binary_read'
      INTEGER(KIND=file_offset), PARAMETER               :: bof = 0_file_offset
      INTEGER, PARAMETER                                 :: metalen = 2

      CHARACTER                                          :: matrix_type
      CHARACTER(LEN=default_string_length)               :: matrix_name
      CHARACTER(LEN=name_v_1_0_len)                      :: matrix_name_v_1_0
      CHARACTER(LEN=version_len)                         :: version
      INTEGER :: bcol, blk, blk_size, blockdata_length, brow, csize, cur_blks, darea_size, &
         data_type, dst_p, end_index, globalinfo_size, group_list_size, handle, i, job_count, k, &
         linfo_length, localinfo_length, m_pos, mp_handle, nblkcols_total, nblkrows_total, nblks, &
         order, proc_idx, rsize, share_size, size_of_pgrid, start_index, total_recv_blks, &
         total_recv_vals, total_send_blks, total_send_vals, type_size, v_pos, val_pos, worker_id
      INTEGER(KIND=file_offset) :: bdata_start, blockdata_offset, localinfo_offset, offset, &
         subh1_chars_bytes, subh2_start, subh3_length, subh3_start, sum_nblks
      INTEGER(KIND=file_offset), ALLOCATABLE, &
         DIMENSION(:)                                    :: bdata_offsets, subh3_offsets
      INTEGER(KIND=int_8)                                :: actual_blk, blkp_abs, end_ind, start_ind
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: extra_darea_size, extra_nblks
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bdata_lens, blkp, coli, ginfo_vec, linfo_lens, &
         linfo_vec, recv_meta, rmeta_cnt, rmeta_displ, rowp, rval_cnt, rval_displ, send_meta, &
         smeta_cnt, smeta_displ, smeta_pos, sval_cnt, sval_displ, sval_pos
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: recv_counts, send_counts, val_data
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: col_blk_size, row_blk_size
      INTEGER, DIMENSION(2)                              :: subh2_pair
      INTEGER, DIMENSION(4)                              :: hdr1_vals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rdp, recv_data, send_data
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk2d
      REAL(KIND=real_4), ALLOCATABLE, DIMENSION(:)       :: rsp
      TYPE(dbcsr_distribution_type)                      :: dist_local
      TYPE(mp_comm_type)                                 :: mp_group
      TYPE(mp_file_type)                                 :: thefile

      CALL timeset(routineN, handle)

      CALL dbcsr_distribution_get(distribution, group=mp_handle, &
                                  mynode=worker_id, numnodes=group_list_size)
      CALL mp_group%set_handle(mp_handle)

      CALL thefile%open(mp_group, filepath, file_amode_rdonly)

      CALL thefile%read_at_all(bof, version)
      IF (version /= version_v_1_0) THEN
         CPABORT("Unsupported DBCSR binary file version in "//TRIM(filepath))
      END IF

      CALL thefile%read_at_all(bof + INT(version_len*mpi_character_size, KIND=file_offset), &
                               matrix_name_v_1_0)
      matrix_name = matrix_name_v_1_0

      CALL thefile%read_at_all(bof + INT((version_len + name_v_1_0_len)*mpi_character_size, KIND=file_offset), &
                               matrix_type)

      subh1_chars_bytes = INT(char_count, KIND=file_offset)*INT(mpi_character_size, KIND=file_offset)
      CALL thefile%read_at_all(bof + subh1_chars_bytes, hdr1_vals)
      size_of_pgrid = hdr1_vals(1)
      data_type = hdr1_vals(2)
      nblkrows_total = hdr1_vals(3)
      nblkcols_total = hdr1_vals(4)

      IF (data_type /= dbcsr_type_real_8_id .AND. data_type /= dbcsr_type_real_4_id) THEN
         CPABORT("Unsupported data_type in DBCSR binary file: "//TRIM(filepath))
      END IF

      globalinfo_size = nblkrows_total + nblkcols_total
      ALLOCATE (ginfo_vec(MAX(1, globalinfo_size)))
      IF (globalinfo_size > 0) THEN
         CALL thefile%read_at_all(bof + subh1_chars_bytes + INT(4*mpi_integer_size, KIND=file_offset), &
                                  ginfo_vec(1:globalinfo_size), msglen=globalinfo_size)
      END IF

      ALLOCATE (row_blk_size(nblkrows_total), col_blk_size(nblkcols_total))
      IF (nblkrows_total > 0) row_blk_size(:) = ginfo_vec(1:nblkrows_total)
      IF (nblkcols_total > 0) col_blk_size(:) = ginfo_vec(nblkrows_total + 1:globalinfo_size)
      DEALLOCATE (ginfo_vec)

      subh2_start = subh1_chars_bytes + INT(4 + globalinfo_size, KIND=file_offset)*INT(mpi_integer_size, KIND=file_offset)
      subh3_start = subh2_start + INT(2*size_of_pgrid, KIND=file_offset)*INT(mpi_integer_size, KIND=file_offset)

      share_size = CEILING(REAL(size_of_pgrid, KIND=dp)/REAL(group_list_size, KIND=dp))
      ALLOCATE (val_data(3, share_size))
      val_data(:, :) = 0

      DO i = 1, share_size
         proc_idx = worker_id + (i - 1)*group_list_size
         IF (proc_idx < size_of_pgrid) THEN
            offset = subh2_start + INT(2*proc_idx, KIND=file_offset)*INT(mpi_integer_size, KIND=file_offset)
            CALL thefile%read_at_all(offset, subh2_pair, msglen=2)
            val_data(1, i) = subh2_pair(1)
            val_data(2, i) = subh2_pair(2)
            val_data(3, i) = proc_idx + 1
         ELSE
            CALL thefile%read_at_all(bof, subh2_pair, msglen=0)
         END IF
      END DO

      nblks = SUM(val_data(1, :))
      darea_size = SUM(val_data(2, :))
      job_count = COUNT(val_data(3, :) /= 0)

      ALLOCATE (linfo_lens(size_of_pgrid), subh3_offsets(size_of_pgrid))
      linfo_lens(:) = 0
      DO i = 1, share_size
         order = val_data(3, i)
         IF (order > 0) THEN
            linfo_lens(order) = 1 + nblkrows_total + 2*val_data(1, i)
         END IF
      END DO
      CALL mp_group%sum(linfo_lens)

      offset = subh3_start
      DO i = 1, size_of_pgrid
         subh3_offsets(i) = offset
         offset = offset + INT(linfo_lens(i), KIND=file_offset)*INT(mpi_integer_size, KIND=file_offset)
      END DO

      sum_nblks = INT(nblks, KIND=file_offset)
      CALL mp_group%sum(sum_nblks)
      subh3_length = INT(size_of_pgrid, KIND=file_offset)*INT(1 + nblkrows_total, KIND=file_offset) + 2_file_offset*sum_nblks

      linfo_length = nblkrows_total + 1 + 2*MAX(0, MAXVAL(val_data(1, :)))
      ALLOCATE (linfo_vec(MAX(1, linfo_length)))
      ALLOCATE (rowp((nblkrows_total + 1)*job_count))
      ALLOCATE (coli(nblks), blkp(nblks))

      DO i = 1, share_size
         order = val_data(3, i)
         cur_blks = val_data(1, i)
         IF (order == 0) THEN
            CALL thefile%read_at_all(bof, linfo_vec, msglen=0)
         ELSE
            localinfo_offset = subh3_offsets(order)
            localinfo_length = linfo_lens(order)
            CALL thefile%read_at_all(localinfo_offset, linfo_vec(1:localinfo_length), msglen=localinfo_length)
            rowp((i - 1)*(nblkrows_total + 1) + 1:i*(nblkrows_total + 1)) = linfo_vec(1:nblkrows_total + 1)
            IF (cur_blks > 0) THEN
               start_index = SUM(val_data(1, 1:i - 1)) + 1
               end_index = start_index + cur_blks - 1
               coli(start_index:end_index) = linfo_vec(nblkrows_total + 2:nblkrows_total + 1 + cur_blks)
               blkp(start_index:end_index) = linfo_vec(nblkrows_total + 2 + cur_blks:nblkrows_total + 1 + 2*cur_blks)
            END IF
         END IF
      END DO
      DEALLOCATE (linfo_lens, subh3_offsets, linfo_vec)

      ALLOCATE (bdata_lens(size_of_pgrid), bdata_offsets(size_of_pgrid))
      bdata_lens(:) = 0
      DO i = 1, share_size
         order = val_data(3, i)
         IF (order > 0) bdata_lens(order) = val_data(2, i)
      END DO
      CALL mp_group%sum(bdata_lens)

      IF (data_type == dbcsr_type_real_8_id) THEN
         type_size = real_8_size
      ELSE
         type_size = real_4_size
      END IF

      bdata_start = subh3_start + subh3_length*INT(mpi_integer_size, KIND=file_offset)
      offset = bdata_start
      DO i = 1, size_of_pgrid
         bdata_offsets(i) = offset
         offset = offset + INT(bdata_lens(i), KIND=file_offset)*INT(type_size, KIND=file_offset)
      END DO

      ALLOCATE (rdp(MAX(1, darea_size)))
      rdp(:) = 0.0_dp
      IF (data_type == dbcsr_type_real_8_id) THEN
         DO i = 1, share_size
            order = val_data(3, i)
            IF (order == 0) THEN
               CALL thefile%read_at_all(bof, rdp, msglen=0)
            ELSE
               blockdata_length = bdata_lens(order)
               IF (blockdata_length > 0) THEN
                  start_index = SUM(val_data(2, 1:i - 1)) + 1
                  end_index = start_index + blockdata_length - 1
                  blockdata_offset = bdata_offsets(order)
                  CALL thefile%read_at_all(blockdata_offset, rdp(start_index:end_index), msglen=blockdata_length)
               ELSE
                  CALL thefile%read_at_all(bof, rdp, msglen=0)
               END IF
            END IF
         END DO
      ELSE
         ALLOCATE (rsp(MAX(1, darea_size)))
         rsp(:) = 0.0_real_4
         DO i = 1, share_size
            order = val_data(3, i)
            IF (order == 0) THEN
               CALL thefile%read_at_all(bof, rsp, msglen=0)
            ELSE
               blockdata_length = bdata_lens(order)
               IF (blockdata_length > 0) THEN
                  start_index = SUM(val_data(2, 1:i - 1)) + 1
                  end_index = start_index + blockdata_length - 1
                  blockdata_offset = bdata_offsets(order)
                  CALL thefile%read_at_all(blockdata_offset, rsp(start_index:end_index), msglen=blockdata_length)
               ELSE
                  CALL thefile%read_at_all(bof, rsp, msglen=0)
               END IF
            END IF
         END DO
         IF (darea_size > 0) rdp(1:darea_size) = REAL(rsp(1:darea_size), KIND=dp)
         DEALLOCATE (rsp)
      END IF

      CALL thefile%close()
      DEALLOCATE (bdata_lens, bdata_offsets)

      ! Hold distribution before releasing any pre-existing matrix_new that might own it
      dist_local = distribution
      CALL dbcsr_distribution_hold(dist_local)
      IF (dbcsr_valid_index(matrix_new)) THEN
         CALL dbcsr_release(matrix_new)
      END IF
      CALL dbcsr_create(matrix_new, name=matrix_name, dist=dist_local, matrix_type=matrix_type, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_distribution_release(dist_local)

      ! Redistribute blocks read by this rank to their target owners in matrix_new
      ALLOCATE (extra_nblks(job_count), extra_darea_size(job_count))
      DO i = 1, job_count
         extra_nblks(i) = SUM(INT(val_data(1, 1:i - 1), KIND=int_8))
         extra_darea_size(i) = SUM(INT(val_data(2, 1:i - 1), KIND=int_8))
      END DO

      ALLOCATE (send_counts(2, group_list_size), recv_counts(2, group_list_size))
      send_counts(:, :) = 0
      recv_counts(:, :) = 0

      DO i = 1, job_count
         DO brow = 1, nblkrows_total
            rsize = row_blk_size(brow)
            DO blk = rowp((i - 1)*(nblkrows_total + 1) + brow) + 1, &
               rowp((i - 1)*(nblkrows_total + 1) + brow + 1)
               actual_blk = INT(blk, KIND=int_8) + extra_nblks(i)
               IF (blkp(actual_blk) == 0) CYCLE
               bcol = coli(actual_blk)
               csize = col_blk_size(bcol)
               CALL dbcsr_get_stored_coordinates(matrix_new, brow, bcol, dst_p)
               send_counts(1, dst_p + 1) = send_counts(1, dst_p + 1) + 1
               send_counts(2, dst_p + 1) = send_counts(2, dst_p + 1) + rsize*csize
            END DO
         END DO
      END DO

      CALL mp_group%alltoall(send_counts, recv_counts, 2)

      ALLOCATE (smeta_cnt(group_list_size), smeta_displ(group_list_size), smeta_pos(group_list_size))
      ALLOCATE (rmeta_cnt(group_list_size), rmeta_displ(group_list_size))
      ALLOCATE (sval_cnt(group_list_size), sval_displ(group_list_size), sval_pos(group_list_size))
      ALLOCATE (rval_cnt(group_list_size), rval_displ(group_list_size))

      DO i = 1, group_list_size
         smeta_cnt(i) = metalen*send_counts(1, i)
         rmeta_cnt(i) = metalen*recv_counts(1, i)
         sval_cnt(i) = send_counts(2, i)
         rval_cnt(i) = recv_counts(2, i)
      END DO

      smeta_displ(1) = 0
      rmeta_displ(1) = 0
      sval_displ(1) = 0
      rval_displ(1) = 0
      DO i = 2, group_list_size
         smeta_displ(i) = smeta_displ(i - 1) + smeta_cnt(i - 1)
         rmeta_displ(i) = rmeta_displ(i - 1) + rmeta_cnt(i - 1)
         sval_displ(i) = sval_displ(i - 1) + sval_cnt(i - 1)
         rval_displ(i) = rval_displ(i - 1) + rval_cnt(i - 1)
      END DO
      smeta_pos(:) = smeta_displ(:)
      sval_pos(:) = sval_displ(:)

      total_send_blks = SUM(send_counts(1, :))
      total_send_vals = SUM(send_counts(2, :))
      total_recv_blks = SUM(recv_counts(1, :))
      total_recv_vals = SUM(recv_counts(2, :))

      ALLOCATE (send_meta(MAX(1, metalen*total_send_blks)), recv_meta(MAX(1, metalen*total_recv_blks)))
      ALLOCATE (send_data(MAX(1, total_send_vals)), recv_data(MAX(1, total_recv_vals)))

      DO i = 1, job_count
         DO brow = 1, nblkrows_total
            rsize = row_blk_size(brow)
            DO blk = rowp((i - 1)*(nblkrows_total + 1) + brow) + 1, &
               rowp((i - 1)*(nblkrows_total + 1) + brow + 1)
               actual_blk = INT(blk, KIND=int_8) + extra_nblks(i)
               IF (blkp(actual_blk) == 0) CYCLE
               bcol = coli(actual_blk)
               csize = col_blk_size(bcol)
               blk_size = rsize*csize
               blkp_abs = INT(ABS(blkp(actual_blk)), KIND=int_8)
               start_ind = blkp_abs + extra_darea_size(i)
               end_ind = start_ind + INT(blk_size, KIND=int_8) - 1_int_8

               CALL dbcsr_get_stored_coordinates(matrix_new, brow, bcol, dst_p)
               m_pos = smeta_pos(dst_p + 1)
               send_meta(m_pos + 1) = brow
               send_meta(m_pos + 2) = bcol
               smeta_pos(dst_p + 1) = m_pos + metalen

               IF (blk_size > 0) THEN
                  v_pos = sval_pos(dst_p + 1)
                  IF (blkp(actual_blk) < 0) THEN
                     send_data(v_pos + 1:v_pos + blk_size) = &
                        RESHAPE(TRANSPOSE(RESHAPE(rdp(start_ind:end_ind), [csize, rsize])), [blk_size])
                  ELSE
                     send_data(v_pos + 1:v_pos + blk_size) = rdp(start_ind:end_ind)
                  END IF
                  sval_pos(dst_p + 1) = v_pos + blk_size
               END IF
            END DO
         END DO
      END DO

      DEALLOCATE (val_data, extra_nblks, extra_darea_size, rowp, coli, blkp, rdp)

      CALL mp_group%alltoall(send_meta, smeta_cnt, smeta_displ, recv_meta, rmeta_cnt, rmeta_displ)
      CALL mp_group%alltoall(send_data, sval_cnt, sval_displ, recv_data, rval_cnt, rval_displ)

      DEALLOCATE (send_counts, recv_counts, send_meta, send_data)
      DEALLOCATE (smeta_cnt, smeta_displ, smeta_pos, rmeta_cnt, rmeta_displ)
      DEALLOCATE (sval_cnt, sval_displ, sval_pos, rval_cnt, rval_displ)

      CALL dbcsr_work_create(matrix_new, nblks_guess=total_recv_blks, sizedata_guess=total_recv_vals)
      val_pos = 1
      DO k = 1, total_recv_blks
         brow = recv_meta(metalen*(k - 1) + 1)
         bcol = recv_meta(metalen*(k - 1) + 2)
         rsize = row_blk_size(brow)
         csize = col_blk_size(bcol)
         blk_size = rsize*csize
         ALLOCATE (blk2d(rsize, csize))
         IF (blk_size > 0) THEN
            blk2d(:, :) = RESHAPE(recv_data(val_pos:val_pos + blk_size - 1), [rsize, csize])
         ELSE
            blk2d(:, :) = 0.0_dp
         END IF
         CALL dbcsr_put_block(matrix_new, brow, bcol, blk2d)
         DEALLOCATE (blk2d)
         val_pos = val_pos + blk_size
      END DO
      CALL dbcsr_finalize(matrix_new)

      DEALLOCATE (recv_meta, recv_data, row_blk_size, col_blk_size)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_binary_read

END MODULE cp_dbcsr_io

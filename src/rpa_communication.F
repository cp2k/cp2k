!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary rountines necessary to redistribute an fm_matrix from a
!>        given blacs_env to another
!> \par History
!>      12.2012 created [Mauro Del Ben]
! **************************************************************************************************
MODULE rpa_communication
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                              cp_dbcsr_init,&
                                              cp_dbcsr_p_type,&
                                              cp_dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_indxg2l,&
                                              cp_fm_indxg2p,&
                                              cp_fm_indxl2g,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_wait,&
                                              mp_waitall
   USE mp2_ri_grad_util,                ONLY: fm2array
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE index_map
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: map
   END TYPE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_communication'

   PUBLIC :: initialize_buffer, &
             fm_redistribute, &
             release_buffer, &
             gamma_fm_to_dbcsr

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_source ...
!> \param fm_mat_dest ...
!> \param RPA_proc_map ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param number_of_rec ...
!> \param number_of_send ...
!> \param map_send_size ...
!> \param map_rec_size ...
!> \param local_size_source ...
!> \param para_env_RPA ...
! **************************************************************************************************
   SUBROUTINE initialize_buffer(fm_mat_source, fm_mat_dest, RPA_proc_map, buffer_rec, buffer_send, &
                                number_of_rec, number_of_send, &
                                map_send_size, map_rec_size, local_size_source, para_env_RPA)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_source, fm_mat_dest
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: RPA_proc_map
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      INTEGER                                            :: number_of_rec, number_of_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: map_send_size, map_rec_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_size_source
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_buffer', &
         routineP = moduleN//':'//routineN

      INTEGER :: dummy_proc, handle, handle2, i, i_global, i_local, iiB, iii, j_global, j_local, &
         jjB, mypcol_d, mypcol_s, myprow_d, myprow_s, ncol_block_d, ncol_block_s, ncol_local_d, &
         ncol_local_s, npcol_d, npcol_s, nprow_d, nprow_s, nrow_block_d, nrow_block_s, &
         nrow_local_d, nrow_local_s, proc_receive, proc_send, proc_shift, rec_counter, &
         rec_local_col, rec_local_row, rec_pcol, rec_prow, ref_rec_pcol, ref_rec_prow, &
         send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: index_counter, proc2counter
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_d, col_indices_s, &
                                                            row_indices_d, row_indices_s

      CALL timeset(routineN, handle)

      ! create the RPA proc_map
      IF (.NOT. (ALLOCATED(RPA_proc_map))) THEN
         ALLOCATE (RPA_proc_map(-para_env_RPA%num_pe:2*para_env_RPA%num_pe-1))
         RPA_proc_map = 0
         DO i = 0, para_env_RPA%num_pe-1
            RPA_proc_map(i) = i
            RPA_proc_map(-i-1) = para_env_RPA%num_pe-i-1
            RPA_proc_map(para_env_RPA%num_pe+i) = i
         END DO
      END IF

      ! get info for the source
      CALL cp_fm_get_info(matrix=fm_mat_source, &
                          nrow_local=nrow_local_s, &
                          ncol_local=ncol_local_s, &
                          row_indices=row_indices_s, &
                          col_indices=col_indices_s, &
                          nrow_block=nrow_block_s, &
                          ncol_block=ncol_block_s)
      myprow_s = fm_mat_source%matrix_struct%context%mepos(1)
      mypcol_s = fm_mat_source%matrix_struct%context%mepos(2)
      nprow_s = fm_mat_source%matrix_struct%context%num_pe(1)
      npcol_s = fm_mat_source%matrix_struct%context%num_pe(2)

      ! get info for the dest
      CALL cp_fm_get_info(matrix=fm_mat_dest, &
                          nrow_local=nrow_local_d, &
                          ncol_local=ncol_local_d, &
                          row_indices=row_indices_d, &
                          col_indices=col_indices_d, &
                          nrow_block=nrow_block_d, &
                          ncol_block=ncol_block_d)
      myprow_d = fm_mat_dest%matrix_struct%context%mepos(1)
      mypcol_d = fm_mat_dest%matrix_struct%context%mepos(2)
      nprow_d = fm_mat_dest%matrix_struct%context%num_pe(1)
      npcol_d = fm_mat_dest%matrix_struct%context%num_pe(2)

      ! 0) create the map for the local sizes
      ALLOCATE (local_size_source(2, 0:para_env_RPA%num_pe-1))
      local_size_source = 0
      local_size_source(1, para_env_RPA%mepos) = nrow_local_s
      local_size_source(2, para_env_RPA%mepos) = ncol_local_s
      CALL mp_sum(local_size_source, para_env_RPA%group)

      ! 1) loop over my local data and define a map for the proc to send data
      ALLOCATE (map_send_size(0:para_env_RPA%num_pe-1))
      map_send_size = 0
      DO jjB = 1, ncol_local_s
         j_global = col_indices_s(jjB)
         send_pcol = cp_fm_indxg2p(j_global, ncol_block_d, dummy_proc, &
                                   fm_mat_dest%matrix_struct%first_p_pos(2), npcol_d)
         DO iiB = 1, nrow_local_s
            i_global = row_indices_s(iiB)
            send_prow = cp_fm_indxg2p(i_global, nrow_block_d, dummy_proc, &
                                      fm_mat_dest%matrix_struct%first_p_pos(1), nprow_d)
            proc_send = fm_mat_dest%matrix_struct%context%blacs2mpi(send_prow, send_pcol)
            map_send_size(proc_send) = map_send_size(proc_send)+1
         END DO
      END DO

      ! 2) loop over my local data of fm_mat_S and define a map for the proc from which rec data
      ALLOCATE (map_rec_size(0:para_env_RPA%num_pe-1))
      map_rec_size = 0
      DO jjB = 1, ncol_local_d
         j_global = col_indices_d(jjB)
         rec_pcol = cp_fm_indxg2p(j_global, ncol_block_s, dummy_proc, &
                                  fm_mat_source%matrix_struct%first_p_pos(2), npcol_s)
         DO iiB = 1, nrow_local_d
            i_global = row_indices_d(iiB)
            rec_prow = cp_fm_indxg2p(i_global, nrow_block_s, dummy_proc, &
                                     fm_mat_source%matrix_struct%first_p_pos(1), nprow_s)
            proc_receive = fm_mat_source%matrix_struct%context%blacs2mpi(rec_prow, rec_pcol)
            map_rec_size(proc_receive) = map_rec_size(proc_receive)+1
         END DO
      END DO

      ! 3) calculate the number of messages to send and allocate the send buffer
      number_of_send = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send+1
         END IF
      END DO

      ALLOCATE (buffer_send(number_of_send))

      ! 3.5) prepare the index map
      CALL timeset(routineN//"_bS", handle2)
      ALLOCATE (proc2counter(0:para_env_RPA%num_pe-1))
      proc2counter = 0
      ! allocate buffer for sending
      send_counter = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
         size_send_buffer = map_send_size(proc_send)

         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter+1
            ! prepare the sending buffer
            ALLOCATE (buffer_send(send_counter)%indx(2, size_send_buffer))
            buffer_send(send_counter)%indx = 0

            proc2counter(proc_send) = send_counter
         END IF
      END DO

      ALLOCATE (index_counter(0:para_env_RPA%num_pe-1))
      index_counter = 0
      DO iiB = 1, nrow_local_s
         i_global = row_indices_s(iiB)
         send_prow = cp_fm_indxg2p(i_global, nrow_block_d, dummy_proc, &
                                   fm_mat_dest%matrix_struct%first_p_pos(1), nprow_d)
         DO jjB = 1, ncol_local_s
            j_global = col_indices_s(jjB)
            send_pcol = cp_fm_indxg2p(j_global, ncol_block_d, dummy_proc, &
                                      fm_mat_dest%matrix_struct%first_p_pos(2), npcol_d)
            iii = fm_mat_dest%matrix_struct%context%blacs2mpi(send_prow, send_pcol)
            IF (iii == para_env_RPA%mepos) CYCLE
            index_counter(iii) = index_counter(iii)+1
            send_counter = proc2counter(iii)
            buffer_send(send_counter)%indx(1, index_counter(iii)) = iiB
            buffer_send(send_counter)%indx(2, index_counter(iii)) = jjB
            ! buffer_send(send_counter)%msg(index_counter(iii))=fm_mat_source%local_data(iiB,jjB)
         END DO
      END DO

      DEALLOCATE (index_counter)
      CALL timestop(handle2)

      ! 4) calculate the number of messages to receive and allocate the rec buffer
      number_of_rec = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec+1
         END IF
      END DO

      ALLOCATE (buffer_rec(number_of_rec))

      ! 4.5) prepare the index map
      CALL timeset(routineN//"_bR", handle2)
      rec_counter = 0
      proc2counter = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)

         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter+1
            ! allocate the auxilliary index structure
            ALLOCATE (buffer_rec(rec_counter)%indx(2, size_rec_buffer))
            buffer_rec(rec_counter)%indx = 0

            proc2counter(proc_receive) = rec_counter
         END IF
      END DO

!$OMP     PARALLEL DO DEFAULT(NONE) SHARED(para_env_RPA,RPA_proc_map,map_rec_size,&
!$OMP                                         local_size_source,fm_mat_source,fm_mat_dest,&
!$OMP                                         nrow_block_s,nprow_s,nrow_block_d,nprow_d,dummy_proc,&
!$OMP                                         ncol_block_s,npcol_s,ncol_block_d,npcol_d,&
!$OMP                                         myprow_d,mypcol_d,&
!$OMP                                         buffer_rec,proc2counter)&
!$OMP     PRIVATE(proc_receive,size_rec_buffer,rec_local_row,rec_local_col,ref_rec_prow,ref_rec_pcol,&
!$OMP             rec_counter,iii,i_global,rec_prow,i_local,j_global,rec_pcol,j_local)
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)

         rec_local_row = local_size_source(1, proc_receive)
         rec_local_col = local_size_source(2, proc_receive)

         ref_rec_prow = fm_mat_source%matrix_struct%context%mpi2blacs(1, proc_receive)
         ref_rec_pcol = fm_mat_source%matrix_struct%context%mpi2blacs(2, proc_receive)

         IF (map_rec_size(proc_receive) > 0) THEN

            rec_counter = proc2counter(proc_receive)
            iii = 0
            DO iiB = 1, rec_local_row
               i_global = cp_fm_indxl2g(iiB, nrow_block_s, ref_rec_prow, &
                                        fm_mat_source%matrix_struct%first_p_pos(1), nprow_s)
               rec_prow = cp_fm_indxg2p(i_global, nrow_block_d, dummy_proc, &
                                        fm_mat_dest%matrix_struct%first_p_pos(1), nprow_d)
               IF (rec_prow == myprow_d) THEN
                  i_local = cp_fm_indxg2l(i_global, nrow_block_d, dummy_proc, &
                                          fm_mat_dest%matrix_struct%first_p_pos(1), nprow_d)
                  DO jjB = 1, rec_local_col
                     j_global = cp_fm_indxl2g(jjB, ncol_block_s, ref_rec_pcol, &
                                              fm_mat_source%matrix_struct%first_p_pos(2), npcol_s)
                     rec_pcol = cp_fm_indxg2p(j_global, ncol_block_d, dummy_proc, &
                                              fm_mat_dest%matrix_struct%first_p_pos(2), npcol_d)
                     IF (rec_pcol == mypcol_d) THEN
                        j_local = cp_fm_indxg2l(j_global, ncol_block_d, dummy_proc, &
                                                fm_mat_dest%matrix_struct%first_p_pos(2), npcol_d)

                        iii = iii+1
                        buffer_rec(rec_counter)%indx(1, iii) = i_local
                        buffer_rec(rec_counter)%indx(2, iii) = j_local
                     END IF
                  END DO
               END IF
            END DO

         END IF
      END DO
!$OMP     END PARALLEL DO

      DEALLOCATE (proc2counter)
      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE initialize_buffer

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_source ...
!> \param fm_mat_dest ...
!> \param RPA_proc_map ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param number_of_send ...
!> \param map_send_size ...
!> \param map_rec_size ...
!> \param local_size_source ...
!> \param para_env_RPA ...
! **************************************************************************************************
   SUBROUTINE fm_redistribute(fm_mat_source, fm_mat_dest, RPA_proc_map, buffer_rec, buffer_send, &
                              number_of_send, &
                              map_send_size, map_rec_size, local_size_source, para_env_RPA)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_source, fm_mat_dest
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: RPA_proc_map
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      INTEGER                                            :: number_of_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: map_send_size, map_rec_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_size_source
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_redistribute', &
         routineP = moduleN//':'//routineN

      INTEGER :: dummy_proc, handle, handle2, i_global, i_local, iiB, iii, j_global, j_local, jjB, &
         mypcol_d, mypcol_s, myprow_d, myprow_s, ncol_block_d, ncol_block_s, ncol_local_d, &
         ncol_local_s, npcol_d, npcol_s, nprow_d, nprow_s, nrow_block_d, nrow_block_s, &
         nrow_local_d, nrow_local_s, proc_receive, proc_send, proc_shift, rec_counter, &
         rec_local_col, rec_local_row, rec_pcol, rec_prow, ref_rec_pcol, ref_rec_prow, &
         send_counter, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: req_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_d, col_indices_s, &
                                                            row_indices_d, row_indices_s

      CALL timeset(routineN, handle)

      ! get info for the source
      CALL cp_fm_get_info(matrix=fm_mat_source, &
                          nrow_local=nrow_local_s, &
                          ncol_local=ncol_local_s, &
                          row_indices=row_indices_s, &
                          col_indices=col_indices_s, &
                          nrow_block=nrow_block_s, &
                          ncol_block=ncol_block_s)
      myprow_s = fm_mat_source%matrix_struct%context%mepos(1)
      mypcol_s = fm_mat_source%matrix_struct%context%mepos(2)
      nprow_s = fm_mat_source%matrix_struct%context%num_pe(1)
      npcol_s = fm_mat_source%matrix_struct%context%num_pe(2)

      ! get info for the dest
      CALL cp_fm_get_info(matrix=fm_mat_dest, &
                          nrow_local=nrow_local_d, &
                          ncol_local=ncol_local_d, &
                          row_indices=row_indices_d, &
                          col_indices=col_indices_d, &
                          nrow_block=nrow_block_d, &
                          ncol_block=ncol_block_d)
      myprow_d = fm_mat_dest%matrix_struct%context%mepos(1)
      mypcol_d = fm_mat_dest%matrix_struct%context%mepos(2)
      nprow_d = fm_mat_dest%matrix_struct%context%num_pe(1)
      npcol_d = fm_mat_dest%matrix_struct%context%num_pe(2)

      ! 0) check if the local data has to be stored in the new fm_mat_S
      CALL timeset(routineN//"_loc", handle2)
      IF (map_rec_size(para_env_RPA%mepos) > 0) THEN
         DO jjB = 1, ncol_local_s
            j_global = col_indices_s(jjB)
            rec_pcol = cp_fm_indxg2p(j_global, ncol_block_d, dummy_proc, &
                                     fm_mat_dest%matrix_struct%first_p_pos(2), npcol_d)
            IF (rec_pcol == mypcol_d) THEN
               j_local = cp_fm_indxg2l(j_global, ncol_block_d, dummy_proc, &
                                       fm_mat_dest%matrix_struct%first_p_pos(2), npcol_d)
               DO iiB = 1, nrow_local_s
                  i_global = row_indices_s(iiB)
                  rec_prow = cp_fm_indxg2p(i_global, nrow_block_d, dummy_proc, &
                                           fm_mat_dest%matrix_struct%first_p_pos(1), nprow_d)
                  IF (rec_prow == myprow_d) THEN
                     i_local = cp_fm_indxg2l(i_global, nrow_block_d, dummy_proc, &
                                             fm_mat_dest%matrix_struct%first_p_pos(1), nprow_d)
                     fm_mat_dest%local_data(i_local, j_local) = fm_mat_source%local_data(iiB, jjB)
                  END IF
               END DO
            END IF
         END DO
      END IF
      CALL timestop(handle2)

      ! 1) prepare receiving buffer
      CALL timeset(routineN//"_post", handle2)
      rec_counter = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)

         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter+1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%msg = 0.0_dp
            buffer_rec(rec_counter)%proc = proc_receive

            ! post the receiving message
            CALL mp_irecv(buffer_rec(rec_counter)%msg, proc_receive, para_env_RPA%group, buffer_rec(rec_counter)%msg_req)

         END IF
      END DO
      CALL timestop(handle2)

      ! 2) prepare sending buffer
      CALL timeset(routineN//"_bS", handle2)
      ! allocate buffer for sending, fill the buffer, send the message
      ALLOCATE (req_send(number_of_send))
      send_counter = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
         size_send_buffer = map_send_size(proc_send)

         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter+1
            ! allocate
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%msg = 0.0_dp
            buffer_send(send_counter)%proc = proc_send

            ! fill
!$OMP         PARALLEL DO DEFAULT(NONE) PRIVATE(iii,iiB,jjB) &
!$OMP                                   SHARED(size_send_buffer,buffer_send,&
!$OMP                                          send_counter,fm_mat_source)
            DO iii = 1, size_send_buffer
               iiB = buffer_send(send_counter)%indx(1, iii)
               jjB = buffer_send(send_counter)%indx(2, iii)
               buffer_send(send_counter)%msg(iii) = fm_mat_source%local_data(iiB, jjB)
            END DO
!$OMP         END PARALLEL DO

            ! send
            CALL mp_isend(buffer_send(send_counter)%msg, proc_send, para_env_RPA%group, buffer_send(send_counter)%msg_req)
            req_send(send_counter) = buffer_send(send_counter)%msg_req

         END IF
      END DO
      CALL timestop(handle2)

      ! 3) fill the fm_mat_dest matrix with the received data
      CALL timeset(routineN//"_fill", handle2)
      rec_counter = 0
      DO proc_shift = 1, para_env_RPA%num_pe-1
         proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)

         rec_local_row = local_size_source(1, proc_receive)
         rec_local_col = local_size_source(2, proc_receive)

         ref_rec_prow = fm_mat_source%matrix_struct%context%mpi2blacs(1, proc_receive)
         ref_rec_pcol = fm_mat_source%matrix_struct%context%mpi2blacs(2, proc_receive)

         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter+1

            CALL mp_wait(buffer_rec(rec_counter)%msg_req)

            ! fill the destination matrix
!$OMP         PARALLEL DO DEFAULT(NONE) PRIVATE(iii,i_local,j_local) &
!$OMP                                   SHARED(size_rec_buffer,proc_receive,&
!$OMP                                          rec_counter,fm_mat_dest,buffer_rec)
            DO iii = 1, size_rec_buffer
               i_local = buffer_rec(rec_counter)%indx(1, iii)
               j_local = buffer_rec(rec_counter)%indx(2, iii)
               fm_mat_dest%local_data(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
            END DO
!$OMP         END PARALLEL DO

            DEALLOCATE (buffer_rec(rec_counter)%msg)

         END IF
      END DO
      CALL timestop(handle2)

      ! 5) wait all
      CALL timeset(routineN//"_wA", handle2)
      CALL mp_waitall(req_send(:))
      DEALLOCATE (req_send)
      DO send_counter = 1, number_of_send
         DEALLOCATE (buffer_send(send_counter)%msg)
      END DO
      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE fm_redistribute

! **************************************************************************************************
!> \brief ...
!> \param RPA_proc_map ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param number_of_rec ...
!> \param number_of_send ...
!> \param map_send_size ...
!> \param map_rec_size ...
!> \param local_size_source ...
! **************************************************************************************************
   SUBROUTINE release_buffer(RPA_proc_map, buffer_rec, buffer_send, &
                             number_of_rec, number_of_send, &
                             map_send_size, map_rec_size, local_size_source)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: RPA_proc_map
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      INTEGER                                            :: number_of_rec, number_of_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: map_send_size, map_rec_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_size_source

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_buffer', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, rec_counter, send_counter

      CALL timeset(routineN, handle)

      IF (ALLOCATED(RPA_proc_map)) DEALLOCATE (RPA_proc_map)
      IF (ALLOCATED(map_send_size)) DEALLOCATE (map_send_size)
      IF (ALLOCATED(map_rec_size)) DEALLOCATE (map_rec_size)
      IF (ALLOCATED(local_size_source)) DEALLOCATE (local_size_source)

      IF (ALLOCATED(buffer_send)) THEN
         DO send_counter = 1, number_of_send
            IF (ASSOCIATED(buffer_send(send_counter)%msg)) DEALLOCATE (buffer_send(send_counter)%msg)
            IF (ASSOCIATED(buffer_send(send_counter)%indx)) DEALLOCATE (buffer_send(send_counter)%indx)
         END DO
         DEALLOCATE (buffer_send)
      END IF

      IF (ALLOCATED(buffer_rec)) THEN
         DO rec_counter = 1, number_of_rec
            IF (ASSOCIATED(buffer_rec(rec_counter)%msg)) DEALLOCATE (buffer_rec(rec_counter)%msg)
            IF (ASSOCIATED(buffer_rec(rec_counter)%indx)) DEALLOCATE (buffer_rec(rec_counter)%indx)
         END DO
         DEALLOCATE (buffer_rec)
      END IF

      CALL timestop(handle)

   END SUBROUTINE release_buffer

! **************************************************************************************************
!> \brief Redistribute RPA-AXK Gamma_3 density matrices: from fm to dbcsr
!> \param fm_mat_Gamma_3 ... ia*dime_RI sized density matrix (fm type on para_env_RPA)
!> \param dbcsr_Gamma_3 ...  redistributed Gamma_3 (dbcsr array): dimen_RI of i*a: i*a on subgroup, L distributed in RPA_group
!> \param para_env_RPA ...
!> \param para_env_sub ...
!> \param homo ...
!> \param virtual ...
!> \param mo_coeff_o ...   dbcsr on a subgroup
!> \param ngroup ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param dimen_RI ...
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE gamma_fm_to_dbcsr(fm_mat_Gamma_3, dbcsr_Gamma_3, para_env_RPA, para_env_sub, &
                                homo, virtual, mo_coeff_o, ngroup, my_group_L_start, my_group_L_end, &
                                my_group_L_size, dimen_RI)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Gamma_3
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: dbcsr_Gamma_3
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA, para_env_sub
      INTEGER                                            :: homo, virtual
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_o
      INTEGER                                            :: ngroup, my_group_L_start, &
                                                            my_group_L_end, my_group_L_size, &
                                                            dimen_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_fm_to_dbcsr', &
         routineP = moduleN//':'//routineN

      INTEGER :: dimen_ia, dummy_proc, handle, i, i_global, i_local, iaia, iib, iii, iproc, &
         itmp(2), j_global, j_local, jjb, jjj, kkb, my_ia_end, my_ia_size, my_ia_start, mypcol, &
         myprow, ncol_block, ncol_local, npcol, nprow, nrow_block, nrow_local, number_of_rec, &
         number_of_send, proc_receive, proc_send, proc_shift, rec_counter, rec_iaia_end, &
         rec_iaia_size, rec_iaia_start, rec_pcol, rec_prow, ref_send_pcol, ref_send_prow, &
         send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_ia, iii_vet, map_rec_size, &
                                                            map_send_size, proc_map, req_send, &
                                                            sizes_ia, starts_ia, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: grid_2_mepos, grid_ref_2_send_pos, &
         group_grid_2_mepos, indices_map_my, mepos_2_grid, mepos_2_grid_group, sizes
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: part_ia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Gamma_2D
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_ia
      TYPE(index_map), ALLOCATABLE, DIMENSION(:)         :: indices_rec
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      dimen_ia = virtual*homo

      ! Prepare sizes for a 2D array
      ALLOCATE (sizes_ia(0:para_env_sub%num_pe-1))
      sizes_ia = 0
      ALLOCATE (starts_ia(0:para_env_sub%num_pe-1))
      starts_ia = 0
      ALLOCATE (ends_ia(0:para_env_sub%num_pe-1))
      ends_ia = 0

      DO iproc = 0, para_env_sub%num_pe-1
         itmp = get_limit(dimen_ia, para_env_sub%num_pe, iproc)
         starts_ia(iproc) = itmp(1)
         ends_ia(iproc) = itmp(2)
         sizes_ia(iproc) = itmp(2)-itmp(1)+1
      END DO

      my_ia_start = starts_ia(para_env_sub%mepos)
      my_ia_end = ends_ia(para_env_sub%mepos)
      my_ia_size = sizes_ia(para_env_sub%mepos)

      ! Make a 2D array intermediate

      CALL prepare_redistribution(para_env_RPA, para_env_sub, ngroup, &
                                  proc_map, sizes, group_grid_2_mepos, mepos_2_grid_group, &
                                  fm_mat_Gamma_3)

      CALL fm2array(Gamma_2D, dimen_ia, dimen_RI, para_env_RPA, proc_map, &
                    my_ia_size, my_ia_start, my_ia_end, &
                    my_group_L_size, my_group_L_start, my_group_L_end, &
                    sizes, group_grid_2_mepos, mepos_2_grid_group, &
                    para_env_sub%num_pe, ngroup, &
                    fm_mat_Gamma_3)

      ! create sub blacs env
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_sub)

      ! create the fm_ia buffer matrix
      NULLIFY (fm_ia)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=homo, &
                               ncol_global=virtual, para_env=para_env_sub)
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
                          col_indices=col_indices, &
                          nrow_block=nrow_block, &
                          ncol_block=ncol_block)
      myprow = fm_ia%matrix_struct%context%mepos(1)
      mypcol = fm_ia%matrix_struct%context%mepos(2)
      nprow = fm_ia%matrix_struct%context%num_pe(1)
      npcol = fm_ia%matrix_struct%context%num_pe(2)

      ! 0) create array containing the processes position and supporting infos
      ALLOCATE (grid_2_mepos(0:nprow-1, 0:npcol-1))
      grid_2_mepos = 0
      ALLOCATE (mepos_2_grid(0:para_env_sub%num_pe-1, 2))
      mepos_2_grid = 0
      ! fill the info array
      grid_2_mepos(myprow, mypcol) = para_env_sub%mepos
      mepos_2_grid(para_env_sub%mepos, 1) = myprow
      mepos_2_grid(para_env_sub%mepos, 2) = mypcol
      ! sum infos
      CALL mp_sum(grid_2_mepos, para_env_sub%group)
      CALL mp_sum(mepos_2_grid, para_env_sub%group)

      ! loop over local index range and define the sending map
      ALLOCATE (map_send_size(0:para_env_sub%num_pe-1))
      map_send_size = 0
      DO iaia = my_ia_start, my_ia_end
         i_global = (iaia-1)/virtual+1
         j_global = MOD(iaia-1, virtual)+1
         send_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                   fm_ia%matrix_struct%first_p_pos(1), nprow)
         send_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                   fm_ia%matrix_struct%first_p_pos(2), npcol)
         proc_send = grid_2_mepos(send_prow, send_pcol)
         map_send_size(proc_send) = map_send_size(proc_send)+1
      END DO

      ! loop over local data of fm_ia and define the receiving map
      ALLOCATE (map_rec_size(0:para_env_sub%num_pe-1))
      map_rec_size = 0
      part_ia = REAL(dimen_ia, KIND=dp)/REAL(para_env_sub%num_pe, KIND=dp)

      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            iaia = (i_global-1)*virtual+j_global
            proc_receive = INT(REAL(iaia-1, KIND=dp)/part_ia)
            proc_receive = MAX(0, proc_receive)
            proc_receive = MIN(proc_receive, para_env_sub%num_pe-1)
            DO
               itmp = get_limit(dimen_ia, para_env_sub%num_pe, proc_receive)
               IF (iaia >= itmp(1) .AND. iaia <= itmp(2)) EXIT
               IF (iaia < itmp(1)) proc_receive = proc_receive-1
               IF (iaia > itmp(2)) proc_receive = proc_receive+1
            END DO
            map_rec_size(proc_receive) = map_rec_size(proc_receive)+1
         END DO
      END DO

      ! create the sub_proc_map
      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
      DO i = 0, para_env_sub%num_pe-1
         sub_proc_map(i) = i
         sub_proc_map(-i-1) = para_env_sub%num_pe-i-1
         sub_proc_map(para_env_sub%num_pe+i) = i
      END DO
      ! allocate the buffer for sending data
      number_of_send = 0
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send+1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be sent
      ALLOCATE (buffer_send(number_of_send))
      ! and the map from the grid of processess to the message position
      ALLOCATE (grid_ref_2_send_pos(0:nprow-1, 0:npcol-1))
      grid_ref_2_send_pos = 0
      ! finally allocate each message
      send_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
         size_send_buffer = map_send_size(proc_send)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter+1
            ! allocate the sending buffer (msg)
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%proc = proc_send
            ! get the pointer to prow, pcol of the process that has
            ! to receive this message
            ref_send_prow = mepos_2_grid(proc_send, 1)
            ref_send_pcol = mepos_2_grid(proc_send, 2)
            ! save the rank of the process that has to receive this message
            grid_ref_2_send_pos(ref_send_prow, ref_send_pcol) = send_counter
         END IF
      END DO

      ! allocate the buffer for receiving data
      number_of_rec = 0
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec+1
         END IF
      END DO

      ! allocate the structure that will hold the messages to be received
      ! and relative indeces
      ALLOCATE (buffer_rec(number_of_rec))
      ALLOCATE (indices_rec(number_of_rec))
      ! finally allocate each message and fill the array of indeces
      rec_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter+1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%proc = proc_receive
            ! create the indices array
            ALLOCATE (indices_rec(rec_counter)%map(2, size_rec_buffer))
            indices_rec(rec_counter)%map = 0
            rec_iaia_start = starts_ia(proc_receive)
            rec_iaia_end = ends_ia(proc_receive)
            rec_iaia_size = sizes_ia(proc_receive)
            iii = 0
            DO iaia = rec_iaia_start, rec_iaia_end
               i_global = (iaia-1)/virtual+1
               j_global = MOD(iaia-1, virtual)+1
               rec_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                        fm_ia%matrix_struct%first_p_pos(1), nprow)
               rec_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                        fm_ia%matrix_struct%first_p_pos(2), npcol)
               IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
               iii = iii+1
               i_local = cp_fm_indxg2l(i_global, nrow_block, dummy_proc, &
                                       fm_ia%matrix_struct%first_p_pos(1), nprow)
               j_local = cp_fm_indxg2l(j_global, ncol_block, dummy_proc, &
                                       fm_ia%matrix_struct%first_p_pos(2), npcol)
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
         DO iaia = my_ia_start, my_ia_end
            i_global = (iaia-1)/virtual+1
            j_global = MOD(iaia-1, virtual)+1
            rec_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                     fm_ia%matrix_struct%first_p_pos(1), nprow)
            rec_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                     fm_ia%matrix_struct%first_p_pos(2), npcol)
            IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
            iii = iii+1
            i_local = cp_fm_indxg2l(i_global, nrow_block, dummy_proc, &
                                    fm_ia%matrix_struct%first_p_pos(1), nprow)
            j_local = cp_fm_indxg2l(j_global, ncol_block, dummy_proc, &
                                    fm_ia%matrix_struct%first_p_pos(2), npcol)
            indices_map_my(1, iii) = i_local
            indices_map_my(2, iii) = j_local
         END DO
      END IF

      ! Allocate dbcsr_Gamma_3
      NULLIFY (dbcsr_Gamma_3)

      !CALL cp_dbcsr_allocate_matrix_set(dbcsr_Gamma_3, ncol_local)
      CALL cp_dbcsr_allocate_matrix_set(dbcsr_Gamma_3, my_group_L_size)

      ! auxiliary vector of indices for the send buffer
      ALLOCATE (iii_vet(number_of_send))
      ! vector for the send requests
      ALLOCATE (req_send(number_of_send))
      ! loop over auxiliary basis function and redistribute into a fm
      ! and then compy the fm into a dbcsr matrix

      !DO kkB = 1, ncol_local
      DO kkB = 1, my_group_L_size
         ! zero the matries of the buffers and post the messages to be received
         CALL cp_fm_set_all(matrix=fm_ia, alpha=0.0_dp)
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter+1
               buffer_rec(rec_counter)%msg = 0.0_dp
               CALL mp_irecv(buffer_rec(rec_counter)%msg, proc_receive, para_env_sub%group, &
                             buffer_rec(rec_counter)%msg_req)
            END IF
         END DO
         ! fill the sending buffer and send the messages
         DO send_counter = 1, number_of_send
            buffer_send(send_counter)%msg = 0.0_dp
         END DO
         iii_vet = 0
         jjj = 0
         DO iaia = my_ia_start, my_ia_end
            i_global = (iaia-1)/virtual+1
            j_global = MOD(iaia-1, virtual)+1
            send_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                      fm_ia%matrix_struct%first_p_pos(1), nprow)
            send_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                      fm_ia%matrix_struct%first_p_pos(2), npcol)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            ! we don't need to send to ourselves
            IF (grid_2_mepos(send_prow, send_pcol) == para_env_sub%mepos) THEN
               ! filling fm_ia with local data
               jjj = jjj+1
               i_local = indices_map_my(1, jjj)
               j_local = indices_map_my(2, jjj)
               fm_ia%local_data(i_local, j_local) = &
                  Gamma_2D(iaia-my_ia_start+1, kkB)

            ELSE
               send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
               iii_vet(send_counter) = iii_vet(send_counter)+1
               iii = iii_vet(send_counter)
               buffer_send(send_counter)%msg(iii) = &
                  Gamma_2D(iaia-my_ia_start+1, kkB)
            END IF
         END DO
         req_send = 0
         send_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
            IF (map_send_size(proc_send) > 0) THEN
               send_counter = send_counter+1
               CALL mp_isend(buffer_send(send_counter)%msg, proc_send, para_env_sub%group, &
                             buffer_send(send_counter)%msg_req)
               req_send(send_counter) = buffer_send(send_counter)%msg_req
            END IF
         END DO

         ! receive the messages and fill the fm_ia
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)
            size_rec_buffer = map_rec_size(proc_receive)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter+1
               ! wait for the message
               CALL mp_wait(buffer_rec(rec_counter)%msg_req)
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
         ALLOCATE (dbcsr_Gamma_3(kkB)%matrix)
         CALL cp_dbcsr_init(dbcsr_Gamma_3(kkB)%matrix)
         CALL cp_dbcsr_m_by_n_from_template(dbcsr_Gamma_3(kkB)%matrix, &
                                            template=mo_coeff_o, &
                                            m=homo, n=virtual, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_ia, dbcsr_Gamma_3(kkB)%matrix, keep_sparsity=.FALSE.)

      END DO

      ! Deallocate memory

      DEALLOCATE (Gamma_2d)
      DEALLOCATE (iii_vet)
      DEALLOCATE (req_send)
      IF (map_rec_size(para_env_sub%mepos) > 0) THEN
         DEALLOCATE (indices_map_my)
      END IF
      DO rec_counter = 1, number_of_rec
         DEALLOCATE (indices_rec(rec_counter)%map)
         DEALLOCATE (buffer_rec(rec_counter)%msg)
      END DO
      DEALLOCATE (indices_rec)
      DEALLOCATE (buffer_rec)
      DO send_counter = 1, number_of_send
         DEALLOCATE (buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE (buffer_send)
      DEALLOCATE (map_send_size)
      DEALLOCATE (map_rec_size)
      DEALLOCATE (sub_proc_map)
      DEALLOCATE (grid_2_mepos)
      DEALLOCATE (mepos_2_grid)

      ! release buffer matrix
      CALL cp_fm_release(fm_ia)
      ! release auxiliary gamma
      CALL cp_fm_release(fm_mat_Gamma_3)

      CALL timestop(handle)

   END SUBROUTINE gamma_fm_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ngroup ...
!> \param proc_map ...
!> \param sizes ...
!> \param group_grid_2_mepos ...
!> \param mepos_2_grid_group ...
!> \param fm_mat_Gamma_3 ...
! **************************************************************************************************
   SUBROUTINE prepare_redistribution(para_env, para_env_sub, ngroup, &
                                     proc_map, sizes, group_grid_2_mepos, mepos_2_grid_group, &
                                     fm_mat_Gamma_3)
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER                                            :: ngroup
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: sizes, group_grid_2_mepos, &
                                                            mepos_2_grid_group
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Gamma_3

      INTEGER                                            :: i, ncol_local, nrow_local, pos_group, &
                                                            pos_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pos_info

      ALLOCATE (proc_map(-para_env%num_pe:2*para_env%num_pe-1))
      proc_map = 0
      ALLOCATE (pos_info(0:para_env%num_pe-1))
      pos_info = 0
      pos_info(para_env%mepos) = para_env_sub%mepos
      CALL mp_sum(pos_info, para_env%group)

      ALLOCATE (group_grid_2_mepos(0:para_env_sub%num_pe-1, 0:ngroup-1))
      group_grid_2_mepos = 0
      ALLOCATE (mepos_2_grid_group(0:para_env%num_pe-1, 2))
      mepos_2_grid_group = 0

      DO i = 0, para_env%num_pe-1
         proc_map(i) = i
         proc_map(-i-1) = para_env%num_pe-i-1
         proc_map(para_env%num_pe+i) = i
         ! calculate postition of the group
         pos_group = i/para_env_sub%num_pe
         ! calculate postition in the subgroup
         pos_sub = pos_info(i)
         ! fill the map from the grid of groups to process
         group_grid_2_mepos(pos_sub, pos_group) = i
         ! and the opposite, from the global pos to the grid pos
         mepos_2_grid_group(i, 1) = pos_sub
         mepos_2_grid_group(i, 2) = pos_group
      ENDDO

      CALL cp_fm_get_info(matrix=fm_mat_Gamma_3, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local)

      ALLOCATE (sizes(2, 0:para_env%num_pe-1))
      sizes = 0
      sizes(1, para_env%mepos) = nrow_local
      sizes(2, para_env%mepos) = ncol_local
      CALL mp_sum(sizes, para_env%group)

   END SUBROUTINE prepare_redistribution

END MODULE rpa_communication

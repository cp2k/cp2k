!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines necessary to redistribute an fm_matrix from a
!>        given blacs_env to another
!> \par History
!>      12.2012 created [Mauro Del Ben]
! **************************************************************************************************
MODULE rpa_communication
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_template,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_indxg2l,&
                                              cp_fm_indxg2p,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              release_group_dist
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_wait,&
                                              mp_waitall
   USE mp2_ri_grad_util,                ONLY: fm2array,&
                                              prepare_redistribution
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE index_map
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: map
   END TYPE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_communication'

   PUBLIC :: gamma_fm_to_dbcsr, &
             communicate_buffer

CONTAINS

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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_3
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA, para_env_sub
      INTEGER, INTENT(IN)                                :: homo, virtual
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o
      INTEGER, INTENT(IN)                                :: ngroup, my_group_L_start, &
                                                            my_group_L_end, my_group_L_size, &
                                                            dimen_RI

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gamma_fm_to_dbcsr'

      INTEGER :: dimen_ia, dummy_proc, handle, i, i_global, i_local, iaia, iib, iii, itmp(2), &
         j_global, j_local, jjb, jjj, kkb, my_ia_end, my_ia_size, my_ia_start, mypcol, myprow, &
         ncol_block, ncol_local, npcol, nprow, nrow_block, nrow_local, number_of_rec, &
         number_of_send, proc_receive, proc_send, proc_shift, rec_counter, rec_iaia_end, &
         rec_iaia_size, rec_iaia_start, rec_pcol, rec_prow, ref_send_pcol, ref_send_prow, &
         send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii_vet, map_rec_size, map_send_size, &
                                                            req_send, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: grid_2_mepos, grid_ref_2_send_pos, &
                                                            group_grid_2_mepos, indices_map_my, &
                                                            mepos_2_grid, mepos_2_grid_group
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: part_ia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Gamma_2D
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_ia
      TYPE(group_dist_d1_type)                           :: gd_ia
      TYPE(index_map), ALLOCATABLE, DIMENSION(:)         :: indices_rec
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      dimen_ia = virtual*homo

      ! Prepare sizes for a 2D array
      CALL create_group_dist(gd_ia, para_env_sub%num_pe, dimen_ia)
      CALL get_group_dist(gd_ia, para_env_sub%mepos, my_ia_start, my_ia_end, my_ia_size)

      ! Make a 2D array intermediate

      CALL prepare_redistribution(para_env_RPA, para_env_sub, ngroup, &
                                  group_grid_2_mepos, mepos_2_grid_group)

      CALL fm2array(Gamma_2D, dimen_ia, dimen_RI, para_env_RPA, &
                    my_ia_size, my_ia_start, my_ia_end, &
                    my_group_L_size, my_group_L_start, my_group_L_end, &
                    group_grid_2_mepos, mepos_2_grid_group, &
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
      ALLOCATE (grid_2_mepos(0:nprow - 1, 0:npcol - 1))
      grid_2_mepos = 0
      ALLOCATE (mepos_2_grid(2, 0:para_env_sub%num_pe - 1))
      ! fill the info array
      grid_2_mepos(myprow, mypcol) = para_env_sub%mepos
      ! sum infos
      CALL mp_sum(grid_2_mepos, para_env_sub%group)
      CALL mp_allgather([myprow, mypcol], mepos_2_grid, para_env_sub%group)

      ! loop over local index range and define the sending map
      ALLOCATE (map_send_size(0:para_env_sub%num_pe - 1))
      map_send_size = 0
      dummy_proc = 0
      DO iaia = my_ia_start, my_ia_end
         i_global = (iaia - 1)/virtual + 1
         j_global = MOD(iaia - 1, virtual) + 1
         send_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                   fm_ia%matrix_struct%first_p_pos(1), nprow)
         send_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                   fm_ia%matrix_struct%first_p_pos(2), npcol)
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
            iaia = (i_global - 1)*virtual + j_global
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

      ! create the sub_proc_map
      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe - 1))
      DO i = 0, para_env_sub%num_pe - 1
         sub_proc_map(i) = i
         sub_proc_map(-i - 1) = para_env_sub%num_pe - i - 1
         sub_proc_map(para_env_sub%num_pe + i) = i
      END DO
      ! allocate the buffer for sending data
      number_of_send = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
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
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
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
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec + 1
         END IF
      END DO

      ! allocate the structure that will hold the messages to be received
      ! and relative indeces
      ALLOCATE (buffer_rec(number_of_rec))
      ALLOCATE (indices_rec(number_of_rec))
      ! finally allocate each message and fill the array of indeces
      rec_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%proc = proc_receive
            ! create the indices array
            ALLOCATE (indices_rec(rec_counter)%map(2, size_rec_buffer))
            indices_rec(rec_counter)%map = 0
            CALL get_group_dist(gd_ia, proc_receive, rec_iaia_start, rec_iaia_end, rec_iaia_size)
            iii = 0
            DO iaia = rec_iaia_start, rec_iaia_end
               i_global = (iaia - 1)/virtual + 1
               j_global = MOD(iaia - 1, virtual) + 1
               rec_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                        fm_ia%matrix_struct%first_p_pos(1), nprow)
               rec_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                        fm_ia%matrix_struct%first_p_pos(2), npcol)
               IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
               iii = iii + 1
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
            i_global = (iaia - 1)/virtual + 1
            j_global = MOD(iaia - 1, virtual) + 1
            rec_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                     fm_ia%matrix_struct%first_p_pos(1), nprow)
            rec_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                     fm_ia%matrix_struct%first_p_pos(2), npcol)
            IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
            iii = iii + 1
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

      !CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_3, ncol_local)
      CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_3, my_group_L_size)

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
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter + 1
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
            i_global = (iaia - 1)/virtual + 1
            j_global = MOD(iaia - 1, virtual) + 1
            send_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                      fm_ia%matrix_struct%first_p_pos(1), nprow)
            send_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                      fm_ia%matrix_struct%first_p_pos(2), npcol)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            ! we don't need to send to ourselves
            IF (grid_2_mepos(send_prow, send_pcol) == para_env_sub%mepos) THEN
               ! filling fm_ia with local data
               jjj = jjj + 1
               i_local = indices_map_my(1, jjj)
               j_local = indices_map_my(2, jjj)
               fm_ia%local_data(i_local, j_local) = &
                  Gamma_2D(iaia - my_ia_start + 1, kkB)

            ELSE
               send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
               iii_vet(send_counter) = iii_vet(send_counter) + 1
               iii = iii_vet(send_counter)
               buffer_send(send_counter)%msg(iii) = &
                  Gamma_2D(iaia - my_ia_start + 1, kkB)
            END IF
         END DO
         req_send = 0
         send_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
            IF (map_send_size(proc_send) > 0) THEN
               send_counter = send_counter + 1
               CALL mp_isend(buffer_send(send_counter)%msg, proc_send, para_env_sub%group, &
                             buffer_send(send_counter)%msg_req)
               req_send(send_counter) = buffer_send(send_counter)%msg_req
            END IF
         END DO

         ! receive the messages and fill the fm_ia
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
            size_rec_buffer = map_rec_size(proc_receive)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter + 1
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
      CALL release_group_dist(gd_ia)

      ! release buffer matrix
      CALL cp_fm_release(fm_ia)
      ! release auxiliary gamma
      CALL cp_fm_release(fm_mat_Gamma_3)

      CALL timestop(handle)

   END SUBROUTINE gamma_fm_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param num_entries_rec ...
!> \param num_entries_send ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param req_array ...
!> \param do_indx ...
!> \param do_msg ...
! **************************************************************************************************
   SUBROUTINE communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, &
                                 req_array, do_indx, do_msg)

      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: num_entries_rec, num_entries_send
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: buffer_rec, buffer_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_indx, do_msg

      CHARACTER(LEN=*), PARAMETER :: routineN = 'communicate_buffer'

      INTEGER                                            :: handle, imepos, rec_counter, send_counter
      LOGICAL                                            :: my_do_indx, my_do_msg

      CALL timeset(routineN, handle)

      my_do_indx = .TRUE.
      IF (PRESENT(do_indx)) my_do_indx = do_indx
      my_do_msg = .TRUE.
      IF (PRESENT(do_msg)) my_do_msg = do_msg

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter + 1
               IF (my_do_indx) THEN
                  CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
               END IF
               IF (my_do_msg) THEN
                  CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
               END IF
            END IF
         END DO

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter + 1
               IF (my_do_indx) THEN
                  CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
               END IF
               IF (my_do_msg) THEN
                  CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
               END IF
            END IF
         END DO

         IF (my_do_indx) THEN
            CALL mp_waitall(req_array(1:send_counter, 1))
            CALL mp_waitall(req_array(1:rec_counter, 3))
         END IF

         IF (my_do_msg) THEN
            CALL mp_waitall(req_array(1:send_counter, 2))
            CALL mp_waitall(req_array(1:rec_counter, 4))
         END IF

      ELSE

         buffer_rec(0)%indx(:, :) = buffer_send(0)%indx
         buffer_rec(0)%msg(:) = buffer_send(0)%msg

      END IF

      CALL timestop(handle)

   END SUBROUTINE communicate_buffer

END MODULE rpa_communication

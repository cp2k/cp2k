!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Auxiliary rountines necessary to redistribute an fm_matrix from a
!>        given blacs_env to another
!> \par History
!>      12.2012 created [Mauro Del Ben]
! *****************************************************************************
MODULE rpa_communication
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_indxg2l,&
                                             cp_fm_indxg2p,&
                                             cp_fm_indxl2g,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_wait,&
                                             mp_waitall
  USE mp2_types,                       ONLY: integ_mat_buffer_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_communication'

  PUBLIC :: initialize_buffer,&
            fm_redistribute, &
            release_buffer

  CONTAINS

  SUBROUTINE initialize_buffer(fm_mat_source,fm_mat_dest,RPA_proc_map,buffer_rec,buffer_send,&
                               number_of_rec,number_of_send,&
                               map_send_size,map_rec_size,local_size_source,para_env_RPA,error)
    TYPE(cp_fm_type), POINTER                :: fm_mat_source, fm_mat_dest
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: RPA_proc_map
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send
    INTEGER                                  :: number_of_rec, number_of_send
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map_send_size, map_rec_size
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_size_source
    TYPE(cp_para_env_type), POINTER          :: para_env_RPA
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_buffer', &
      routineP = moduleN//':'//routineN

    INTEGER :: dummy_proc, handle, i, i_global, iiB, j_global, jjB, mypcol_d, &
      mypcol_s, myprow_d, myprow_s, ncol_block, ncol_block_d, ncol_block_s, &
      ncol_local, ncol_local_d, ncol_local_s, npcol_d, npcol_s, nprow_d, &
      nprow_s, nrow_block, nrow_block_d, nrow_block_s, nrow_local, &
      nrow_local_d, nrow_local_s, proc_receive, proc_send, proc_shift, &
      rec_pcol, rec_prow, send_pcol, send_prow, stat
    INTEGER, DIMENSION(:), POINTER           :: col_indices_d, col_indices_s, &
                                                row_indices_d, row_indices_s
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! create the RPA proc_map
    IF(.NOT.(ALLOCATED(RPA_proc_map))) THEN
      ALLOCATE(RPA_proc_map(-para_env_RPA%num_pe:2*para_env_RPA%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      RPA_proc_map=0
      DO i=0,para_env_RPA%num_pe-1
        RPA_proc_map(i)=i
        RPA_proc_map(-i-1)=para_env_RPA%num_pe-i-1
        RPA_proc_map(para_env_RPA%num_pe+i)=i
      END DO
    END IF

    ! get info for the source
    CALL cp_fm_get_info(matrix=fm_mat_source,&
                        nrow_local=nrow_local_s,&
                        ncol_local=ncol_local_s,&
                        row_indices=row_indices_s,&
                        col_indices=col_indices_s,&
                        nrow_block=nrow_block_s,&
                        ncol_block=ncol_block_s,&
                        error=error)
    myprow_s=fm_mat_source%matrix_struct%context%mepos(1)
    mypcol_s=fm_mat_source%matrix_struct%context%mepos(2)
    nprow_s =fm_mat_source%matrix_struct%context%num_pe(1)
    npcol_s =fm_mat_source%matrix_struct%context%num_pe(2)

    ! get info for the dest
    CALL cp_fm_get_info(matrix=fm_mat_dest,&
                        nrow_local=nrow_local_d,&
                        ncol_local=ncol_local_d,&
                        row_indices=row_indices_d,&
                        col_indices=col_indices_d,&
                        nrow_block=nrow_block_d,&
                        ncol_block=ncol_block_d,&
                        error=error)
    myprow_d=fm_mat_dest%matrix_struct%context%mepos(1)
    mypcol_d=fm_mat_dest%matrix_struct%context%mepos(2)
    nprow_d =fm_mat_dest%matrix_struct%context%num_pe(1)
    npcol_d =fm_mat_dest%matrix_struct%context%num_pe(2)

    ! 0) create the map for the local sizes
    ALLOCATE(local_size_source(2,0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    local_size_source=0
    local_size_source(1,para_env_RPA%mepos)=nrow_local_s
    local_size_source(2,para_env_RPA%mepos)=ncol_local_s
    CALL mp_sum(local_size_source,para_env_RPA%group)

    ! 1) loop over my local data and define a map for the proc to send data
    ALLOCATE(map_send_size(0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_send_size=0
    DO jjB=1, ncol_local_s
      j_global=col_indices_s(jjB)
      send_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
                              fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
      DO iiB=1, nrow_local_s
        i_global=row_indices_s(iiB)
        send_prow=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
                                fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
        proc_send=fm_mat_dest%matrix_struct%context%blacs2mpi(send_prow,send_pcol)
        map_send_size(proc_send)=map_send_size(proc_send)+1
      END DO
    END DO

    ! 2) loop over my local data of fm_mat_S and define a map for the proc from which rec data
    ALLOCATE(map_rec_size(0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    map_rec_size=0
    DO jjB=1, ncol_local_d
      j_global=col_indices_d(jjB)
      rec_pcol=cp_fm_indxg2p(j_global,ncol_block_s,dummy_proc,&
                             fm_mat_source%matrix_struct%first_p_pos(2),npcol_s)
      DO iiB=1, nrow_local_d
        i_global=row_indices_d(iiB)
        rec_prow=cp_fm_indxg2p(i_global,nrow_block_s,dummy_proc,&
                               fm_mat_source%matrix_struct%first_p_pos(1),nprow_s)
        proc_receive=fm_mat_source%matrix_struct%context%blacs2mpi(rec_prow,rec_pcol)
        map_rec_size(proc_receive)=map_rec_size(proc_receive)+1
      END DO
    END DO

    ! 3) calculate the number of messages to send and allocate the send buffer
    number_of_send=0
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        number_of_send=number_of_send+1
      END IF
    END DO

    ALLOCATE(buffer_send(number_of_send),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! 4) calculate the number of messages to receive and allocate the rec buffer
    number_of_rec=0
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)
      IF(map_rec_size(proc_receive)>0) THEN
        number_of_rec=number_of_rec+1
      END IF
    END DO

    ALLOCATE(buffer_rec(number_of_rec),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE initialize_buffer

  SUBROUTINE fm_redistribute(fm_mat_source,fm_mat_dest,RPA_proc_map,buffer_rec,buffer_send,&
                             number_of_rec,number_of_send,&
                             map_send_size,map_rec_size,local_size_source,para_env_RPA,error)
    TYPE(cp_fm_type), POINTER                :: fm_mat_source, fm_mat_dest
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: RPA_proc_map
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send
    INTEGER                                  :: number_of_rec, number_of_send
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map_send_size, map_rec_size
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_size_source
    TYPE(cp_para_env_type), POINTER          :: para_env_RPA
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_redistribute', &
      routineP = moduleN//':'//routineN

    INTEGER :: dummy_proc, handle, handle2, i_global, i_local, iiB, iii, &
      j_global, j_local, jjB, mypcol_d, mypcol_s, myprow_d, myprow_s, &
      ncol_block, ncol_block_d, ncol_block_s, ncol_local, ncol_local_d, &
      ncol_local_s, npcol_d, npcol_s, nprow_d, nprow_s, nrow_block, &
      nrow_block_d, nrow_block_s, nrow_local, nrow_local_d, nrow_local_s, &
      proc_receive, proc_send, proc_shift, rec_counter, rec_local_col, &
      rec_local_row, rec_pcol, rec_prow, ref_rec_pcol, ref_rec_prow, &
      send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer, &
      stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index_counter, proc2counter, &
                                                req_send
    INTEGER, DIMENSION(:), POINTER           :: col_indices_d, col_indices_s, &
                                                row_indices_d, row_indices_s
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get info for the source
    CALL cp_fm_get_info(matrix=fm_mat_source,&
                        nrow_local=nrow_local_s,&
                        ncol_local=ncol_local_s,&
                        row_indices=row_indices_s,&
                        col_indices=col_indices_s,&
                        nrow_block=nrow_block_s,&
                        ncol_block=ncol_block_s,&
                        error=error)
    myprow_s=fm_mat_source%matrix_struct%context%mepos(1)
    mypcol_s=fm_mat_source%matrix_struct%context%mepos(2)
    nprow_s =fm_mat_source%matrix_struct%context%num_pe(1)
    npcol_s =fm_mat_source%matrix_struct%context%num_pe(2)

    ! get info for the dest
    CALL cp_fm_get_info(matrix=fm_mat_dest,&
                        nrow_local=nrow_local_d,&
                        ncol_local=ncol_local_d,&
                        row_indices=row_indices_d,&
                        col_indices=col_indices_d,&
                        nrow_block=nrow_block_d,&
                        ncol_block=ncol_block_d,&
                        error=error)
    myprow_d=fm_mat_dest%matrix_struct%context%mepos(1)
    mypcol_d=fm_mat_dest%matrix_struct%context%mepos(2)
    nprow_d =fm_mat_dest%matrix_struct%context%num_pe(1)
    npcol_d =fm_mat_dest%matrix_struct%context%num_pe(2)

    ! 0) check if the local data has to be stored in the new fm_mat_S
    CALL timeset(routineN//"_loc",handle2)
    IF(map_rec_size(para_env_RPA%mepos)>0) THEN
      DO jjB=1, ncol_local_s
        j_global=col_indices_s(jjB)
        rec_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
                               fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)   
        IF(rec_pcol==mypcol_d) THEN
          j_local=cp_fm_indxg2l(j_global,ncol_block_d,dummy_proc,&
                                fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
          DO iiB=1, nrow_local_s
            i_global=row_indices_s(iiB)
            rec_prow=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
                                   fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
            IF(rec_prow==myprow_d) THEN
              i_local=cp_fm_indxg2l(i_global,nrow_block_d,dummy_proc,&
                                    fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
              fm_mat_dest%local_data(i_local,j_local)=fm_mat_source%local_data(iiB,jjB)
            END IF
          END DO
        END IF
      END DO
    END IF
    CALL timestop(handle2)

    ! 1) prepare sending buffer
    CALL timeset(routineN//"_bS",handle2)
    ! new version
    ALLOCATE(proc2counter(0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    proc2counter=0
    ! allocate buffer for sending
    send_counter=0
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
      size_send_buffer=map_send_size(proc_send)

      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        ! prepare the sending buffer
        ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        buffer_send(send_counter)%msg=0.0_dp
        buffer_send(send_counter)%proc=proc_send
    
        proc2counter(proc_send)=send_counter

      END IF
    END DO

    ! fill the buffer
    ALLOCATE(index_counter(0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    index_counter=0
    DO iiB=1, nrow_local_s
      i_global=row_indices_s(iiB)
      send_prow=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
                              fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
      DO jjB=1, ncol_local_s
        j_global=col_indices_s(jjB)
        send_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
                                fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
        iii=fm_mat_dest%matrix_struct%context%blacs2mpi(send_prow,send_pcol)
        IF(iii==para_env_RPA%mepos) CYCLE
        index_counter(iii)=index_counter(iii)+1
        send_counter=proc2counter(iii)
        buffer_send(send_counter)%msg(index_counter(iii))=fm_mat_source%local_data(iiB,jjB)
      END DO
    END DO

    DEALLOCATE(proc2counter)
    DEALLOCATE(index_counter)

    ! ! old version
    ! send_counter=0
    ! DO proc_shift=1, para_env_RPA%num_pe-1
    !   proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
    !   size_send_buffer=map_send_size(proc_send)

    !   IF(map_send_size(proc_send)>0) THEN
    !     send_counter=send_counter+1
    !     ! prepare the sending buffer
    !     ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
    !     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !     buffer_send(send_counter)%msg=0.0_dp
    !     buffer_send(send_counter)%proc=proc_send

    !     ref_send_prow=fm_mat_dest%matrix_struct%context%mpi2blacs(1,proc_send)
    !     ref_send_pcol=fm_mat_dest%matrix_struct%context%mpi2blacs(2,proc_send)

    !     ! iii=0
    !     ! DO jjB=1, ncol_local_s
    !     !   j_global=col_indices_s(jjB)
    !     !   send_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
    !     !                           fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
    !     !   IF(ref_send_pcol==send_pcol) THEN
    !     !     DO iiB=1, nrow_local_s
    !     !       i_global=row_indices_s(iiB)
    !     !       send_pcol=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
    !     !                               fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
    !     !       IF(ref_send_prow==send_pcol) THEN
    !     iii=0
    !     DO iiB=1, nrow_local_s
    !       i_global=row_indices_s(iiB)
    !       send_pcol=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
    !                               fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
    !       IF(ref_send_prow==send_pcol) THEN
    !         DO jjB=1, ncol_local_s
    !           j_global=col_indices_s(jjB)
    !           send_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
    !                                   fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
    !           IF(ref_send_pcol==send_pcol) THEN
    !             iii=iii+1
    !             buffer_send(send_counter)%msg(iii)=fm_mat_source%local_data(iiB,jjB)
    !           END IF
    !         END DO
    !       END IF
    !     END DO

    !   END IF
    ! END DO
    CALL timestop(handle2)

    ! 2) prepare receiving buffer
    CALL timeset(routineN//"_post",handle2)
    rec_counter=0
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)

      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1
        ! prepare the buffer for receive
        ALLOCATE(buffer_rec(rec_counter)%msg(size_rec_buffer),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        buffer_rec(rec_counter)%msg=0.0_dp
        buffer_rec(rec_counter)%proc=proc_receive

        ! post the receiving message
        CALL mp_irecv(buffer_rec(rec_counter)%msg,proc_receive,para_env_RPA%group,buffer_rec(rec_counter)%msg_req)
      END IF
    END DO

    ! 3) send data
    ALLOCATE(req_send(number_of_send),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    send_counter=0
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        CALL mp_isend(buffer_send(send_counter)%msg,proc_send,para_env_RPA%group,buffer_send(send_counter)%msg_req)
        req_send(send_counter)=buffer_send(send_counter)%msg_req
      END IF
    END DO
    CALL timestop(handle2)

    ! 4) fill the fm_mat_dest matrix with the received data
    CALL timeset(routineN//"_fill",handle2)
    rec_counter=0
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)

      rec_local_row=local_size_source(1,proc_receive)
      rec_local_col=local_size_source(2,proc_receive)

      ref_rec_prow=fm_mat_source%matrix_struct%context%mpi2blacs(1,proc_receive)
      ref_rec_pcol=fm_mat_source%matrix_struct%context%mpi2blacs(2,proc_receive)

      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1

        CALL mp_wait(buffer_rec(rec_counter)%msg_req)

        ! iii=0
        ! DO jjB=1, rec_local_col
        !   j_global=cp_fm_indxl2g(jjB,ncol_block_s,ref_rec_pcol,&
        !                          fm_mat_source%matrix_struct%first_p_pos(2),npcol_s)
        !   rec_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
        !                          fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
        !   IF(rec_pcol==mypcol_d) THEN
        !     j_local=cp_fm_indxg2l(j_global,ncol_block_d,dummy_proc,&
        !                           fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
        !     DO iiB=1, rec_local_row
        !       i_global=cp_fm_indxl2g(iiB,nrow_block_s,ref_rec_prow,&
        !                              fm_mat_source%matrix_struct%first_p_pos(1),nprow_s) 
        !       rec_prow=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
        !                              fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
        !       IF(rec_prow==myprow_d) THEN
        !         i_local=cp_fm_indxg2l(i_global,nrow_block_d,dummy_proc,&
        !                               fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
        iii=0
        DO iiB=1, rec_local_row
          i_global=cp_fm_indxl2g(iiB,nrow_block_s,ref_rec_prow,&
                                 fm_mat_source%matrix_struct%first_p_pos(1),nprow_s)
          rec_prow=cp_fm_indxg2p(i_global,nrow_block_d,dummy_proc,&
                                 fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)
          IF(rec_prow==myprow_d) THEN
            i_local=cp_fm_indxg2l(i_global,nrow_block_d,dummy_proc,&
                                  fm_mat_dest%matrix_struct%first_p_pos(1),nprow_d)  
            DO jjB=1, rec_local_col
              j_global=cp_fm_indxl2g(jjB,ncol_block_s,ref_rec_pcol,&
                                     fm_mat_source%matrix_struct%first_p_pos(2),npcol_s)
              rec_pcol=cp_fm_indxg2p(j_global,ncol_block_d,dummy_proc,&
                                     fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)
              IF(rec_pcol==mypcol_d) THEN
                j_local=cp_fm_indxg2l(j_global,ncol_block_d,dummy_proc,&
                                      fm_mat_dest%matrix_struct%first_p_pos(2),npcol_d)

                iii=iii+1
                fm_mat_dest%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
              END IF
            END DO
          END IF
        END DO

        DEALLOCATE(buffer_rec(rec_counter)%msg)

      END IF
    END DO
    CALL timestop(handle2)

    ! 5) wait all
    CALL timeset(routineN//"_wA",handle2)
    CALL mp_waitall(req_send(:))
    DEALLOCATE(req_send)
    DO send_counter=1, number_of_send
      DEALLOCATE(buffer_send(send_counter)%msg)
    END DO
    CALL timestop(handle2)

    CALL timestop(handle)

  END SUBROUTINE fm_redistribute

  SUBROUTINE release_buffer(RPA_proc_map,buffer_rec,buffer_send,&
                            number_of_rec,number_of_send,&
                            map_send_size,map_rec_size,local_size_source,error)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: RPA_proc_map
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send
    INTEGER                                  :: number_of_rec, number_of_send
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map_send_size, map_rec_size
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_size_source
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_buffer', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, rec_counter, &
                                                send_counter
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    IF(ALLOCATED(RPA_proc_map))   DEALLOCATE(RPA_proc_map)
    IF(ALLOCATED(map_send_size))  DEALLOCATE(map_send_size)
    IF(ALLOCATED(map_rec_size))   DEALLOCATE(map_rec_size)
    IF(ALLOCATED(local_size_source)) DEALLOCATE(local_size_source)

    IF(ALLOCATED(buffer_send)) THEN
      DO send_counter=1, number_of_send
        IF(ASSOCIATED(buffer_send(send_counter)%msg)) DEALLOCATE(buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE(buffer_send)
    END IF

    IF(ALLOCATED(buffer_rec)) THEN
      DO rec_counter=1, number_of_rec
        IF(ASSOCIATED(buffer_rec(rec_counter)%msg)) DEALLOCATE(buffer_rec(rec_counter)%msg)
      END DO
      DEALLOCATE(buffer_rec)
    END IF

    CALL timestop(handle)

  END SUBROUTINE release_buffer

END MODULE rpa_communication

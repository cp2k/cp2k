!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sm_fm_interactions [1.0] *
!!
!!   NAME
!!     cp_sm_fm_interactions
!!
!!   FUNCTION
!!     used for all operations involving cp_fm_types and sparse matrices (real_matrix) 
!!
!!   NOTES
!!     first version : most routines imported 
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
MODULE cp_sm_fm_interactions
  
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size


  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type


  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_generate_filename,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
   
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_type

  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             sparse_plus_loc_loct,&
                                             sparse_times_local

  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_gather,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync, &
                                             mp_alltoall, &
                                             mp_shift



  USE timings,                         ONLY: timeset,&
                                             timestop

  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: copy_fm_to_sm, &
            copy_sm_to_fm, &
            cp_sm_fm_multiply, &
            cp_sm_plus_fm_fm_t

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_sm_fm_interactions'

!****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE copy_fm_to_sm(fm,real_matrix)

!   Purpose: Copy a BLACS matrix to a sparse matrix. 
!   History: - Creation (06.06.2001, Matthias Krack)
!   08.2002 adapted to local_data, could be optimized
!   08.2003 imported form qs_blacs
! ***************************************************************************

    TYPE(cp_fm_type), POINTER          :: fm
    TYPE(real_matrix_type), POINTER    :: real_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_fm_to_sm (MODULE cp_sm_fm_interactions)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: fm_block,real_matrix_block
    TYPE(cp_blacs_env_type), POINTER  :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_fm_to_sm","I","",handle)

    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context

    CALL get_matrix_info(matrix=real_matrix,&
                         nblock_row=nblock_row,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block

#if defined(__SCALAPACK)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context%group,iprow,ipcol)

        IF (ipe /= mype) THEN
          nrow_local = fm%matrix_struct%nrow_locals(iprow)
          ncol_local = fm%matrix_struct%ncol_locals(ipcol)
          ALLOCATE (fm_block(nrow_local,ncol_local), STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"fm_block",nrow_local*ncol_local)
        ELSE
          fm_block => fm%local_data
        END IF


        CALL mp_bcast(fm_block,ipe,group)

        DO iblock_row=1,nblock_row

          block_node => first_block_node(matrix=real_matrix,&
                                         block_row=iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                BLOCK=real_matrix_block)

            icol = 1

            DO icol_global=first_col(iblock_col),last_col(iblock_col)

              jpcol = indxg2p(icol_global,ncol_block,mypcol,&
                   fm%matrix_struct%first_p_pos(2),npcol)

              IF (jpcol == ipcol) THEN

                icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                     fm%matrix_struct%first_p_pos(2),npcol)

                irow = 1

                DO irow_global=first_row(iblock_row),last_row(iblock_row)

                  jprow = indxg2p(irow_global,nrow_block,myprow,&
                       fm%matrix_struct%first_p_pos(1),nprow)

                  IF (jprow == iprow) THEN

                     irow_local = indxg2l(irow_global,nrow_block,myprow,&
                          fm%matrix_struct%first_p_pos(1),nprow)

                    real_matrix_block(irow,icol) = fm_block(irow_local,&
                                                          icol_local)

                  END IF

                  irow = irow + 1

                END DO

              END IF

              icol = icol + 1

            END DO

            block_node => next_block_node(block_node)

          END DO

        END DO

        IF (ipe /= mype) THEN
          DEALLOCATE (fm_block,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"fm_block")
        END IF

      END DO
    END DO

#else

    fm_block => fm%local_data

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=real_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=real_matrix_block)

        icol = 1

        DO icol_global=first_col(iblock_col),last_col(iblock_col)

          irow = 1

          DO irow_global=first_row(iblock_row),last_row(iblock_row)

            real_matrix_block(irow,icol) = fm_block(irow_global,icol_global)

            irow = irow + 1

          END DO

          icol = icol + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_fm_to_sm

#if 1
! *****************************************************************************
  SUBROUTINE copy_sm_to_fm(real_matrix,fm,error)

!   Purpose: Copy a real_matrix to a fm. It assumes that a block is at most
!   present once in the real_matrix. 
!   Rewrite of the original routines for improved efficiency (Joost VandeVondele, 2003-08)
!   imported from qs_blacs

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER    :: real_matrix ! the matrix to copy FROM
    TYPE(cp_fm_type), POINTER :: fm ! the matrix to copy TO
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_sm_to_fm (MODULE cp_sm_fm_interactions)",&
      routineN="copy_sm_to_fm"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,nblock_col,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local

    INTEGER :: ncol_small_block, nrow_small_block,nrow_global,ncol_global, &
               iblock,jblock,boundary_sparse,boundary_full,isblock

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row, last_col, last_row
    REAL(wp), DIMENSION(:,:), POINTER :: fm_block, real_matrix_block
    TYPE(cp_blacs_env_type), POINTER  :: context

    INTEGER :: i,j,k
    INTEGER, DIMENSION(:), POINTER   ::  number_of_blocks_row, cum_num_of_blocks_row
    INTEGER, DIMENSION(:,:), POINTER :: block_info_row
    INTEGER                          :: total_number_of_blocks_row 
    INTEGER, DIMENSION(:), POINTER   ::  number_of_blocks_col, cum_num_of_blocks_col
    INTEGER, DIMENSION(:,:), POINTER :: block_info_col
    INTEGER                          :: total_number_of_blocks_col 
    INTEGER, DIMENSION(:,:), POINTER :: number_of_blocks
    INTEGER, DIMENSION(:,:), POINTER :: total_size_of_blocks
    INTEGER :: total_send_blocks,total_recv_blocks
    INTEGER :: total_send_data,total_recv_data
    INTEGER, DIMENSION(:), POINTER :: block_send,block_recv
    REAL(wp), DIMENSION(:), POINTER :: data_send, data_recv
    INTEGER, DIMENSION(:), POINTER :: block_offset
    INTEGER, DIMENSION(:), POINTER :: data_offset
    INTEGER :: sb_offset_i,sb_offset_j,smblock_row,smblock_col
    INTEGER, DIMENSION(:), POINTER :: send_offset,send_count
    INTEGER, DIMENSION(:), POINTER :: recv_offset,recv_count
    LOGICAL :: failure

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p
#endif
!   ---------------------------------------------------------------------------
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(real_matrix),cp_failure_level,routine,error,failure)
    IF (failure) RETURN
    

    CALL timeset("copy_sm_to_fm","I","",handle)
    NULLIFY(number_of_blocks_row, cum_num_of_blocks_row,block_info_row,number_of_blocks_col,cum_num_of_blocks_col, &
            block_info_col,number_of_blocks,total_size_of_blocks,block_send,block_recv,data_send,data_recv,block_offset, &
            data_offset,send_offset,send_count,recv_offset,recv_count)


    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    ! info about the real matrix
    CALL get_matrix_info(matrix=real_matrix,&
                         nblock_row=nblock_row,&
                         nblock_col=nblock_col,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)


#if defined(__SCALAPACK)
    ! find the blocks that are commensurate with both sparse_matrix and cp_fm blocks 
    ! (i.e. they fully live on just one cpu in both the layouts)

    ! rows
    ! first the counting
    ALLOCATE(number_of_blocks_row(nblock_row),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"number_of_blocks_row",nblock_row)
    ALLOCATE(cum_num_of_blocks_row(nblock_row),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cum_num_of_blocks_row",nblock_row)
    number_of_blocks_row=0
    iblock=1
    boundary_sparse=last_row(iblock)
    boundary_full=nrow_block
    DO
       IF (boundary_full .eq. boundary_sparse) boundary_full=boundary_full+nrow_block
       number_of_blocks_row(iblock)=number_of_blocks_row(iblock)+1
       IF (boundary_sparse .lt. boundary_full) THEN
          iblock=iblock+1
       ELSE
          boundary_full=boundary_full+nrow_block
       ENDIF
       if (iblock.gt.nblock_row) EXIT
       boundary_sparse=last_row(iblock)
    ENDDO
    total_number_of_blocks_row=SUM(number_of_blocks_row)
    cum_num_of_blocks_row(1)=0
    DO iblock=2,nblock_row
       cum_num_of_blocks_row(iblock)=cum_num_of_blocks_row(iblock-1)+number_of_blocks_row(iblock-1) 
    ENDDO

    ! collect properties
    ALLOCATE(block_info_row(3,total_number_of_blocks_row),STAT=istat) !start,end,blacs_pe 
    IF (istat /= 0) CALL stop_memory(routine,"block_info_row",3*total_number_of_blocks_row)
    iblock=1
    isblock=0
    boundary_sparse=last_row(iblock)
    boundary_full=nrow_block
    DO
       isblock=isblock+1
       IF (boundary_full .eq. boundary_sparse) boundary_full=boundary_full+nrow_block
       IF (boundary_sparse .lt. boundary_full) THEN
          block_info_row(2,isblock)=boundary_sparse
          iblock=iblock+1
       ELSE
          block_info_row(2,isblock)=boundary_full
          boundary_full=boundary_full+nrow_block
       ENDIF
       block_info_row(3,isblock)= indxg2p(block_info_row(2,isblock),nrow_block,myprow,&
                                          fm%matrix_struct%first_p_pos(1),nprow)
       if (iblock.gt.nblock_row) EXIT
       boundary_sparse=last_row(iblock)
    ENDDO
    block_info_row(1,1)=1
    DO iblock=2,total_number_of_blocks_row
       block_info_row(1,iblock)=block_info_row(2,iblock-1)+1 
    ENDDO

    ! cols
    ! first the counting
    ALLOCATE(number_of_blocks_col(nblock_col),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"number_of_blocks_col",nblock_col)
    ALLOCATE(cum_num_of_blocks_col(nblock_col),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cum_num_of_blocks_col",nblock_col)
    number_of_blocks_col=0
    iblock=1
    boundary_sparse=last_col(iblock)
    boundary_full=ncol_block
    DO
       IF (boundary_full .eq. boundary_sparse) boundary_full=boundary_full+ncol_block
       number_of_blocks_col(iblock)=number_of_blocks_col(iblock)+1
       IF (boundary_sparse .lt. boundary_full) THEN
          iblock=iblock+1
       ELSE
          boundary_full=boundary_full+ncol_block
       ENDIF
       if (iblock.gt.nblock_col) EXIT
       boundary_sparse=last_col(iblock)
    ENDDO
    total_number_of_blocks_col=SUM(number_of_blocks_col)
    cum_num_of_blocks_col(1)=0
    DO iblock=2,nblock_col
       cum_num_of_blocks_col(iblock)=cum_num_of_blocks_col(iblock-1)+number_of_blocks_col(iblock-1) 
    ENDDO

    ! collect properties
    ALLOCATE(block_info_col(3,total_number_of_blocks_col),STAT=istat) !start,end,blacs_pe 
    IF (istat /= 0) CALL stop_memory(routine,"block_info_col",3*total_number_of_blocks_col)
    iblock=1
    isblock=0
    boundary_sparse=last_col(iblock)
    boundary_full=ncol_block
    DO
       isblock=isblock+1
       IF (boundary_full .eq. boundary_sparse) boundary_full=boundary_full+ncol_block
       IF (boundary_sparse .lt. boundary_full) THEN
          block_info_col(2,isblock)=boundary_sparse
          iblock=iblock+1
       ELSE
          block_info_col(2,isblock)=boundary_full
          boundary_full=boundary_full+ncol_block
       ENDIF
       block_info_col(3,isblock)= indxg2p(block_info_col(2,isblock),ncol_block,mypcol,&
                                          fm%matrix_struct%first_p_pos(2),npcol)
       if (iblock.gt.nblock_col) EXIT
       boundary_sparse=last_col(iblock)
    ENDDO
    block_info_col(1,1)=1
    DO iblock=2,total_number_of_blocks_col
       block_info_col(1,iblock)=block_info_col(2,iblock-1)+1 
    ENDDO

! now count how much / what data has to go to what cpu
    ALLOCATE(number_of_blocks(0:npe-1,0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"number_of_blocks",npe*npe)
    ALLOCATE(total_size_of_blocks(0:npe-1,0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"total_size_of_blocks",npe*npe)

    number_of_blocks=0
    total_size_of_blocks=0
    DO iblock_row=1,nblock_row
       block_node => first_block_node(matrix=real_matrix, block_row=iblock_row)
       DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,block_col=iblock_col)
          DO jblock=1,number_of_blocks_col(iblock_col)
             DO iblock=1,number_of_blocks_row(iblock_row)
                smblock_row=iblock+cum_num_of_blocks_row(iblock_row)
                smblock_col=jblock+cum_num_of_blocks_col(iblock_col)
                ! ipe is the target processor for this block
                ipe = blacs_pnum(context%group,block_info_row(3,smblock_row), &
                                               block_info_col(3,smblock_col))
                number_of_blocks(ipe,mype)=number_of_blocks(ipe,mype)+1
                total_size_of_blocks(ipe,mype)=total_size_of_blocks(ipe,mype)+ &
                         (block_info_row(2,smblock_row)-block_info_row(1,smblock_row)+1)* &
                         (block_info_col(2,smblock_col)-block_info_col(1,smblock_col)+1)
             ENDDO
          ENDDO
          block_node => next_block_node(block_node)
       END DO
    END DO

! get all info everywhere using the usual ugly sum trick
    CALL mp_sum(number_of_blocks,group)
    CALL mp_sum(total_size_of_blocks,group)

! packing data
    total_send_blocks = SUM(number_of_blocks(:,mype))
    total_recv_blocks = SUM(number_of_blocks(mype,:))
    total_send_data = SUM(total_size_of_blocks(:,mype))
    total_recv_data = SUM(total_size_of_blocks(mype,:))


    ALLOCATE(block_send(2*total_send_blocks),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"block_send",2*total_send_blocks)
    ALLOCATE(block_recv(2*total_recv_blocks),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"block_send",2*total_recv_blocks)
    ALLOCATE(data_send(total_send_data),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"data_send",total_send_data)
    ALLOCATE(data_recv(total_recv_data),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"data_recv",total_recv_data)

    ALLOCATE(block_offset(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"block_offset",npe)
    ALLOCATE(data_offset(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"data_offset",npe)
    block_offset(0)=0
    data_offset(0)=0
    DO ipe=1,npe-1
       block_offset(ipe)=block_offset(ipe-1)+2*number_of_blocks(ipe-1,mype)
       data_offset(ipe) =data_offset(ipe-1) +total_size_of_blocks(ipe-1,mype)
    ENDDO
    data_send=-1
    data_recv=-2
    DO iblock_row=1,nblock_row
       block_node => first_block_node(matrix=real_matrix, block_row=iblock_row)
       DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,block_col=iblock_col,block=real_matrix_block)
          DO jblock=1,number_of_blocks_col(iblock_col)
             DO iblock=1,number_of_blocks_row(iblock_row)
                smblock_row=iblock+cum_num_of_blocks_row(iblock_row)
                smblock_col=jblock+cum_num_of_blocks_col(iblock_col)
                ! ipe is the target processor for this block
                ipe = blacs_pnum(context%group,block_info_row(3,smblock_row), &
                                               block_info_col(3,smblock_col))
                block_send(block_offset(ipe)+1)=smblock_row
                block_send(block_offset(ipe)+2)=smblock_col
                block_offset(ipe)=block_offset(ipe)+2
                sb_offset_i = block_info_row(1,smblock_row) - first_row(iblock_row)
                sb_offset_j = block_info_col(1,smblock_col) - first_col(iblock_col)
                DO j=1,block_info_col(2,smblock_col)-block_info_col(1,smblock_col)+1
                   DO i=1,block_info_row(2,smblock_row)-block_info_row(1,smblock_row)+1
                      data_offset(ipe)=data_offset(ipe)+1 
                      data_send(data_offset(ipe))=real_matrix_block(sb_offset_i+i,sb_offset_j+j)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          block_node => next_block_node(block_node)
       END DO
    END DO

! send around the data
    ALLOCATE(send_offset(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"send_offset",npe)
    ALLOCATE(recv_offset(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"recv_offset",npe)
    ALLOCATE(send_count(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"send_count",npe)
    ALLOCATE(recv_count(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"recv_count",npe)

    send_count(:)=2*number_of_blocks(:,mype)
    recv_count(:)=2*number_of_blocks(mype,:)
    send_offset(0)=0
    recv_offset(0)=0
    DO ipe=1,npe-1
       send_offset(ipe)=send_offset(ipe-1)+send_count(ipe-1)
       recv_offset(ipe)=recv_offset(ipe-1)+recv_count(ipe-1)
    ENDDO
    CALL mp_alltoall(block_send,send_count,send_offset,block_recv,recv_count,recv_offset,group)

    send_count(:)=total_size_of_blocks(:,mype)
    recv_count(:)=total_size_of_blocks(mype,:)
    send_offset(0)=0
    recv_offset(0)=0
    DO ipe=1,npe-1
       send_offset(ipe)=send_offset(ipe-1)+send_count(ipe-1)
       recv_offset(ipe)=recv_offset(ipe-1)+recv_count(ipe-1)
    ENDDO
    ! write(6,*) "IN",mype,data_send
    CALL mp_alltoall(data_send,send_count,send_offset,data_recv,recv_count,recv_offset,group)
    ! write(6,*) "OUT",mype,data_recv

! Finally, unpack the data filling in the fm%local_data
    fm%local_data=0.0_wp
    k=0
    DO iblock=1,total_recv_blocks
       smblock_row=block_recv((iblock-1)*2+1)
       smblock_col=block_recv((iblock-1)*2+2)
       sb_offset_i = indxg2l(block_info_row(1,smblock_row),nrow_block,myprow,&
                     fm%matrix_struct%first_p_pos(1),nprow)-1
       sb_offset_j = indxg2l(block_info_col(1,smblock_col),ncol_block,mypcol,&
                     fm%matrix_struct%first_p_pos(2),npcol)-1
       DO j=1,block_info_col(2,smblock_col)-block_info_col(1,smblock_col)+1
          DO i=1,block_info_row(2,smblock_row)-block_info_row(1,smblock_row)+1
             k=k+1
             fm%local_data(sb_offset_i+i,sb_offset_j+j)=data_recv(k)
          ENDDO
       ENDDO
    ENDDO

! clean up all allocations
    DEALLOCATE(number_of_blocks_row, cum_num_of_blocks_row,block_info_row,number_of_blocks_col,cum_num_of_blocks_col, &
            block_info_col,number_of_blocks,total_size_of_blocks,block_send,block_recv,data_send,data_recv,block_offset, &
            data_offset,send_offset,send_count,recv_offset,recv_count,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"deallocate",0)

#else

    fm_block => fm%local_data

    fm_block(:,:) = 0.0_wp

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=real_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=real_matrix_block)

        icol = 1

        DO icol_global=first_col(iblock_col),last_col(iblock_col)

          irow = 1

          DO irow_global=first_row(iblock_row),last_row(iblock_row)

            fm_block(irow_global,icol_global) = real_matrix_block(irow,icol)

            irow = irow + 1

          END DO

          icol = icol + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_sm_to_fm

#endif


!
!
! keep to old routine alive as to allow for easy testing. Should however be
! removed quite soon. 
!
!
#if 0
! *****************************************************************************
  SUBROUTINE copy_sm_to_fm(real_matrix,fm)

!   Purpose: Copy a real_matrix to a fm. The fm blocks
!            are allocated during the copy procedure.

!   temporary old version for reference only (2003-08-02)
!   imported from qs_blacs
!   ***************************************************************************

    TYPE(real_matrix_type), POINTER    :: real_matrix ! the matrix to copy FROM
    TYPE(cp_fm_type), POINTER :: fm ! the matrix to copy TO

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_sm_to_fm (MODULE cp_sm_fm_interactions)",&
      routineN="copy_sm_to_fm"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row, last_col, last_row
    REAL(wp), DIMENSION(:,:), POINTER :: fm_block, real_matrix_block
    TYPE(cp_blacs_env_type), POINTER  :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_sm_to_fm","I","",handle)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block

    ! info about the real matrix
    CALL get_matrix_info(matrix=real_matrix,&
                         nblock_row=nblock_row,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)
#if defined(__SCALAPACK)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context%group,iprow,ipcol)

        IF (ipe /= mype) THEN
          nrow_local = fm%matrix_struct%nrow_locals(iprow)
          ncol_local = fm%matrix_struct%ncol_locals(ipcol)
          ALLOCATE (fm_block(nrow_local,ncol_local),&
                    STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"alloc fm block",&
                             nrow_local*ncol_local*wp_size)
          END IF
        ELSE
           fm_block => fm%local_data
        END IF

        fm_block(:,:) = 0.0_wp

        DO iblock_row=1,nblock_row

          block_node => first_block_node(matrix=real_matrix,&
                                         block_row=iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                BLOCK=real_matrix_block)

            icol = 1

            DO icol_global=first_col(iblock_col),last_col(iblock_col)

               jpcol = indxg2p(icol_global,ncol_block,mypcol,&
                    fm%matrix_struct%first_p_pos(2),npcol)

              IF (jpcol == ipcol) THEN

                icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                     fm%matrix_struct%first_p_pos(2),npcol)

                irow = 1

                DO irow_global=first_row(iblock_row),last_row(iblock_row)

                  jprow = indxg2p(irow_global,nrow_block,myprow,&
                       fm%matrix_struct%first_p_pos(1),nprow)

                  IF (jprow == iprow) THEN

                    irow_local = indxg2l(irow_global,nrow_block,myprow,&
                         fm%matrix_struct%first_p_pos(1),nprow)

                    fm_block(irow_local,icol_local) = real_matrix_block(irow,&
                                                                      icol)

                  END IF

                  irow = irow + 1

                END DO

              END IF

              icol = icol + 1

            END DO

            block_node => next_block_node(block_node)

          END DO

        END DO

        CALL mp_sum(fm_block,ipe,group)

        IF (ipe /= mype) THEN
          DEALLOCATE (fm_block,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"dealloc fm block")
          END IF
        END IF

      END DO
    END DO

#else


    fm_block => fm%local_data

    fm_block(:,:) = 0.0_wp

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=real_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=real_matrix_block)

        icol = 1

        DO icol_global=first_col(iblock_col),last_col(iblock_col)

          irow = 1

          DO irow_global=first_row(iblock_row),last_row(iblock_row)

            fm_block(irow_global,icol_global) = real_matrix_block(irow,icol)

            irow = irow + 1

          END DO

          icol = icol + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_sm_to_fm

#endif

!******************************************************************************
!
! new version of the cp_sm_fm_multiply, optimized communication, 
! imported from cp_fm_basic_linalg
!
!*******************************************************************************
  SUBROUTINE cp_sm_fm_multiply(sparse_matrix,v_in,v_out,ncol,para_env, &
       alpha,beta,error)
! to do: better interface
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_fm_type) , POINTER :: v_in
    TYPE(cp_fm_type) , POINTER :: v_out
    INTEGER, INTENT(IN)               :: ncol
    TYPE(cp_para_env_type), POINTER :: para_env ! of the sparse matrix
    REAL(kind=wp), INTENT(in), OPTIONAL :: alpha, beta
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine = "cp_sm_fm_multiply"

!   *** Local variables ***
    LOGICAL :: failure
    REAL(wp), DIMENSION(:,:), POINTER :: local_v_in, local_v_out
    INTEGER :: i,j,info,nrow_global,handle,num_pe,mepos,ncol_max,stat,ishift
    INTEGER, DIMENSION(:), POINTER :: ncol_loc
    INTEGER :: ictxt_loc, im,in, myprow,mypcol,nprow,npcol,desc(9),rb,cb
    REAL(wp), POINTER, DIMENSION(:,:) :: newdat
    REAL(wp) :: alpha_l,beta_l
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: NUMROC
#endif
    
    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_sm_fm_multiply", &
                                   routineP = moduleN//"/"//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"cp_sm_fm_multiply would be most happy to use scalapack libraries"//&
     CPSourceFileRef,&
     error)
#endif
    CALL timeset("cp_sm_fm_multiply","I","",handle)

    CPPrecondition(ASSOCIATED(sparse_matrix),cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(v_in),cp_failure_level,routine,error,failure)
    CPPrecondition(v_in%ref_count>0,cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(v_out),cp_failure_level,routine,error,failure)
    CPPrecondition(v_out%ref_count>0,cp_failure_level,routine,error,failure)
    CPPrecondition(ncol.le.v_in%matrix_struct%ncol_global,cp_failure_level,routine,error,failure)
    CPPrecondition(ncol.le.v_out%matrix_struct%ncol_global,cp_failure_level,routine,error,failure)

    IF (PRESENT(alpha)) THEN
       alpha_l=alpha
    ELSE
       alpha_l=1.0_wp
    ENDIF
    IF (PRESENT(beta)) THEN
       beta_l=beta
    ELSE
       beta_l=0.0_wp
    ENDIF

    failure=.FALSE.
    NULLIFY(local_v_in,local_v_out,newdat)

    CALL cp_fm_get_info(v_in,nrow_global=nrow_global)
    num_pe=para_env%num_pe
    mepos =para_env%mepos
    rb=nrow_global
    cb=1
    im=nrow_global
    ALLOCATE(ncol_loc(0:num_pe-1))
    ncol_loc=0
#if defined(__SCALAPACK)
    CALL blacs_get(-1,0,ictxt_loc)
    CALL blacs_gridinit(ictxt_loc,'R',1,num_pe)
    CALL blacs_gridinfo(ictxt_loc,nprow,npcol,myprow,mypcol)
    CALL descinit(desc,nrow_global,ncol,rb,cb,0,0,ictxt_loc,nrow_global,info)
    CPPostcondition(info==0,cp_failure_level,routine,error,failure)
    in=numroc(ncol,cb,mypcol,0,npcol)

    ALLOCATE(newdat(im,MAX(1,in)))

    ! do the actual scalapack to cols reordering
    CALL pdgemr2d(nrow_global,ncol,v_in%local_data(1,1),1,1,&
                  v_in%matrix_struct%descriptor, &
                  newdat(1,1),1,1,desc,ictxt_loc)
 
    ! obtain data in transposed structures
    ncol_loc(mepos)=in
    CALL mp_sum(ncol_loc,para_env%group)
    ncol_max=MAXVAL(ncol_loc)

    ALLOCATE(local_v_in(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    ALLOCATE(local_v_out(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    DO j=1,in
       DO i=1,im
              local_v_in(j,i)=newdat(i,j)
       END DO
    END DO
#else
    in=ncol
    ncol_loc(mepos)=in
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(local_v_in(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    ALLOCATE(local_v_out(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    DO j=1,in
       DO i=1,im
              local_v_in(j,i)=v_in%local_data(i,j)
       END DO
    END DO
#endif

    IF (beta_l.EQ.0.0_wp) THEN
       CALL dcopy(SIZE(local_v_out,1)*SIZE(local_v_out,2),&
            0.0_wp,0,local_v_out(1,1),1) ! zero output vector
    ELSE ! this rather expensive and should probably not be used. Anyway, not tested
#if defined(__SCALAPACK)
         CALL pdgemr2d(nrow_global,ncol,v_out%local_data(1,1),1,1,&
                       v_out%matrix_struct%descriptor, &
                       newdat(1,1),1,1,desc,ictxt_loc)
         DO j=1,in
            DO i=1,im
                local_v_out(j,i)=beta_l*newdat(i,j)
            END DO
         END DO
#else
       DO j=1,in
          DO i=1,im
                 local_v_out(j,i)=beta_l*v_out%local_data(i,j)
          END DO
       END DO
#endif
    ENDIF

    ! now we shift around the data in a ring, multiplying the v_in with the local sparse matrix
    ! adding the results to v_out, we need to do num_pe shifts

    DO ishift=0,num_pe-1
       CALL sparse_times_local(sparse_matrix,local_v_in,local_v_out,&
            ncol_loc(MODULO(mepos-ishift,num_pe)), ncol_max, para_env=para_env,alpha=alpha_l)
       CALL mp_shift(local_v_in,para_env%group)
       CALL mp_shift(local_v_out,para_env%group)
    ENDDO

#if defined(__SCALAPACK)
    ! copy the result back into the scalapack data structure
    DO j=1,in
       DO i=1,im
              newdat(i,j)=local_v_out(j,i)
       END DO
    END DO
    ! shuffle around the data
    CALL pdgemr2d(nrow_global,ncol, &
                  newdat(1,1),1,1,desc, &
                  v_out%local_data(1,1),1,1,v_out%matrix_struct%descriptor, &
                  ictxt_loc)
    ! give the grid back
    CALL blacs_gridexit(ictxt_loc)
    DEALLOCATE(newdat)
#else
    DO j=1,in
       DO i=1,im
              v_out%local_data(i,j)=local_v_out(j,i)
       END DO
    END DO
#endif

    DEALLOCATE(local_v_in,local_v_out,ncol_loc)

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_fm_multiply

! *****************************************************************************
! computes sparse=sparse+alpha*v*g^T 
! imported from cp_fm_basic_linalg
! *****************************************************************************
  SUBROUTINE cp_sm_plus_fm_fm_t(sparse_matrix,matrix_v,matrix_g,last_col,&
       first_col,alpha,error)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_fm_type) , POINTER :: matrix_v
    TYPE(cp_fm_type) , POINTER, OPTIONAL :: matrix_g
    INTEGER, INTENT(IN)               :: last_col
    REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
    INTEGER, OPTIONAL, INTENT(IN)     :: first_col
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE cp_sm_plus_fm_fm_t (MODULE cp_fm_basic_linalg)"

!   *** Local variables ***
    REAL(wp), DIMENSION(:,:), POINTER :: loc_v,loc_g,newdat
    REAL(wp) :: prefactor
    INTEGER :: col1,i,j,handle,nrow_global,num_pe,mepos,ictxt_loc,desc(9)
    INTEGER :: nprow,npcol,myprow,mypcol,rb,cb,ncol_max,ishift,in,im,info,handle2,ncol
    INTEGER, DIMENSION(:), POINTER :: ncol_loc
    LOGICAL g_present,FAILURE

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_sm_plus_fm_fm_t", &
                                   routineP = moduleN//"/"//routineN

!   ---------------------------------------------------------------------------
    failure=.FALSE.
#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"cp_sm_plus_fm_fm_t would be most happy to use scalapack libraries"//&
     CPSourceFileRef,&
     error)
#endif
    CALL timeset("cp_sm_plus_fm_fm_t","I","",handle)
    NULLIFY(newdat,ncol_loc)

    CPPrecondition(ASSOCIATED(sparse_matrix),cp_failure_level,routine,error,failure)
    CPPrecondition(ASSOCIATED(matrix_v),cp_failure_level,routine,error,failure)
    CPPrecondition(matrix_v%ref_count>0,cp_failure_level,routine,error,failure)

    IF (PRESENT(first_col)) THEN
       col1 = first_col
    ELSE
       col1 = 1
    END IF

    ncol=last_col-col1+1

    g_present=.FALSE.
    IF (PRESENT(matrix_g)) THEN
       CPPrecondition(ASSOCIATED(matrix_g),cp_failure_level,routine,error,failure)
       CPPrecondition(matrix_g%ref_count>0,cp_failure_level,routine,error,failure)
       g_present=.TRUE. 
    END IF
    ! we start with te easiest case

    CALL cp_fm_get_info(matrix_v,nrow_global=nrow_global)
    num_pe=matrix_v%matrix_struct%para_env%num_pe
    mepos =matrix_v%matrix_struct%para_env%mepos
    rb=nrow_global
    cb=1
    im=nrow_global
    ALLOCATE(ncol_loc(0:num_pe-1))

#if defined(__SCALAPACK)
    CALL blacs_get(-1,0,ictxt_loc)
    CALL blacs_gridinit(ictxt_loc,'R',1,num_pe)
    CALL blacs_gridinfo(ictxt_loc,nprow,npcol,myprow,mypcol)
    CALL descinit(desc,nrow_global,ncol,rb,cb,0,0,ictxt_loc,nrow_global,info)
    CPPostcondition(info==0,cp_failure_level,routine,error,failure)
    in=numroc(ncol,cb,mypcol,0,npcol)
    ALLOCATE(newdat(im,MAX(1,in)))
    ncol_loc=0
    ncol_loc(mepos)=in
    CALL mp_sum(ncol_loc,matrix_v%matrix_struct%para_env%group)
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(loc_v(ncol_max,nrow_global))

    ! do the actual scalapack to cols reordering
    CALL timeset("cp_sm_plus_fm_fm_t_pdgemr2d","I","",handle2)
    CALL pdgemr2d(nrow_global,ncol,matrix_v%local_data(1,1),1,col1, &
                  matrix_v%matrix_struct%descriptor, &
                  newdat(1,1),1,1,desc,ictxt_loc)
    CALL timestop(0.0_wp,handle2)
    DO j=1,in
       DO i=1,im
          loc_v(j,i)=newdat(i,j)
       END DO
    END DO

    IF (g_present) THEN
       ALLOCATE(loc_g(ncol_max,nrow_global))
       ! do the actual scalapack to cols reordering
       CALL timeset("cp_sm_plus_fm_fm_t_pdgemr2d","I","",handle2)
       CALL pdgemr2d(nrow_global,ncol,matrix_g%local_data(1,1),1,col1, &
                     matrix_g%matrix_struct%descriptor, &
                     newdat(1,1),1,1,desc,ictxt_loc)
       CALL timestop(0.0_wp,handle2)
       DO j=1,in
          DO i=1,im
             loc_g(j,i)=newdat(i,j)
          END DO
       END DO
    ELSE
       loc_g=>loc_v
    ENDIF
    DEALLOCATE(newdat)
    ! EXIT THE GRID 
    CALL blacs_gridexit(ictxt_loc)
#else
    in=ncol
    ncol_loc=0
    ncol_loc(mepos)=in
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(loc_v(ncol_max,nrow_global))
    DO j=1,in
       DO i=1,im
          loc_v(j,i)=matrix_v%local_data(i,col1-1+j)
       END DO
    END DO
    IF (g_present) THEN
       ALLOCATE(loc_g(ncol_max,nrow_global))
       DO j=1,in
          DO i=1,im
             loc_g(j,i)=matrix_g%local_data(i,col1-1+j)
          END DO
       END DO
    ELSE
       loc_g=>loc_v
    ENDIF
#endif

    IF (PRESENT(alpha)) THEN
       prefactor = alpha
    ELSE
       prefactor = 1.0_wp
    END IF

    DO ishift=0,num_pe-1
       CALL sparse_plus_loc_loct(sparse_matrix,loc_v,loc_g,&
                                 ncol_loc(MODULO(mepos-ishift,num_pe)), &
                                 ncol_max,alpha=prefactor)
       IF (g_present) THEN
          CALL mp_shift(loc_g,matrix_g%matrix_struct%para_env%group)
          CALL mp_shift(loc_v,matrix_v%matrix_struct%para_env%group)
       ELSE
          CALL mp_shift(loc_v,matrix_v%matrix_struct%para_env%group)
       ENDIF
    ENDDO

    IF (g_present) THEN
       DEALLOCATE(loc_v,loc_g)
    ELSE
       DEALLOCATE(loc_v)
    ENDIF

    DEALLOCATE(ncol_loc)
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE cp_sm_plus_fm_fm_t

END MODULE cp_sm_fm_interactions

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
!!     used for all operations involving cp_fm_types and sparse matrices
!!     (real_matrix) 
!!
!!   NOTES
!!     first version : most routines imported 
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
MODULE cp_sm_fm_interactions
  USE cp_array_i_utils,                ONLY: cp_2d_i_write
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type,&
                                             cp_2d_r_write
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_equivalent,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type,&
                                             cp_fm_struct_write
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_maxval,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_element,&
                                             cp_fm_type,&
                                             cp_fm_write
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_generate_filename,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_release,&
                                             distribution_2d_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: lcm
  USE message_passing,                 ONLY: mp_alltoall,&
                                             mp_bcast,&
                                             mp_gather,&
                                             mp_isendrecv,&
                                             mp_max,&
                                             mp_sendrecv,&
                                             mp_shift,&
                                             mp_sum,&
                                             mp_sync,&
                                             mp_wait,&
                                             mp_waitall
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             cp_sm_maxabsval,&
                                             cp_sm_scale_and_add,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             replicate_matrix,&
                                             sparse_plus_loc_loct,&
                                             sparse_times_local
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: copy_fm_to_sm, &
       copy_sm_to_fm, &
       cp_sm_fm_multiply, &
       cp_sm_plus_fm_fm_t, &
       cp_sm_fm_multiply_tests

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
       IF (boundary_full .EQ. boundary_sparse) boundary_full=boundary_full+nrow_block
       number_of_blocks_row(iblock)=number_of_blocks_row(iblock)+1
       IF (boundary_sparse .LT. boundary_full) THEN
          iblock=iblock+1
       ELSE
          boundary_full=boundary_full+nrow_block
       ENDIF
       IF (iblock.gt.nblock_row) EXIT
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
       IF (boundary_full .EQ. boundary_sparse) boundary_full=boundary_full+nrow_block
       IF (boundary_sparse .LT. boundary_full) THEN
          block_info_row(2,isblock)=boundary_sparse
          iblock=iblock+1
       ELSE
          block_info_row(2,isblock)=boundary_full
          boundary_full=boundary_full+nrow_block
       ENDIF
       block_info_row(3,isblock)= indxg2p(block_info_row(2,isblock),nrow_block,myprow,&
            fm%matrix_struct%first_p_pos(1),nprow)
       IF (iblock.gt.nblock_row) EXIT
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
       IF (boundary_full .EQ. boundary_sparse) boundary_full=boundary_full+ncol_block
       number_of_blocks_col(iblock)=number_of_blocks_col(iblock)+1
       IF (boundary_sparse .LT. boundary_full) THEN
          iblock=iblock+1
       ELSE
          boundary_full=boundary_full+ncol_block
       ENDIF
       IF (iblock.gt.nblock_col) EXIT
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
       IF (boundary_full .EQ. boundary_sparse) boundary_full=boundary_full+ncol_block
       IF (boundary_sparse .LT. boundary_full) THEN
          block_info_col(2,isblock)=boundary_sparse
          iblock=iblock+1
       ELSE
          block_info_col(2,isblock)=boundary_full
          boundary_full=boundary_full+ncol_block
       ENDIF
       block_info_col(3,isblock)= indxg2p(block_info_col(2,isblock),ncol_block,mypcol,&
            fm%matrix_struct%first_p_pos(2),npcol)
       IF (iblock.gt.nblock_col) EXIT
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
          CALL get_block_node(block_node=block_node,block_col=iblock_col,BLOCK=real_matrix_block)
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


  !******************************************************************************
  !
  ! new version of the cp_sm_fm_multiply, optimized communication, 
  ! imported from cp_fm_basic_linalg
  !
  !*******************************************************************************
  SUBROUTINE cp_sm_fm_multiply_general(sparse_matrix,v_in,v_out,ncol,para_env, &
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
    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_sm_fm_multiply_general", &
         routineP = moduleN//"/"//routineN

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


    !   ---------------------------------------------------------------------------
    failure=.FALSE.
#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"cp_sm_fm_multiply would be most happy to use scalapack libraries"//&
     CPSourceFileRef,&
     error)
#endif
    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(sparse_matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(v_in%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_out),cp_failure_level,routineP,error,failure)
    CPPrecondition(v_out%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ncol.le.v_in%matrix_struct%ncol_global,cp_failure_level,routineP,error,failure)
    CPPrecondition(ncol.le.v_out%matrix_struct%ncol_global,cp_failure_level,routineP,error,failure)

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
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
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
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(local_v_out(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(local_v_out(ncol_max,nrow_global),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
            ncol_loc(MODULO(mepos-ishift,num_pe)), ncol_max, alpha=alpha_l)
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
  END SUBROUTINE cp_sm_fm_multiply_general

  ! *****************************************************************************
  ! computes sparse=sparse+alpha*v*g^T 
  ! imported from cp_fm_basic_linalg
  ! *****************************************************************************
  SUBROUTINE cp_sm_plus_fm_fm_t_general(sparse_matrix,matrix_v,matrix_g,ncol,&
       alpha,error)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(cp_fm_type) , POINTER :: matrix_v
    TYPE(cp_fm_type) , POINTER, OPTIONAL :: matrix_g
    INTEGER, INTENT(IN), OPTIONAL :: ncol
    REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    REAL(wp), DIMENSION(:,:), POINTER :: loc_v,loc_g,newdat
    REAL(wp) :: prefactor
    INTEGER :: col1,i,j,handle,nrow_global,num_pe,mepos,ictxt_loc,desc(9), my_ncol
    INTEGER :: nprow,npcol,myprow,mypcol,rb,cb,ncol_max,ishift,in,im,info,handle2
    INTEGER, DIMENSION(:), POINTER :: ncol_loc
    LOGICAL g_present,FAILURE

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp_sm_plus_fm_fm_t_general", &
         routineP = moduleN//"/"//routineN

    !   ---------------------------------------------------------------------------
    failure=.FALSE.
#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"cp_sm_plus_fm_fm_t would be most happy to use scalapack libraries"//&
     CPSourceFileRef,&
     error)
#endif
    CALL timeset(routineN,"I","",handle)
    NULLIFY(newdat,ncol_loc)

    CPPrecondition(ASSOCIATED(sparse_matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_v),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix_v%ref_count>0,cp_failure_level,routineP,error,failure)

    col1 = 1

    g_present=.FALSE.
    IF (PRESENT(matrix_g)) THEN
       CPPrecondition(ASSOCIATED(matrix_g),cp_failure_level,routineP,error,failure)
       CPPrecondition(matrix_g%ref_count>0,cp_failure_level,routineP,error,failure)
       g_present=.TRUE. 
    END IF
    ! we start with te easiest case

    CALL cp_fm_get_info(matrix_v,nrow_global=nrow_global, ncol_global=my_ncol)
    IF (PRESENT(ncol)) my_ncol=ncol
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
    CALL descinit(desc,nrow_global,my_ncol,rb,cb,0,0,ictxt_loc,nrow_global,info)
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    in=numroc(my_ncol,cb,mypcol,0,npcol)
    ALLOCATE(newdat(im,MAX(1,in)))
    ncol_loc=0
    ncol_loc(mepos)=in
    CALL mp_sum(ncol_loc,matrix_v%matrix_struct%para_env%group)
    ncol_max=MAXVAL(ncol_loc)
    ALLOCATE(loc_v(ncol_max,nrow_global))

    ! do the actual scalapack to cols reordering
    CALL timeset("cp_sm_plus_fm_fm_t_pdgemr2d","I","",handle2)
    CALL pdgemr2d(nrow_global,my_ncol,matrix_v%local_data(1,1),1,col1, &
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
       CALL pdgemr2d(nrow_global,my_ncol,matrix_g%local_data(1,1),1,col1, &
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
    in=my_ncol
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
  END SUBROUTINE cp_sm_plus_fm_fm_t_general

!!****f* cp_sm_fm_interactions/cp_sm_fm_multiply [1.0] *
!!
!!   NAME
!!     cp_sm_fm_multiply
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_fm_multiply(sparse_matrix, v_in, v_out, ncol,&
!!         para_env, alpha, beta, error)
!!       Type(real_matrix_type), Pointer:: sparse_matrix
!!       Type(cp_fm_type), Pointer:: v_in
!!       Type(cp_fm_type), Pointer:: v_out
!!       Integer, Intent (IN):: ncol
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Real(Kind=wp), Intent (IN), Optional:: alpha, beta
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_sm_fm_multiply
!!
!!   FUNCTION
!!     multiplies a sparse matrix with a full matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed & Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [f&j]
!!
!!*** **********************************************************************
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

    LOGICAL :: failure
    TYPE(cp_logger_type), POINTER :: logger
    CHARACTER(len=*), PARAMETER :: routineN='cp_sm_fm_multiply',&
         routineP=moduleN//':'//routineN
    REAL(kind=wp) :: my_alpha,my_beta
!FM  REAL(kind=wp) :: diff
!FM  TYPE(cp_fm_type), POINTER :: tmp_m

    logger=> cp_error_get_logger(error)

    my_alpha=1.0_wp
    my_beta=0.0_wp
    IF (PRESENT(alpha)) my_alpha=alpha
    IF (PRESENT(beta)) my_beta=beta
    failure=.FALSE.
!FM  NULLIFY(tmp_m)
    IF (ASSOCIATED(sparse_matrix%distribution_2d)) THEN
       IF (v_in%matrix_struct%para_env%num_pe>2) THEN
!FM     CALL cp_fm_create(tmp_m,v_out%matrix_struct,error=error)
!FM     CALL cp_sm_fm_multiply_general(sparse_matrix,v_in,tmp_m,ncol,para_env, &
!FM          alpha,beta,error)
          SELECT CASE(sparse_matrix%symmetry)
          CASE("symmetric")
             CALL cp_sm_fm_multiply_2d(sparse_matrix,v_in,v_out,ncol=ncol,&
                  alpha=my_alpha,beta=my_beta,transpose_sm=.FALSE.,&
                  do_diagonal_blocks=.TRUE.,antisymmetric=.FALSE.,error=error)
             CALL cp_sm_fm_multiply_2d(sparse_matrix,v_in,v_out,ncol=ncol,&
                  alpha=my_alpha,beta=1.0_wp,transpose_sm=.TRUE.,&
                  do_diagonal_blocks=.FALSE.,antisymmetric=.FALSE.,error=error)
          CASE("antisymmetric")
             CALL cp_sm_fm_multiply_2d(sparse_matrix,v_in,v_out,ncol=ncol,&
                  alpha=my_alpha,beta=my_beta,transpose_sm=.FALSE.,&
                  do_diagonal_blocks=.TRUE.,antisymmetric=.TRUE.,error=error)
             CALL cp_sm_fm_multiply_2d(sparse_matrix,v_in,v_out,ncol=ncol,&
                  alpha=my_alpha,beta=1.0_wp,transpose_sm=.TRUE.,&
                  do_diagonal_blocks=.FALSE.,antisymmetric=.TRUE.,error=error)
          CASE("none","no symmetry")
             CALL cp_sm_fm_multiply_2d(sparse_matrix,v_in,v_out,ncol=ncol,&
                  alpha=my_alpha,beta=my_beta,transpose_sm=.FALSE.,&
                  do_diagonal_blocks=.TRUE.,antisymmetric=.FALSE.,error=error)
          CASE DEFAULT
             ! possibly just OK (i.e. general matrix) just provide an empty slot
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "wrong matrix symmetry '"//sparse_matrix%symmetry//"' specification in "//&
                  CPSourceFileRef,&
                  error=error,failure=failure)
          END SELECT
!FM    CALL cp_fm_scale_and_add(alpha=1.0_wp,matrix_a=tmp_m,matrix_b=v_out,&
!FM         beta=-1.0_wp,error=error)
!FM    CALL cp_fm_maxval(tmp_m,diff)
!FM    IF (logger%para_env%source==logger%para_env%mepos) THEN
!FM       WRITE(unit=cp_logger_get_default_unit_nr(logger),fmt="(a,e16.8)")&
!FM            "matrix_multiply diff=",diff
!FM    END IF
!FM    CPPostcondition(diff<1.0e-14,cp_warning_level,routineP,error,failure)
       ELSE
          CALL cp_sm_fm_multiply_general(sparse_matrix,v_in,v_out,ncol,para_env, &
               alpha,beta,error)
       END IF
 ELSE
    CALL cp_log(logger, level=cp_warning_level, fromWhere=routineP , &
         message="Old matrix multiply called", local=.FALSE.)
    CALL cp_log(logger, level=cp_warning_level, fromWhere=routineP, &
         message=sparse_matrix%name)
    CALL cp_sm_fm_multiply_general(sparse_matrix,v_in,v_out,ncol,para_env, &
       alpha,beta,error)
 END IF

END SUBROUTINE cp_sm_fm_multiply
!***************************************************************************

!!****f* cp_sm_fm_interactions/cp_sm_plus_fm_fm_t [1.0] *
!!
!!   NAME
!!     cp_sm_plus_fm_fm_t
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_plus_fm_fm_t(sparse_matrix, matrix_v, matrix_g,&
!!         ncol, alpha, error)
!!       Type(real_matrix_type), Pointer:: sparse_matrix
!!       Type(cp_fm_type), Pointer:: matrix_v
!!       Type(cp_fm_type), Pointer, Optional:: matrix_g
!!       Integer, Intent (IN):: ncol
!!       Real(Kind=wp), Optional, Intent (IN):: alpha
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_plus_fm_fm_t
!!
!!   FUNCTION
!!     multiplies a sparse matrix with a full matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed & Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [f&j]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_plus_fm_fm_t(sparse_matrix,matrix_v,matrix_g,ncol,&
     alpha,error)
  TYPE(real_matrix_type), POINTER   :: sparse_matrix
  TYPE(cp_fm_type) , POINTER :: matrix_v
  TYPE(cp_fm_type) , POINTER, OPTIONAL :: matrix_g
  INTEGER, INTENT(IN)               :: ncol
  REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  TYPE(cp_logger_type), POINTER :: logger
  CHARACTER(len=*), PARAMETER :: routineN='cp_sm_plus_fm_fm_t',&
       routineP=moduleN//':'//routineN
!FM  REAL(kind=wp) :: diff
!FM  TYPE(real_matrix_type), POINTER :: tmp_m

  logger=> cp_error_get_logger(error)

  failure=.FALSE.
!FM  NULLIFY(tmp_m)
  IF (ASSOCIATED(sparse_matrix%distribution_2d)) THEN
     IF (matrix_v%matrix_struct%para_env%num_pe>1) THEN
!FM     CALL replicate_matrix(sparse_matrix,tmp_m,target_name="tmp_m",allocate_blocks=.true.)
!FM     call cp_sm_plus_fm_fm_t_general(sparse_matrix=tmp_m,matrix_v=matrix_v,&
!FM          matrix_g=matrix_g,ncol=ncol,alpha=alpha,error=error)

        CALL cp_sm_plus_fm_fm_t_2d(sm=sparse_matrix,matrix_v=matrix_v,&
             matrix_g=matrix_g,ncol=ncol,alpha=alpha, error=error)

!FM     CALL cp_sm_scale_and_add(alpha=1.0_wp,matrix_a=tmp_m,matrix_b=sparse_matrix,&
!FM          beta=-1.0_wp,error=error)
!FM     CALL cp_sm_maxabsval(tmp_m,diff,matrix_v%matrix_struct%para_env)
!FM     IF (logger%para_env%source==logger%para_env%mepos) THEN
!FM        WRITE(unit=cp_logger_get_default_unit_nr(logger),fmt="(a,e16.8)")&
!FM             "matrix_fm_fm_t diff=",diff
!FM     END IF
!FM     CPPostcondition(diff<1.0e-14,cp_warning_level,routineP,error,failure)
!FM     call deallocate_matrix(tmp_m)
     ELSE
        CALL cp_sm_plus_fm_fm_t_general(sparse_matrix=sparse_matrix,&
             matrix_v=matrix_v,&
             matrix_g=matrix_g,ncol=ncol,alpha=alpha,error=error)
     END IF
  ELSE
     CALL cp_log(logger, level=cp_warning_level, fromWhere=routineP , &
          message="Old matrix fm_fm_t called", local=.FALSE.)
     CALL cp_log(logger, level=cp_warning_level, fromWhere=routineP, &
          message=sparse_matrix%name)

     CALL cp_sm_plus_fm_fm_t_general(sparse_matrix=sparse_matrix,&
          matrix_v=matrix_v,&
          matrix_g=matrix_g,ncol=ncol,alpha=alpha,error=error)
  END IF

END SUBROUTINE cp_sm_plus_fm_fm_t
!***************************************************************************

!!****f* cp_sm_fm_interactions/cp_sm_fm_multiply_2d [1.0] *
!!
!!   NAME
!!     cp_sm_fm_multiply_2d
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele & Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [j&f]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sm_fm_multiply_2d(sm,fm_in,fm_out, alpha,beta,transpose_sm, &
       do_diagonal_blocks,antisymmetric, ncol, error)

    TYPE(real_matrix_type), POINTER :: sm
    TYPE(cp_fm_type), POINTER :: fm_in, fm_out
    LOGICAL, INTENT(in) :: transpose_sm, do_diagonal_blocks,antisymmetric
    REAL(kind=wp), INTENT(in) :: alpha, beta
    INTEGER, INTENT(in), OPTIONAL :: ncol
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    INTEGER :: fm_ncol_global, handle, handle2
    INTEGER :: nprow,npcol,fm_ncol_local, fm_nrow_local
    INTEGER :: nchunk,chunk_size, gindex
    INTEGER :: nchunk_local_in,nchunk_local_out

    CHARACTER(len=*), PARAMETER :: routineN='cp_sm_fm_multiply_2d',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure, transpose_block
    INTEGER :: sm_ncol_local, myprow,mypcol, stat,&
         mepos,ip,icol,irow,ipcol,iprow, sm_nrow_local,ichunk,&
         ichunk_offset, ichunk_ordered, irow_chunk, iblock_atomic,&
         irow_global, iblock, iblock_global, &
         chunk_in, chunk_out, iblock_row, iblock_col, sbnrow,sbncol,&
         chunk_recv_in, chunk_recv_out, row_dest, row_src, col_dest,&
         col_src, sm_nblock_row, sm_nblock_col, fm_ncol_block, fm_nrow_block,&
         fm_first_prow, fm_first_pcol
    INTEGER, DIMENSION(4) :: reqs
    INTEGER, DIMENSION(:), POINTER :: chunk_of_col, ncol_local_of_chunk,&
         col_offset_in_chunk, row_offset_in_chunk, fm_row_indices, fm_col_indices,&
         local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution, ncol_global_of_chunk
    INTEGER, DIMENSION(:,:), POINTER :: blacs2mpi
    REAL(kind=wp), DIMENSION(:,:), POINTER :: sparse_block
    REAL(kind=wp) :: flops,flops_l,alpha_ji

    TYPE(real_block_node_type), POINTER :: block_node
    TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: v_in, v_out
    TYPE(cp_para_env_type), POINTER :: fm_para_env
    TYPE (cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","Gflops",handle)
    logger => cp_error_get_logger(error)
    NULLIFY(chunk_of_col, ncol_local_of_chunk,&
         col_offset_in_chunk, row_offset_in_chunk, fm_row_indices, fm_col_indices,&
         local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution, ncol_global_of_chunk)
    NULLIFY(blacs2mpi,sparse_block,block_node,v_in,v_out,fm_para_env)
    CPPrecondition(ASSOCIATED(sm),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(fm_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(fm_out),cp_failure_level,routineP,error,failure)
    CALL cp_assert(cp_fm_struct_equivalent(fm_in%matrix_struct,fm_out%matrix_struct,error=error),&
         cp_failure_level, cp_assertion_failed, routineP, "input and output matrixes must be compatible "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.transpose_sm) THEN
       nprow         =fm_in%matrix_struct%context%num_pe(1)
       npcol         =fm_in%matrix_struct%context%num_pe(2)
       myprow        =fm_in%matrix_struct%context%mepos(1)
       mypcol        =fm_in%matrix_struct%context%mepos(2)
       blacs2mpi     => fm_in%matrix_struct%context%blacs2mpi
       
       local_atomic_rows => sm%distribution_2d%local_rows
       local_atomic_cols => sm%distribution_2d%local_cols
       atomic_row_distribution => sm%distribution_2d%row_distribution
       atomic_col_distribution => sm%distribution_2d%col_distribution

       CALL get_matrix_info(sm,nblock_row=sm_nblock_row,&
            nblock_col=sm_nblock_col, last_row=sm_last_row,&
            last_col=sm_last_col, first_row=sm_first_row,&
            first_col=sm_first_col)
    ELSE
       npcol         =fm_in%matrix_struct%context%num_pe(1)
       nprow         =fm_in%matrix_struct%context%num_pe(2)
       mypcol        =fm_in%matrix_struct%context%mepos(1)
       myprow        =fm_in%matrix_struct%context%mepos(2)
       ALLOCATE(blacs2mpi(0:SIZE(fm_in%matrix_struct%context%blacs2mpi,2)-1,&
            0:SIZE(fm_in%matrix_struct%context%blacs2mpi,1)-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iprow=0,nprow-1
          DO ipcol=0,npcol-1
             blacs2mpi(iprow,ipcol)=fm_in%matrix_struct%context%blacs2mpi(ipcol,iprow)
          END DO
       END DO
       
       local_atomic_cols => sm%distribution_2d%local_rows
       local_atomic_rows => sm%distribution_2d%local_cols
       atomic_col_distribution => sm%distribution_2d%row_distribution
       atomic_row_distribution => sm%distribution_2d%col_distribution

       CALL get_matrix_info(sm,nblock_col=sm_nblock_row,&
            nblock_row=sm_nblock_col, last_col=sm_last_row,&
            last_row=sm_last_col, first_col=sm_first_row,&
            first_row=sm_first_col)
    END IF

    CALL cp_fm_get_info(fm_in,ncol_global=fm_ncol_global,&
         ncol_local=fm_ncol_local, nrow_local=fm_nrow_local,&
         row_indices=fm_row_indices, col_indices=fm_col_indices,&
         para_env=fm_para_env,error=error)
    IF (PRESENT(ncol)) fm_ncol_global=ncol

    nchunk        =lcm(nprow,npcol)
    IF (nchunk==nprow.OR.nchunk==npcol)nchunk=2*nchunk
    chunk_size = (fm_ncol_global+nchunk-1)/ nchunk
    nchunk_local_in = nchunk / nprow
    nchunk_local_out= nchunk / npcol
    CPPrecondition(nchunk_local_in>1,cp_failure_level,routineP,error,failure)
    CPPrecondition(nchunk_local_out>1,cp_failure_level,routineP,error,failure)

    ! number of global fm columns belonging to a chunk
    ALLOCATE(ncol_global_of_chunk(nchunk),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ncol_global_of_chunk=0
    DO ichunk=1,nchunk
       ncol_global_of_chunk(ichunk)=fm_ncol_global/nchunk
       IF (ichunk<=MODULO(fm_ncol_global,nchunk)) THEN
          ncol_global_of_chunk(ichunk)=ncol_global_of_chunk(ichunk)+1
       END IF
    END DO

    ! where in v_out begins the given atomic block
    ALLOCATE(row_offset_in_chunk(sm_nblock_row),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    row_offset_in_chunk=-HUGE(0)
    sm_nrow_local=0
    DO iblock=1,SIZE(local_atomic_rows)
       iblock_global=local_atomic_rows(iblock)
       row_offset_in_chunk(iblock_global)=sm_nrow_local
       sm_nrow_local=sm_nrow_local+sm_last_row(iblock_global)-sm_first_row(iblock_global)+1
    END DO

    ! where in v_in begins the given atomic block
    ALLOCATE(col_offset_in_chunk(sm_nblock_col),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    col_offset_in_chunk=-HUGE(0)
    sm_ncol_local=0
    DO iblock=1,SIZE(local_atomic_cols)
       iblock_global=local_atomic_cols(iblock)
       col_offset_in_chunk(iblock_global)=sm_ncol_local
       sm_ncol_local=sm_ncol_local+sm_last_col(iblock_global)-&
            sm_first_col(iblock_global)+1
    END DO

    ! alloc & init v_in
    ALLOCATE(v_in (nchunk_local_in+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ichunk=1,nchunk_local_in+1
       ALLOCATE(v_in(ichunk)%array(chunk_size,sm_ncol_local), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    CALL fm2chunk(sm,fm_in,v_in,shift=0,hole_pos=nchunk_local_in+1,&
         sm_transposed=transpose_sm,alpha=alpha,ncol=fm_ncol_global,&
         error=error)
    alpha_ji=1.0_wp
    IF(antisymmetric) alpha_ji=-1.0_wp

    ! alloc & init v_out
    ALLOCATE(v_out(nchunk_local_out+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ichunk=1,nchunk_local_out+1
       ALLOCATE(v_out(ichunk)%array(chunk_size,sm_nrow_local), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       v_out(ichunk)%array=0.0_wp
    END DO

    row_dest=blacs2mpi(MODULO(myprow-1,nprow),mypcol)
    row_src=blacs2mpi(MODULO(myprow+1,nprow),mypcol)
    col_dest=blacs2mpi(myprow,MODULO(mypcol-1,npcol))
    col_src=blacs2mpi(myprow,MODULO(mypcol+1,npcol))
    
    ! let the magic of the distribution play the game
    flops=0.0_wp
    DO ichunk=1,nchunk
       ! go circular
       CALL timeset(routineN//"_local","I","Gflops",handle2)
       flops_l=0.0_wp
       chunk_in =MOD(ichunk-1,nchunk_local_in+1)+1
       chunk_out=MOD(ichunk-1,nchunk_local_out+1)+1
       gindex=MODULO(nchunk_local_out*mypcol+nchunk_local_in*myprow+ichunk-1,&
            nchunk)+1

       !   *** Traverse all block nodes of the sparse matrix ***
       !   no simple parallel do because we write both the iblock_row and iblock_col

       DO iblock_row=1,sm%nblock_row

          block_node => first_block_node(sm,iblock_row)

          DO WHILE (ASSOCIATED(block_node))
             transpose_block=transpose_sm

             CALL get_block_node(block_node=block_node,&
                  block_col=iblock_col,&
                  BLOCK=sparse_block)
             
             IF (iblock_col/=iblock_row.OR.do_diagonal_blocks) THEN
                sbnrow=sm%last_row(iblock_row)- &
                     sm%first_row(iblock_row)+1
                sbncol=sm%last_col(iblock_col)- &
                     sm%first_col(iblock_col)+1

                flops_l=flops_l+2.0_wp*sbnrow*sbncol*ncol_global_of_chunk(gindex)

                IF (.NOT.(sm%distribution_2d%row_distribution(iblock_row)==sm%distribution_2d%blacs_env%mepos(1).AND.&
                     sm%distribution_2d%col_distribution(iblock_col)==sm%distribution_2d%blacs_env%mepos(2))) THEN
                   transpose_block=.NOT.transpose_block
                END IF
                IF (.NOT.transpose_block) THEN
                   CALL DGEMM('N','T',ncol_global_of_chunk(gindex),sbnrow,sbncol,1.0_wp, &
                        v_in(chunk_in)%array(1,col_offset_in_chunk(iblock_col)+1), &
                        chunk_size, &
                        sparse_block(1,1),SIZE(sparse_block,1), &
                        1.0_wp,v_out(chunk_out)%array(1,row_offset_in_chunk(iblock_row)+1), &
                        chunk_size)
                   
                ELSE
                   CALL DGEMM('N','N',ncol_global_of_chunk(gindex),sbncol,sbnrow,alpha_ji, &
                        v_in(chunk_in)%array(1,col_offset_in_chunk(iblock_row)+1), &
                        chunk_size, &
                        sparse_block(1,1),SIZE(sparse_block,1), &
                        1.0_wp,v_out(chunk_out)%array(1,row_offset_in_chunk(iblock_col)+1), &
                        chunk_size)
                END IF
             END IF

             block_node => next_block_node(block_node)

          END DO

       END DO
       flops=flops+flops_l
       CALL timestop(flops_l*1.0e-9,handle2)

       chunk_recv_in=MODULO(chunk_in-2,nchunk_local_in+1)+1
       chunk_recv_out=MODULO(chunk_out-2,nchunk_local_out+1)+1

       IF (ichunk>1) CALL mp_waitall(reqs)
       IF (ichunk==nchunk) EXIT

       CALL mp_isendrecv(msgin=v_in(chunk_in)%array,dest=row_dest,&
            msgout=v_in(chunk_recv_in)%array,source=row_src,&
            comm=fm_para_env%group, tag=1, send_request=reqs(1), &
            recv_request=reqs(2))
       CALL mp_isendrecv(msgin=v_out(chunk_out)%array,dest=col_dest,&
            msgout=v_out(chunk_recv_out)%array,source=col_src,&
            comm=fm_para_env%group, tag=2, send_request=reqs(3), &
            recv_request=reqs(4))
    ENDDO

    DEALLOCATE(ncol_global_of_chunk,row_offset_in_chunk,col_offset_in_chunk,&
         stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    DO ichunk=1,SIZE(v_in)
       DEALLOCATE(v_in(ichunk)%array,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END DO
    DEALLOCATE(v_in, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    ! reshuffle v_out
    CALL chunk2fm(sm=sm,fm=fm_out,chunks=v_out,shift=-1,hole_pos=chunk_recv_out,&
         sm_transposed=.NOT.transpose_sm, beta=beta,ncol=fm_ncol_global,error=error)

    IF (transpose_sm) THEN
       DEALLOCATE(blacs2mpi,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    DO ichunk=1,SIZE(v_out)
       DEALLOCATE(v_out(ichunk)%array,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END DO
    DEALLOCATE(v_out,stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(flops*1.0e-9_wp,handle)

  END SUBROUTINE cp_sm_fm_multiply_2d
!**************************************************************************

!!****f* cp_sm_fm_interactions/cp_sm_plus_fm_fm_t_2d [1.0] *
!!
!!   NAME
!!     cp_sm_plus_fm_fm_t_2d
!!
!!   FUNCTION
!!     computes sparse=sparse+alpha*v*g^T  
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele & Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [j&f]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sm_plus_fm_fm_t_2d(sm,matrix_v,matrix_g,ncol,&
       alpha,error)
    TYPE(real_matrix_type), POINTER   :: sm
    TYPE(cp_fm_type) , POINTER :: matrix_v
    TYPE(cp_fm_type) , POINTER, OPTIONAL :: matrix_g
    INTEGER, INTENT(IN),OPTIONAL :: ncol
    REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    INTEGER :: fm_ncol_global, handle, handle2
    INTEGER :: nprow,npcol,fm_ncol_local, fm_nrow_local
    INTEGER :: nchunk,chunk_size, gindex
    INTEGER :: nchunk_local_in,nchunk_local_out

    CHARACTER(len=*), PARAMETER :: routineN='cp_sm_plus_fm_fm_t_2d',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure, transpose_block
    INTEGER :: sm_ncol_local, myprow,mypcol, stat,&
         mepos,ip,icol,irow,ipcol,iprow, sm_nrow_local,ichunk,&
         ichunk_offset, ichunk_ordered, irow_chunk, iblock_atomic,&
         irow_global, iblock, iblock_global, &
         chunk_in, chunk_out, iblock_row, iblock_col, sbnrow,sbncol,&
         chunk_recv_in, chunk_recv_out, row_dest, row_src, col_dest,&
         col_src, sm_nblock_row, sm_nblock_col, fm_ncol_block, fm_nrow_block,&
         fm_first_prow, fm_first_pcol
    INTEGER, DIMENSION(4) :: reqs
    INTEGER, DIMENSION(:), POINTER :: chunk_of_col, ncol_local_of_chunk,&
         col_offset_in_chunk, row_offset_in_chunk, fm_row_indices, fm_col_indices,&
         local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution, ncol_global_of_chunk
    INTEGER, DIMENSION(:,:), POINTER :: blacs2mpi
    REAL(kind=wp), DIMENSION(:,:), POINTER :: sparse_block
    REAL(kind=wp) :: flops,flops_l,my_alpha
    TYPE(cp_fm_type), POINTER :: my_matrix_g

    TYPE(real_block_node_type), POINTER :: block_node
    TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: v_in, v_out
    TYPE(cp_para_env_type), POINTER :: fm_para_env
    TYPE (cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","Gflops",handle)
    logger => cp_error_get_logger(error)
    NULLIFY(chunk_of_col, ncol_local_of_chunk,&
         col_offset_in_chunk, row_offset_in_chunk, fm_row_indices, fm_col_indices,&
         local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution, ncol_global_of_chunk)
    NULLIFY(blacs2mpi,sparse_block,my_matrix_g,block_node,v_in,v_out,fm_para_env)

    my_alpha=1.0_wp
    IF (PRESENT(alpha)) my_alpha=alpha
    my_matrix_g => matrix_v
    IF (PRESENT(matrix_g)) my_matrix_g => matrix_g

    CPPrecondition(ASSOCIATED(sm),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_v),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(my_matrix_g),cp_failure_level,routineP,error,failure)

    CALL cp_assert(cp_fm_struct_equivalent(matrix_v%matrix_struct,my_matrix_g%matrix_struct,error=error),&
         cp_failure_level, cp_assertion_failed, routineP, "input and output matrixes must be compatible "//&
         CPSourceFileRef,&
         error,failure)
    CALL cp_assert(sm%symmetry=="none".OR.sm%symmetry=="symmetric",cp_failure_level,&
         cp_assertion_failed,routineP,"matrix symmetry incorrect in "//&
         CPSourceFileRef,&
         error,failure) 

    nprow         =matrix_v%matrix_struct%context%num_pe(1)
    npcol         =matrix_v%matrix_struct%context%num_pe(2)
    myprow        =matrix_v%matrix_struct%context%mepos(1)
    mypcol        =matrix_v%matrix_struct%context%mepos(2)
    blacs2mpi     => matrix_v%matrix_struct%context%blacs2mpi

    local_atomic_rows => sm%distribution_2d%local_rows
    local_atomic_cols => sm%distribution_2d%local_cols
    atomic_row_distribution => sm%distribution_2d%row_distribution
    atomic_col_distribution => sm%distribution_2d%col_distribution

    CALL get_matrix_info(sm,nblock_row=sm_nblock_row,&
         nblock_col=sm_nblock_col, last_row=sm_last_row,&
         last_col=sm_last_col, first_row=sm_first_row,&
         first_col=sm_first_col)

    CALL cp_fm_get_info(matrix_v,ncol_global=fm_ncol_global,&
         ncol_local=fm_ncol_local, nrow_local=fm_nrow_local,&
         row_indices=fm_row_indices, col_indices=fm_col_indices,&
         para_env=fm_para_env,error=error)
    IF (PRESENT(ncol)) fm_ncol_global=ncol

    nchunk        =lcm(nprow,npcol)
    IF (nchunk==nprow.OR.nchunk==npcol)nchunk=2*nchunk
    chunk_size = (fm_ncol_global+nchunk-1)/ nchunk
    nchunk_local_in = nchunk / nprow
    nchunk_local_out= nchunk / npcol

    CPPrecondition(nchunk_local_in>1,cp_failure_level,routineP,error,failure)
    CPPrecondition(nchunk_local_out>1,cp_failure_level,routineP,error,failure)

    ! number of global fm columns belonging to a chunk
    ALLOCATE(ncol_global_of_chunk(nchunk),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ncol_global_of_chunk=0
    DO ichunk=1,nchunk
       ncol_global_of_chunk(ichunk)=fm_ncol_global/nchunk
       IF (ichunk<=MODULO(fm_ncol_global,nchunk)) THEN
          ncol_global_of_chunk(ichunk)=ncol_global_of_chunk(ichunk)+1
       END IF
    END DO

    ! where in v_out begins the given atomic block
    ALLOCATE(row_offset_in_chunk(sm_nblock_row),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    row_offset_in_chunk=-HUGE(0)
    sm_nrow_local=0
    DO iblock=1,SIZE(local_atomic_rows)
       iblock_global=local_atomic_rows(iblock)
       row_offset_in_chunk(iblock_global)=sm_nrow_local
       sm_nrow_local=sm_nrow_local+sm_last_row(iblock_global)-sm_first_row(iblock_global)+1
    END DO

    ! where in v_in begins the given atomic block
    ALLOCATE(col_offset_in_chunk(sm_nblock_col),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    col_offset_in_chunk=-HUGE(0)
    sm_ncol_local=0
    DO iblock=1,SIZE(local_atomic_cols)
       iblock_global=local_atomic_cols(iblock)
       col_offset_in_chunk(iblock_global)=sm_ncol_local
       sm_ncol_local=sm_ncol_local+sm_last_col(iblock_global)-&
            sm_first_col(iblock_global)+1
    END DO

    ! alloc & init v_in
    ALLOCATE(v_in (nchunk_local_in+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ichunk=1,nchunk_local_in+1
       ALLOCATE(v_in(ichunk)%array(chunk_size,sm_ncol_local), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    CALL fm2chunk(sm,my_matrix_g,v_in,shift=0,hole_pos=nchunk_local_in+1,&
         sm_transposed=.FALSE.,alpha=my_alpha,ncol=fm_ncol_global,&
         error=error)

    ! alloc & init v_out
    ALLOCATE(v_out(nchunk_local_out+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ichunk=1,nchunk_local_out+1
       ALLOCATE(v_out(ichunk)%array(chunk_size,sm_nrow_local), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    CALL fm2chunk(sm,matrix_v,v_out,shift=0,hole_pos=nchunk_local_out+1,&
         sm_transposed=.TRUE.,alpha=1.0_wp,ncol=fm_ncol_global,&
         error=error)

    row_dest=blacs2mpi(MODULO(myprow-1,nprow),mypcol)
    row_src=blacs2mpi(MODULO(myprow+1,nprow),mypcol)
    col_dest=blacs2mpi(myprow,MODULO(mypcol-1,npcol))
    col_src=blacs2mpi(myprow,MODULO(mypcol+1,npcol))

    ! let the magic of the distribution play the game
    flops=0.0_wp
    DO ichunk=1,nchunk
       ! go circular
       CALL timeset(routineN//"_local","I","Gflops",handle2)
       flops_l=0.0_wp
       chunk_in =MOD(ichunk-1,nchunk_local_in+1)+1
       chunk_out=MOD(ichunk-1,nchunk_local_out+1)+1
       gindex=MODULO(nchunk_local_out*mypcol+nchunk_local_in*myprow+ichunk-1,&
            nchunk)+1

       !   *** Traverse all block nodes of the sparse matrix ***
       !   no simple parallel do because we write both the iblock_row and iblock_col

       DO iblock_row=1,sm%nblock_row

          block_node => first_block_node(sm,iblock_row)

          DO WHILE (ASSOCIATED(block_node))

             CALL get_block_node(block_node=block_node,&
                  block_col=iblock_col,&
                  BLOCK=sparse_block)

             sbnrow=sm%last_row(iblock_row)- &
                  sm%first_row(iblock_row)+1
             sbncol=sm%last_col(iblock_col)- &
                  sm%first_col(iblock_col)+1

             flops_l=flops_l+2.0_wp*sbnrow*sbncol*ncol_global_of_chunk(gindex)

             IF (.NOT.(sm%distribution_2d%row_distribution(iblock_row)==sm%distribution_2d%blacs_env%mepos(1).AND.&
                  sm%distribution_2d%col_distribution(iblock_col)==sm%distribution_2d%blacs_env%mepos(2))) THEN
                transpose_block=.TRUE.
             ELSE
                transpose_block=.FALSE.
             END IF
             IF (.NOT.transpose_block) THEN
                CALL DGEMM('T','N',sbnrow,sbncol,ncol_global_of_chunk(gindex),1.0_wp, &
                     v_out(chunk_out)%array(1,row_offset_in_chunk(iblock_row)+1),&
                     chunk_size,&
                     v_in(chunk_in)%array(1,col_offset_in_chunk(iblock_col)+1), &
                     chunk_size, &
                     1.0_wp,sparse_block(1,1),SIZE(sparse_block,1))

             ELSE
                CALL DGEMM('T','N',sbnrow,sbncol,ncol_global_of_chunk(gindex),1.0_wp, &
                     v_in(chunk_in)%array(1,col_offset_in_chunk(iblock_row)+1), &
                     chunk_size, &
                     v_out(chunk_out)%array(1,row_offset_in_chunk(iblock_col)+1), &
                     chunk_size,&
                     1.0_wp,sparse_block(1,1),SIZE(sparse_block,1))
             END IF

             block_node => next_block_node(block_node)

          END DO

       END DO
       flops=flops+flops_l
       CALL timestop(flops_l*1.0e-9,handle2)

       chunk_recv_in=MODULO(chunk_in-2,nchunk_local_in+1)+1
       chunk_recv_out=MODULO(chunk_out-2,nchunk_local_out+1)+1

       IF (ichunk>1) CALL mp_waitall(reqs)
       IF (ichunk==nchunk) EXIT

       CALL mp_isendrecv(msgin=v_in(chunk_in)%array,dest=row_dest,&
            msgout=v_in(chunk_recv_in)%array,source=row_src,&
            comm=fm_para_env%group, tag=1, send_request=reqs(1), &
            recv_request=reqs(2))
       CALL mp_isendrecv(msgin=v_out(chunk_out)%array,dest=col_dest,&
            msgout=v_out(chunk_recv_out)%array,source=col_src,&
            comm=fm_para_env%group, tag=2, send_request=reqs(3), &
            recv_request=reqs(4))
    ENDDO

    DEALLOCATE(ncol_global_of_chunk,row_offset_in_chunk,col_offset_in_chunk,&
         stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    DO ichunk=1,SIZE(v_in)
       DEALLOCATE(v_in(ichunk)%array,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END DO
    DEALLOCATE(v_in, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    DO ichunk=1,SIZE(v_out)
       DEALLOCATE(v_out(ichunk)%array,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END DO
    DEALLOCATE(v_out,stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(flops*1.0e-9_wp,handle)

  END SUBROUTINE cp_sm_plus_fm_fm_t_2d
!***************************************************************************

!!****f* cp_sm_fm_interactions/chunk_of_fm_index [1.0] *
!!
!!   NAME
!!     chunk_of_fm_index
!!
!!   FUNCTION
!!     returns the chunk (starting at 1) in which the given global index is using the
!!     distribution size_small+1,size_small+1,...,size_small,size_small,...
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fm_index: index of the full matrix (1 based)
!!     - fm_size: number of elements of the full matrix
!!     - nchunk: number of chunks you want
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele & Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [j&f]
!!
!!*** **********************************************************************
  FUNCTION chunk_of_fm_index(fm_index,fm_size,nchunk,error) RESULT(res)
    INTEGER ::res
    INTEGER, INTENT(in) :: fm_index, fm_size, nchunk
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='chunk_of_fm_index',&
         routineP=moduleN//':'//routineN
    INTEGER :: size_small,n_large

    failure=.FALSE.

    CPPrecondition(fm_index<=fm_size,cp_failure_level,routineP,error,failure)
    CPPrecondition(nchunk>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(0<fm_index,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       size_small=fm_size/nchunk
       n_large=MODULO(fm_size,nchunk)
       res=(fm_index-1)/(size_small+1)+1
       IF (res>n_large) res= n_large+(fm_index-1-n_large*(size_small+1))/size_small+1
    ELSE
       res=-HUGE(0)
    END IF
  END FUNCTION chunk_of_fm_index
  !***************************************************************************

!!****f* cp_sm_fm_interactions/fm_reshuffle_create_layout [1.0] *
!!
!!   NAME
!!     fm_reshuffle_create_layout
!!
!!   FUNCTION
!!     returns various info about the layout needed by the multiplication
!!
!!   NOTES
!!     aloocates the arrays
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele & Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [j&f]
!!
!!*** **********************************************************************
  SUBROUTINE fm_reshuffle_create_layout(sm,fm, transpose_layout, shift,&
       ncol,gindex_of_chunk, ncol_global_of_chunk, ipcol_of_chunk_col,&
       chunk_of_col, ncol_local_of_chunk, pcol_of_fm_row, f2c_send_count, &
       f2c_recv_count,chunk_on_prow, error)
    TYPE(real_matrix_type), POINTER :: sm
    TYPE(cp_fm_type), POINTER :: fm
    LOGICAL, INTENT(in) :: transpose_layout
    INTEGER, INTENT(in) :: shift,ncol
    INTEGER, DIMENSION(:), POINTER :: gindex_of_chunk, ncol_global_of_chunk,&
         chunk_of_col, ncol_local_of_chunk,&
         pcol_of_fm_row, f2c_send_count, f2c_recv_count
    INTEGER, DIMENSION(:,:), POINTER :: chunk_on_prow,&
         ipcol_of_chunk_col
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_reshuffle_create_layout',&
         routineP=moduleN//':'//routineN
    INTEGER :: myprow, mypcol,fm_ncol_global, iblock, handle, handle1,icol2,bsize,irow2
    INTEGER :: nprow,npcol,fm_ncol_local, fm_nrow_local
    INTEGER :: nchunk,chunk_size, ipcol, ichunk, stat, ichunk_offset, icol, iblock_atomic
    INTEGER :: nchunk_local_in, nchunk_local_out, irow, irow_global, iprow,&
         first_g_col,fm_npcol,nl,ss
    INTEGER, DIMENSION(:), POINTER :: fm_row_indices, fm_col_indices,&
         sm_last_row, sm_first_row, sm_first_col,sm_last_col,&
         local_atomic_cols, local_atomic_rows, atomic_col_distribution,&
         atomic_row_distribution
    INTEGER, DIMENSION(:,:), POINTER :: blacs2mpi
    INTEGER :: sm_ncol_local, fm_ncol_block, fm_nrow_block,&
         fm_first_prow, fm_first_pcol, sm_nblock_col, sm_nblock_row,&
         sm_nrow_local, i
    TYPE(cp_para_env_type), POINTER :: fm_para_env
#ifdef __SCALAPACK
    INTEGER, EXTERNAL :: indxg2p
#endif

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(fm_row_indices, fm_col_indices,&
         sm_last_row, sm_first_row, sm_first_col,sm_last_col,&
         local_atomic_cols, local_atomic_rows, atomic_col_distribution,&
         atomic_row_distribution)
    NULLIFY(blacs2mpi,fm_para_env)

    IF (.NOT.transpose_layout) THEN
       nprow         =fm%matrix_struct%context%num_pe(1)
       npcol         =fm%matrix_struct%context%num_pe(2)
       myprow        =fm%matrix_struct%context%mepos(1)
       mypcol        =fm%matrix_struct%context%mepos(2)
       blacs2mpi     => fm%matrix_struct%context%blacs2mpi

       local_atomic_rows => sm%distribution_2d%local_rows
       local_atomic_cols => sm%distribution_2d%local_cols
       atomic_row_distribution => sm%distribution_2d%row_distribution
       atomic_col_distribution => sm%distribution_2d%col_distribution

       CALL get_matrix_info(sm,nblock_row=sm_nblock_row,&
            nblock_col=sm_nblock_col, last_row=sm_last_row,&
            last_col=sm_last_col, first_row=sm_first_row,&
            first_col=sm_first_col)
       
    ELSE
       npcol         =fm%matrix_struct%context%num_pe(1)
       nprow         =fm%matrix_struct%context%num_pe(2)
       mypcol        =fm%matrix_struct%context%mepos(1)
       myprow        =fm%matrix_struct%context%mepos(2)
       ALLOCATE(blacs2mpi(0:SIZE(fm%matrix_struct%context%blacs2mpi,2)-1,&
            0:SIZE(fm%matrix_struct%context%blacs2mpi,1)-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iprow=0,nprow-1
          DO ipcol=0,npcol-1
             blacs2mpi(iprow,ipcol)=fm%matrix_struct%context%blacs2mpi(ipcol,iprow)
          END DO
       END DO

       local_atomic_cols => sm%distribution_2d%local_rows
       local_atomic_rows => sm%distribution_2d%local_cols
       atomic_col_distribution => sm%distribution_2d%row_distribution
       atomic_row_distribution => sm%distribution_2d%col_distribution

       CALL get_matrix_info(sm,nblock_col=sm_nblock_row,&
            nblock_row=sm_nblock_col, last_col=sm_last_row,&
            last_row=sm_last_col, first_col=sm_first_row,&
            first_row=sm_first_col)
    END IF

    fm_first_prow=fm%matrix_struct%first_p_pos(1)
    fm_first_pcol=fm%matrix_struct%first_p_pos(2)
    fm_npcol=fm%matrix_struct%context%num_pe(2)
    CALL cp_fm_get_info(fm,&
         ncol_local=fm_ncol_local, nrow_local=fm_nrow_local,&
         row_indices=fm_row_indices, col_indices=fm_col_indices,&
         para_env=fm_para_env,&
         nrow_block=fm_nrow_block,ncol_block=fm_ncol_block,&
         error=error)
    fm_ncol_global=ncol
    sm_nrow_local=0
    DO iblock=1,SIZE(local_atomic_rows)
       sm_nrow_local=sm_nrow_local+&
            sm_last_row(local_atomic_rows(iblock))-&
            sm_first_row(local_atomic_rows(iblock))+1
    END DO
    sm_ncol_local=0
    DO iblock=1,SIZE(local_atomic_cols)
       sm_ncol_local=sm_ncol_local+&
            sm_last_col(local_atomic_cols(iblock))-&
            sm_first_col(local_atomic_cols(iblock))+1
    END DO
       
    nchunk        =lcm(nprow,npcol)
    IF (nchunk==nprow.OR.nchunk==npcol)nchunk=2*nchunk
    chunk_size = (fm_ncol_global+nchunk-1)/ nchunk
    nchunk_local_in = nchunk / nprow
    nchunk_local_out= nchunk / npcol
    
    ! tells you on which processor row is the chunk, depending on pcol
    ALLOCATE(chunk_on_prow(nchunk,0:npcol-1), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ipcol=0,npcol-1
       DO ichunk=1,nchunk
          chunk_on_prow(ichunk,ipcol)=MODULO((2*nchunk+ichunk-1-ipcol*nchunk_local_out-shift)/nchunk_local_in,nprow)
       ENDDO
    ENDDO

    ! global indices of the chunks that are stored locally in v_in
    ALLOCATE(gindex_of_chunk(nchunk_local_in),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ichunk_offset=nchunk_local_out*mypcol+nchunk_local_in*myprow+shift
    DO ichunk=1,nchunk_local_in
       gindex_of_chunk(ichunk)=MODULO(ichunk-1+ichunk_offset,nchunk)+1
    END DO

    ! number of global fm columns belonging to a chunk
    ALLOCATE(ncol_global_of_chunk(nchunk),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ncol_global_of_chunk=0
    DO ichunk=1,nchunk
       ncol_global_of_chunk(ichunk)=fm_ncol_global/nchunk
       IF (ichunk<=MODULO(fm_ncol_global,nchunk)) THEN
          ncol_global_of_chunk(ichunk)=ncol_global_of_chunk(ichunk)+1
       END IF
    END DO

    ! global indices of the columns in the given chunk
    ALLOCATE(ipcol_of_chunk_col(chunk_size,nchunk_local_in),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ipcol_of_chunk_col=-HUGE(0)
    DO ichunk=1,nchunk_local_in
       nl=MODULO(fm_ncol_global,nchunk)
       ss=fm_ncol_global/nchunk
       IF (gindex_of_chunk(ichunk)>nl) THEN
          first_g_col=nl*(ss+1)+(gindex_of_chunk(ichunk)-nl-1)*ss+1
       ELSE
          first_g_col=(gindex_of_chunk(ichunk)-1)*(ss+1)+1
       END IF
!FM       CPPrecondition(first_g_col==SUM(ncol_global_of_chunk(1:gindex_of_chunk(ichunk)))-ncol_global_of_chunk(gindex_of_chunk(ichunk))+1,cp_failure_level,routineP,error,failure)
       bsize=fm_ncol_block-MOD(first_g_col-1,fm_ncol_block)
#ifdef __SCALAPACK
       ipcol = indxg2p(first_g_col,fm_ncol_block,&
            mypcol, fm_first_pcol,fm_npcol)
#else
       ipcol =0
#endif
       DO icol=1,MIN(bsize,ncol_global_of_chunk(gindex_of_chunk(ichunk)))
          ipcol_of_chunk_col(icol,ichunk)=ipcol
       END DO
       DO icol=MIN(bsize,ncol_global_of_chunk(gindex_of_chunk(ichunk)))+1,&
            ncol_global_of_chunk(gindex_of_chunk(ichunk)),fm_ncol_block
          ipcol=MODULO(ipcol+1,fm_npcol)
          DO icol2=icol,MIN(icol+fm_ncol_block-1,ncol_global_of_chunk(gindex_of_chunk(ichunk)))
             ipcol_of_chunk_col(icol2,ichunk)=ipcol
          END DO
       END DO
    END DO

   ! number of local fm columns belonging to a chunk
    ! index of the chunk that owns the given local fm column
    ALLOCATE(chunk_of_col(fm_ncol_local), ncol_local_of_chunk(nchunk),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ncol_local_of_chunk=0
    chunk_of_col=-HUGE(0)
    DO icol=1,fm_ncol_local
       IF (fm_col_indices(icol)>fm_ncol_global) EXIT
       chunk_of_col(icol)=chunk_of_fm_index(fm_col_indices(icol),fm_ncol_global,&
            nchunk,error=error)
       ncol_local_of_chunk(chunk_of_col(icol))=ncol_local_of_chunk(chunk_of_col(icol))+1
    END DO

    ! processor that should receive the given full matrix row
    ALLOCATE(pcol_of_fm_row(fm_nrow_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    iblock_atomic=1
    DO irow=1,fm_nrow_local
       irow_global=fm_row_indices(irow)
       DO WHILE (sm_last_col(iblock_atomic)< irow_global)
          iblock_atomic=iblock_atomic+1
       END DO
       pcol_of_fm_row(irow)=atomic_col_distribution(iblock_atomic)
    END DO

    ! size of the message being sent by a processor to a processor
    ALLOCATE(f2c_send_count(0:fm_para_env%num_pe-1),&
         f2c_recv_count(0:fm_para_env%num_pe-1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    f2c_send_count=0
    DO ichunk=1,nchunk
       DO irow=1,fm_nrow_local
          ipcol=pcol_of_fm_row(irow)
          iprow=chunk_on_prow(ichunk,ipcol)
          f2c_send_count(blacs2mpi(iprow,ipcol))=&
               f2c_send_count(blacs2mpi(iprow,ipcol))+ncol_local_of_chunk(ichunk)
       END DO
    END DO
 
    CALL timeset(routineN//"_c","I","",handle1)
    CALL mp_alltoall( sb=f2c_send_count, rb=f2c_recv_count, count=1, group=fm_para_env%group )
    CALL timestop(0.0_wp,handle1)

    IF (transpose_layout) THEN
       DEALLOCATE(blacs2mpi,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE fm_reshuffle_create_layout
!***************************************************************************

!!****f* cp_sm_fm_interactions/fm2chunk [1.0] *
!!
!!   NAME
!!     fm2chunk
!!
!!   FUNCTION
!!     transfer a full matrix from block cyclic format to chunk format
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele & Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [j&f]
!!
!!*** **********************************************************************
  SUBROUTINE fm2chunk(sm,fm,chunks,shift,hole_pos,sm_transposed,&
       alpha,ncol, error)

    TYPE(real_matrix_type), POINTER :: sm
    TYPE(cp_fm_type), POINTER :: fm
    TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: chunks
    INTEGER, INTENT(in) :: shift, hole_pos,ncol
    LOGICAL,INTENT(in) :: sm_transposed
    REAL(kind=wp), INTENT(in) :: alpha
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    INTEGER :: fm_ncol_global, ichunk_in_chunks
    INTEGER :: nprow,npcol,fm_ncol_local, fm_nrow_local, irow2,bsize
    INTEGER :: nchunk,chunk_size
    INTEGER :: nchunk_local_in,nchunk_local_out, fm_nprow, fm_npcol

    CHARACTER(len=*), PARAMETER :: routineN='fm2chunk',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: myprow,mypcol, stat, handle, handle2,&
         ip,icol,irow,ipcol,iprow, ichunk,&
         ichunk_offset, ichunk_ordered, irow_chunk, iblock_atomic,&
         first_grow,last_grow, row_offset, iblock, iblock_global, col_offset,&
         chunk_in, chunk_out, iblock_row, iblock_col, sbnrow,sbncol,&
         chunk_recv_in, chunk_recv_out, row_dest, row_src, col_dest,&
         col_src, sm_nblock_row, sm_nblock_col, fm_ncol_block, fm_nrow_block,&
         fm_first_prow, fm_first_pcol
    INTEGER, DIMENSION(:), POINTER :: chunk_of_col, ncol_local_of_chunk,&
         fm_row_indices, fm_col_indices,&
         gindex_of_chunk,ncol_global_of_chunk,&
         pcol_of_fm_row, send_count, rcv_offsets, rcv_count, send_c, rcv_c,&
         send_offsets,&
         rcv_pos, local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution
    INTEGER, DIMENSION(:,:), POINTER :: chunk_on_prow, &
         blacs2mpi, ipcol_of_chunk_col, fm_blacs2mpi
    REAL(kind=wp), DIMENSION(:), POINTER :: send_buffer, rcv_buffer
    REAL(kind=wp), DIMENSION(:,:), POINTER :: sparse_block
    REAL(kind=wp) :: flops, alpha_l

    TYPE(real_block_node_type), POINTER :: block_node
    TYPE(cp_para_env_type), POINTER :: fm_para_env

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxg2p
#endif

    CALL timeset(routineN,"I","",handle)
    NULLIFY(chunk_of_col, ncol_local_of_chunk,&
         fm_row_indices, fm_col_indices,&
         gindex_of_chunk,ncol_global_of_chunk,&
         pcol_of_fm_row, send_count, send_c, rcv_c,&
         rcv_offsets, rcv_count, send_offsets,&
         rcv_pos, local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution)
    NULLIFY(chunk_on_prow, &
         blacs2mpi, ipcol_of_chunk_col,send_buffer, rcv_buffer,sparse_block,&
         block_node,fm_para_env)
    CPPrecondition(ASSOCIATED(sm),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(chunks),cp_failure_level,routineP,error,failure)
#ifndef __SCALAPACK
    CPAssert(fm%matrix_struct%para_env%num_pe==1,cp_failure_level,routineP,error,failure)
#endif

    IF (.NOT.sm_transposed) THEN
       nprow         =fm%matrix_struct%context%num_pe(1)
       npcol         =fm%matrix_struct%context%num_pe(2)
       blacs2mpi     => fm%matrix_struct%context%blacs2mpi
       
       local_atomic_rows => sm%distribution_2d%local_rows
       local_atomic_cols => sm%distribution_2d%local_cols
       atomic_row_distribution => sm%distribution_2d%row_distribution
       atomic_col_distribution => sm%distribution_2d%col_distribution
       
       CALL get_matrix_info(sm,nblock_row=sm_nblock_row,&
            nblock_col=sm_nblock_col, last_row=sm_last_row,&
            last_col=sm_last_col, first_row=sm_first_row,&
            first_col=sm_first_col)
    ELSE
       npcol         =fm%matrix_struct%context%num_pe(1)
       nprow         =fm%matrix_struct%context%num_pe(2)

       ALLOCATE(blacs2mpi(0:SIZE(fm%matrix_struct%context%blacs2mpi,2)-1,&
            0:SIZE(fm%matrix_struct%context%blacs2mpi,1)-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iprow=0,nprow-1
          DO ipcol=0,npcol-1
             blacs2mpi(iprow,ipcol)=fm%matrix_struct%context%blacs2mpi(ipcol,iprow)
          END DO
       END DO
       
       local_atomic_cols => sm%distribution_2d%local_rows
       local_atomic_rows => sm%distribution_2d%local_cols
       atomic_col_distribution => sm%distribution_2d%row_distribution
       atomic_row_distribution => sm%distribution_2d%col_distribution
       
       CALL get_matrix_info(sm,nblock_col=sm_nblock_row,&
            nblock_row=sm_nblock_col, last_col=sm_last_row,&
            last_row=sm_last_col, first_col=sm_first_row,&
            first_row=sm_first_col)
    END IF

    myprow        =fm%matrix_struct%context%mepos(1)
    mypcol        =fm%matrix_struct%context%mepos(2)

    fm_nprow=fm%matrix_struct%context%num_pe(1)
    fm_npcol=fm%matrix_struct%context%num_pe(2)
    fm_blacs2mpi => fm%matrix_struct%context%blacs2mpi
    fm_first_prow=fm%matrix_struct%first_p_pos(1)
    fm_first_pcol=fm%matrix_struct%first_p_pos(2)   
    CALL cp_fm_get_info(fm,ncol_global=fm_ncol_global,&
         ncol_local=fm_ncol_local, nrow_local=fm_nrow_local,&
         row_indices=fm_row_indices, col_indices=fm_col_indices,&
         ncol_block=fm_ncol_block, nrow_block=fm_nrow_block,&
         para_env=fm_para_env,error=error)
    CPPrecondition(fm_ncol_global>=ncol,cp_failure_level,routineP,error,failure)
    fm_ncol_global=ncol
    
    nchunk        =lcm(nprow,npcol)
    IF (nchunk==nprow.OR.nchunk==npcol)nchunk=2*nchunk
    chunk_size = (fm_ncol_global+nchunk-1)/ nchunk
    nchunk_local_in = nchunk / nprow
    nchunk_local_out= nchunk / npcol

    CPPrecondition(SIZE(chunks)>=nchunk_local_in+1,cp_failure_level,routineP,error,failure)
    DO ichunk=1,SIZE(chunks)
       CPPrecondition(ASSOCIATED(chunks(ichunk)%array),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(chunks(ichunk)%array,1)==chunk_size,cp_failure_level,routineP,error,failure)
    END DO

    CALL fm_reshuffle_create_layout(sm=sm,fm=fm, transpose_layout=sm_transposed, &
         shift=shift,ncol=fm_ncol_global,&
         gindex_of_chunk=gindex_of_chunk, ncol_global_of_chunk=ncol_global_of_chunk,&
         ipcol_of_chunk_col=ipcol_of_chunk_col,&
         chunk_of_col=chunk_of_col, ncol_local_of_chunk=ncol_local_of_chunk, &
         pcol_of_fm_row=pcol_of_fm_row, f2c_send_count=send_count,&
         f2c_recv_count=rcv_count, chunk_on_prow=chunk_on_prow, error=error)

    ALLOCATE(send_offsets(0:fm_para_env%num_pe-1),rcv_offsets(0:fm_para_env%num_pe-1),&
         send_c(0:fm_para_env%num_pe-1),rcv_c(0:fm_para_env%num_pe-1),&
         stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    send_offsets(0)=0
    rcv_offsets(0)=0
    DO ip=1,fm_para_env%num_pe-1
       send_offsets(ip)=send_offsets(ip-1)+send_count(ip-1)
       rcv_offsets(ip)=rcv_offsets(ip-1)+rcv_count(ip-1)
    END DO

    ALLOCATE(send_buffer(SUM(send_count)),&
         rcv_buffer(SUM(rcv_count)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    send_c=send_count
    IF (alpha/=1.0_wp) THEN
       DO icol=1,fm_ncol_local
          IF (fm_col_indices(icol)>fm_ncol_global) EXIT
          DO irow=1,fm_nrow_local
             ipcol=pcol_of_fm_row(irow)
             iprow=chunk_on_prow(chunk_of_col(icol),ipcol)
             ip=blacs2mpi(iprow,ipcol)
             send_offsets(ip)=send_offsets(ip)+1
             send_buffer(send_offsets(ip))=alpha*fm%local_data(irow,icol)
!FM             send_c(ip)=send_c(ip)-1
          END DO
       END DO
!FM       CPPostcondition(ALL(send_c==0),cp_failure_level,routineP,error,failure)
    ELSE
       DO icol=1,fm_ncol_local
          IF (fm_col_indices(icol)>fm_ncol_global) EXIT
          DO irow=1,fm_nrow_local
             ipcol=pcol_of_fm_row(irow)
             iprow=chunk_on_prow(chunk_of_col(icol),ipcol)
             ip=blacs2mpi(iprow,ipcol)
             send_offsets(ip)=send_offsets(ip)+1
             send_buffer(send_offsets(ip))=fm%local_data(irow,icol)
!FM             send_c(ip)=send_c(ip)-1
          END DO
       END DO
!FM       CPPostcondition(ALL(send_c==0),cp_failure_level,routineP,error,failure)
    END IF

    send_offsets(0)=0
    DO ip=1,fm_para_env%num_pe-1
       send_offsets(ip)=send_offsets(ip-1)+send_count(ip-1)
    END DO

    CALL timeset(routineN//"_all2all","I","",handle2)
    CALL mp_alltoall( sb=send_buffer, scount=send_count, sdispl=send_offsets,&
         rb=rcv_buffer, rcount=rcv_count, rdispl=rcv_offsets, &
         group=fm_para_env%group )
    CALL timestop(0.0_wp,handle2)

    ALLOCATE(rcv_pos(0:fm_para_env%num_pe-1), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    rcv_pos=rcv_offsets
    
    IF (gindex_of_chunk(1)+nchunk_local_in > nchunk) THEN
       ichunk_offset=nchunk-(gindex_of_chunk(1)-1)
    ELSE
       ichunk_offset=0
    END IF

    rcv_c=rcv_count
    DO ichunk_ordered=1,nchunk_local_in
       ichunk=MODULO(ichunk_ordered+ichunk_offset-1,nchunk_local_in)+1
       ichunk_in_chunks=MODULO(ichunk-1+hole_pos,nchunk_local_in+1)+1

       DO icol=1,ncol_global_of_chunk(gindex_of_chunk(ichunk))
          irow_chunk=0
          ipcol = ipcol_of_chunk_col(icol,ichunk)
          DO iblock_atomic=1,SIZE(local_atomic_cols)
             first_grow=sm_first_col(local_atomic_cols(iblock_atomic))
             last_grow=sm_last_col(local_atomic_cols(iblock_atomic))
             bsize=fm_nrow_block-MOD(first_grow-1,fm_nrow_block)
#ifdef __SCALAPACK
             iprow =MOD( fm_first_prow + (first_grow - 1) / fm_nrow_block, fm_nprow )
!FM                iprow = indxg2p(first_grow,fm_nrow_block, &
!FM                     myprow, fm_first_prow,fm%matrix_struct%context%num_pe(1))
#else
             iprow =0
#endif
             ip=fm_blacs2mpi(iprow,ipcol)
             DO irow=first_grow,MIN(first_grow+bsize-1,last_grow)
                irow_chunk=irow_chunk+1
!FM                rcv_c(ip)=rcv_c(ip)-1
                rcv_pos(ip)=rcv_pos(ip)+1
                chunks(ichunk_in_chunks)%array(icol,irow_chunk)=&
                     rcv_buffer(rcv_pos(ip))
             END DO
             DO irow=MIN(first_grow+bsize-1,last_grow)+1,last_grow,fm_nrow_block
                iprow=MOD(iprow+1,fm_nprow)
                ip=fm_blacs2mpi(iprow,ipcol)
                DO irow2=irow,MIN(irow+fm_nrow_block-1,last_grow)
                   irow_chunk=irow_chunk+1
!FM                   rcv_c(ip)=rcv_c(ip)-1
                   rcv_pos(ip)=rcv_pos(ip)+1
                   chunks(ichunk_in_chunks)%array(icol,irow_chunk)=&
                        rcv_buffer(rcv_pos(ip))                      
                END DO
             END DO
          END DO
       END DO
    END DO

!FM    CPPostcondition(ALL(rcv_c==0),cp_failure_level,routineP,error,failure)

    IF (sm_transposed) THEN
       DEALLOCATE(blacs2mpi,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    DEALLOCATE(rcv_buffer,rcv_offsets, rcv_count, rcv_c, rcv_pos, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    DEALLOCATE(send_buffer,send_offsets, send_count, send_c, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    DEALLOCATE(gindex_of_chunk,ncol_global_of_chunk,ipcol_of_chunk_col, &
         chunk_of_col,ncol_local_of_chunk,pcol_of_fm_row,chunk_on_prow,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE fm2chunk
!***************************************************************************


!!****f* cp_sm_fm_interactions/chunk2fm [1.0] *
!!
!!   NAME
!!     chunk2fm
!!
!!   FUNCTION
!!     transfer a full matrix from chunk format to block cyclic format
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele & Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [j&f]
!!
!!*** **********************************************************************
  SUBROUTINE chunk2fm(sm,fm,chunks,shift,hole_pos,sm_transposed,beta,ncol,error)

    TYPE(real_matrix_type), POINTER :: sm
    TYPE(cp_fm_type), POINTER :: fm
    TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: chunks
    INTEGER, INTENT(in) :: shift, hole_pos,ncol
    LOGICAL,INTENT(in) :: sm_transposed
    REAL(kind=wp), INTENT(in) :: beta
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    INTEGER :: fm_ncol_global, ichunk_in_chunks
    INTEGER :: nprow,npcol,fm_ncol_local, fm_nrow_local
    INTEGER :: nchunk,chunk_size
    INTEGER :: nchunk_local_in,nchunk_local_out

    CHARACTER(len=*), PARAMETER :: routineN='chunk2fm',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: myprow,mypcol, stat,handle, handle1,&
         ip,icol,irow,ipcol,iprow, ichunk,&
         ichunk_offset, ichunk_ordered, irow_chunk, iblock_atomic,&
         irow_global, row_offset, iblock, iblock_global, col_offset,&
         chunk_in, chunk_out, iblock_row, iblock_col, sbnrow,sbncol,&
         chunk_recv_in, chunk_recv_out, row_dest, row_src, col_dest,&
         col_src, sm_nblock_row, sm_nblock_col, fm_ncol_block, fm_nrow_block,&
         fm_first_prow, fm_first_pcol, fm_nprow, fm_npcol, first_grow, last_grow,&
         bsize, irow2
    INTEGER, DIMENSION(:), POINTER :: chunk_of_col, ncol_local_of_chunk,&
         fm_row_indices, fm_col_indices,&
         gindex_of_chunk, ncol_global_of_chunk, &
         pcol_of_fm_row, send_count, rcv_offsets,rcv_count,send_offsets,&
         send_c,rcv_c,&
         send_pos, local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution
    INTEGER, DIMENSION(:,:), POINTER :: chunk_on_prow, &
         blacs2mpi, ipcol_of_chunk_col, fm_blacs2mpi
    REAL(kind=wp), DIMENSION(:), POINTER :: send_buffer, rcv_buffer
    REAL(kind=wp), DIMENSION(:,:), POINTER :: sparse_block
    REAL(kind=wp) :: flops, alpha_l

    TYPE(real_block_node_type), POINTER :: block_node
    TYPE(cp_para_env_type), POINTER :: fm_para_env

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxg2p
#endif

    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(sm),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)
#ifndef __SCALAPACK
    CPAssert(fm%matrix_struct%para_env%num_pe==1,cp_failure_level,routineP,error,failure)
#endif
    NULLIFY(chunk_of_col, ncol_local_of_chunk,&
         fm_row_indices, fm_col_indices,&
         gindex_of_chunk, ncol_global_of_chunk, &
         pcol_of_fm_row, send_count,send_c,rcv_c,&
         rcv_offsets, rcv_count, send_offsets,&
         send_pos, local_atomic_rows, local_atomic_cols,sm_last_row,&
         sm_first_row, sm_first_col,sm_last_col, atomic_col_distribution,&
         atomic_row_distribution)
    NULLIFY(chunk_on_prow, &
         blacs2mpi, ipcol_of_chunk_col,send_buffer, rcv_buffer,sparse_block,&
         block_node,fm_para_env)

    IF (.NOT.sm_transposed) THEN
       nprow         =fm%matrix_struct%context%num_pe(1)
       npcol         =fm%matrix_struct%context%num_pe(2)
       blacs2mpi     => fm%matrix_struct%context%blacs2mpi
       
       local_atomic_rows => sm%distribution_2d%local_rows
       local_atomic_cols => sm%distribution_2d%local_cols
       atomic_row_distribution => sm%distribution_2d%row_distribution
       atomic_col_distribution => sm%distribution_2d%col_distribution
       
       CALL get_matrix_info(sm,nblock_row=sm_nblock_row,&
            nblock_col=sm_nblock_col, last_row=sm_last_row,&
            last_col=sm_last_col, first_row=sm_first_row,&
            first_col=sm_first_col)
    ELSE
       npcol         =fm%matrix_struct%context%num_pe(1)
       nprow         =fm%matrix_struct%context%num_pe(2)

       ALLOCATE(blacs2mpi(0:SIZE(fm%matrix_struct%context%blacs2mpi,2)-1,&
            0:SIZE(fm%matrix_struct%context%blacs2mpi,1)-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iprow=0,nprow-1
          DO ipcol=0,npcol-1
             blacs2mpi(iprow,ipcol)=fm%matrix_struct%context%blacs2mpi(ipcol,iprow)
          END DO
       END DO
       
       local_atomic_cols => sm%distribution_2d%local_rows
       local_atomic_rows => sm%distribution_2d%local_cols
       atomic_col_distribution => sm%distribution_2d%row_distribution
       atomic_row_distribution => sm%distribution_2d%col_distribution
       
       CALL get_matrix_info(sm,nblock_col=sm_nblock_row,&
            nblock_row=sm_nblock_col, last_col=sm_last_row,&
            last_row=sm_last_col, first_col=sm_first_row,&
            first_row=sm_first_col)
    END IF

    myprow        =fm%matrix_struct%context%mepos(1)
    mypcol        =fm%matrix_struct%context%mepos(2)

    fm_nprow=fm%matrix_struct%context%num_pe(1)
    fm_npcol=fm%matrix_struct%context%num_pe(2)
    fm_blacs2mpi => fm%matrix_struct%context%blacs2mpi
    fm_first_prow=fm%matrix_struct%first_p_pos(1)
    fm_first_pcol=fm%matrix_struct%first_p_pos(2)    
    CALL cp_fm_get_info(fm,ncol_global=fm_ncol_global,&
         ncol_local=fm_ncol_local, nrow_local=fm_nrow_local,&
         row_indices=fm_row_indices, col_indices=fm_col_indices,&
         ncol_block=fm_ncol_block, nrow_block=fm_nrow_block,&
         para_env=fm_para_env,error=error)
    CPPrecondition(fm_ncol_global>=ncol,cp_failure_level,routineP,error,failure)
    fm_ncol_global=ncol
    
    nchunk        =lcm(nprow,npcol)
    IF (nchunk==nprow.OR.nchunk==npcol)nchunk=2*nchunk
    chunk_size = (fm_ncol_global+nchunk-1)/ nchunk
    nchunk_local_in = nchunk / nprow
    nchunk_local_out= nchunk / npcol

    CALL fm_reshuffle_create_layout(sm=sm,fm=fm, transpose_layout=sm_transposed,&
         shift=shift,ncol=fm_ncol_global,&
         gindex_of_chunk=gindex_of_chunk, ncol_global_of_chunk=ncol_global_of_chunk,&
         ipcol_of_chunk_col=ipcol_of_chunk_col,&
         chunk_of_col=chunk_of_col, ncol_local_of_chunk=ncol_local_of_chunk, &
         pcol_of_fm_row=pcol_of_fm_row, f2c_send_count=rcv_count,&
         f2c_recv_count=send_count,&
         chunk_on_prow=chunk_on_prow, error=error)

    ALLOCATE(send_offsets(0:fm_para_env%num_pe-1),rcv_offsets(0:fm_para_env%num_pe-1),&
         send_c(0:fm_para_env%num_pe-1),rcv_c(0:fm_para_env%num_pe-1),&
         stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    send_offsets(0)=0
    rcv_offsets(0)=0
    DO ip=1,fm_para_env%num_pe-1
       send_offsets(ip)=send_offsets(ip-1)+send_count(ip-1)
       rcv_offsets(ip)=rcv_offsets(ip-1)+rcv_count(ip-1)
    END DO

    ALLOCATE(send_buffer(SUM(send_count)),&
         rcv_buffer(SUM(rcv_count)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(send_pos(0:fm_para_env%num_pe-1), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    send_pos=send_offsets

    IF (gindex_of_chunk(1)+nchunk_local_in > nchunk) THEN
       ichunk_offset=nchunk-(gindex_of_chunk(1)-1)
    ELSE
       ichunk_offset=0
    END IF
    send_c=send_count
    DO ichunk_ordered=1,nchunk_local_in
       ichunk=MODULO(ichunk_ordered+ichunk_offset-1,nchunk_local_in)+1
       ichunk_in_chunks=MODULO(ichunk+hole_pos-1,nchunk_local_in+1)+1

       DO icol=1,ncol_global_of_chunk(gindex_of_chunk(ichunk))
          irow_chunk=0
          ipcol = ipcol_of_chunk_col(icol,ichunk)
          DO iblock_atomic=1,SIZE(local_atomic_cols)
             first_grow=sm_first_col(local_atomic_cols(iblock_atomic))
             last_grow=sm_last_col(local_atomic_cols(iblock_atomic))
             bsize=fm_nrow_block-MOD(first_grow-1,fm_nrow_block)
#ifdef __SCALAPACK
             iprow =MOD( fm_first_prow + (first_grow - 1) / fm_nrow_block, fm_nprow )
!FM                iprow = indxg2p(first_grow,fm_nrow_block, &
!FM                     myprow, fm_first_prow,fm%matrix_struct%context%num_pe(1))
#else
             iprow =0
#endif
             ip=fm_blacs2mpi(iprow,ipcol)
             DO irow=first_grow,MIN(first_grow+bsize-1,last_grow)
                irow_chunk=irow_chunk+1
!FM                send_c(ip)=send_c(ip)-1
                send_pos(ip)=send_pos(ip)+1
                send_buffer(send_pos(ip))=&
                     chunks(ichunk_in_chunks)%array(icol,irow_chunk)
             END DO
             DO irow=MIN(first_grow+bsize-1,last_grow)+1,last_grow,fm_nrow_block
                iprow=MOD(iprow+1,fm_nprow)
                ip=fm_blacs2mpi(iprow,ipcol)
                DO irow2=irow,MIN(irow+fm_nrow_block-1,last_grow)
                irow_chunk=irow_chunk+1
!FM                send_c(ip)=send_c(ip)-1
                send_pos(ip)=send_pos(ip)+1
                send_buffer(send_pos(ip))=&
                     chunks(ichunk_in_chunks)%array(icol,irow_chunk)
                END DO
             END DO
             
          END DO
       END DO
    END DO
!FM    CPPostcondition(ALL(send_c==0),cp_failure_level,routineP,error,failure)
    CALL timeset(routineN//"_all2all","I","",handle1)
    CALL mp_alltoall( sb=send_buffer, scount=send_count, sdispl=send_offsets,&
         rb=rcv_buffer, rcount=rcv_count, rdispl=rcv_offsets, &
         group=fm_para_env%group )
    CALL timestop(0.0_wp,handle1)

    rcv_c=rcv_count
    IF (beta/=0.0_wp) THEN
       DO icol=1,fm_ncol_local
          IF (fm_col_indices(icol)>fm_ncol_global) EXIT
          DO irow=1,fm_nrow_local
             ipcol=pcol_of_fm_row(irow)
             iprow=chunk_on_prow(chunk_of_col(icol),ipcol)
             ip=blacs2mpi(iprow,ipcol)
             rcv_offsets(ip)=rcv_offsets(ip)+1
             fm%local_data(irow,icol)=beta*fm%local_data(irow,icol)+rcv_buffer(rcv_offsets(ip))
             rcv_c(ip)=rcv_c(ip)-1
          END DO
       END DO
       CPPostcondition(ALL(rcv_c==0),cp_failure_level,routineP,error,failure)
    ELSE
       DO icol=1,fm_ncol_local
          IF (fm_col_indices(icol)>fm_ncol_global) EXIT
          DO irow=1,fm_nrow_local
             ipcol=pcol_of_fm_row(irow)
             iprow=chunk_on_prow(chunk_of_col(icol),ipcol)
             ip=blacs2mpi(iprow,ipcol)
             rcv_offsets(ip)=rcv_offsets(ip)+1
             fm%local_data(irow,icol)=rcv_buffer(rcv_offsets(ip))
             rcv_c(ip)=rcv_c(ip)-1
          END DO
       END DO
       CPPostcondition(ALL(rcv_c==0),cp_failure_level,routineP,error,failure)
    END IF

    IF (sm_transposed) THEN
       DEALLOCATE(blacs2mpi,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    DEALLOCATE(rcv_buffer,rcv_offsets, rcv_count, rcv_c, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    DEALLOCATE(send_buffer,send_offsets, send_count,send_c,send_pos, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    DEALLOCATE(chunk_on_prow,ipcol_of_chunk_col,chunk_of_col,ncol_local_of_chunk, &
        gindex_of_chunk,ncol_global_of_chunk,pcol_of_fm_row, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE chunk2fm
!***************************************************************************

!!****f* cp_sm_fm_interactions/cp_sm_fm_multiply_tests [1.0] *
!!
!!   NAME
!!     cp_sm_fm_multiply_tests
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_fm_multiply_tests(blacs_env, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_fm_multiply_tests
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_fm_multiply_tests(blacs_env,error)
  TYPE(cp_blacs_env_type), POINTER :: blacs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_sm_fm_multiply_tests',&
        routineP=moduleN//':'//routineN
  TYPE(cp_fm_struct_type), POINTER :: fm_struct
  TYPE(cp_fm_type), POINTER :: fm1, fm2, fm3
  TYPE(real_matrix_type), POINTER :: sm
  TYPE(distribution_2d_type), POINTER :: atomic_distribution
  INTEGER, PARAMETER :: ncols=7
  INTEGER, DIMENSION(3), PARAMETER :: block_sizes=(/2,1,1/),&
       row_dist=(/1,0,0/), col_dist=(/1,0,1/)
  INTEGER :: nrows, nblocks, unit_nr,iblock_atomic, stat,&
       irow,icol,ichunk,iblock, iblock_row, iblock_col, sm_ncol_local,&
       nchunk, chunk_size, nchunk_local_in, nchunk_local_out
  INTEGER, DIMENSION(:), POINTER :: first_row,first_col, last_row, last_col
  TYPE(cp_logger_type), POINTER :: logger
  TYPE(cp_2d_r_p_type), DIMENSION(:),POINTER :: v_in
  REAL(kind=wp) :: diff
  REAL(kind=wp), DIMENSION(:,:), POINTER :: sm_block

  logger => cp_error_get_logger(error)
  failure=.FALSE.
  nrows=SUM(block_sizes)
  nblocks=SIZE(block_sizes)
  unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
  nchunk        =lcm(blacs_env%num_pe(1),blacs_env%num_pe(2))
  IF (nchunk==blacs_env%num_pe(1).OR.nchunk==blacs_env%num_pe(2))nchunk=2*nchunk
  chunk_size = (ncols+nchunk-1)/ nchunk
  nchunk_local_in = nchunk / blacs_env%num_pe(1)
  nchunk_local_out= nchunk / blacs_env%num_pe(2)

  WRITE(unit=unit_nr,fmt=*) "nrows=",nrows,"ncols=",ncols,"nblocks=",nblocks
  WRITE(unit=unit_nr,fmt=*) "block_sizes=",block_sizes
  WRITE(unit=unit_nr,fmt=*) "row_dist=",row_dist
  WRITE(unit=unit_nr,fmt=*) "col_dist=",col_dist
  CALL m_flush(unit_nr)

  CALL cp_fm_struct_create(fm_struct,context=blacs_env, nrow_global=nrows,&
     ncol_global=ncols,error=error)
  CALL cp_fm_struct_write(fm_struct,unit_nr=unit_nr,long_description=.TRUE.,&
       error=error)
  CALL cp_fm_create(fm1,fm_struct,name="fm1",error=error)
  CALL cp_fm_create(fm2,fm_struct,name="fm2",error=error)
  CALL cp_fm_create(fm3,fm_struct,name="fm3",error=error)
  CALL cp_fm_struct_release(fm_struct,error=error)

  ALLOCATE(first_row(nblocks),first_col(nblocks),last_row(nblocks),&
       last_col(nblocks),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  first_row(1)=1
  first_col(1)=1
  last_row(1)=block_sizes(1)
  last_col(1)=last_row(1)
  DO iblock=2,nblocks
     first_row(iblock)=first_row(iblock-1)+block_sizes(iblock-1)
     first_col(iblock)=first_row(iblock)
     last_row(iblock)=last_row(iblock-1)+block_sizes(iblock)
     last_col(iblock)=last_row(iblock)
  END DO
  
  CALL distribution_2d_create(atomic_distribution, row_distribution=row_dist,&
       col_distribution=col_dist, blacs_env=blacs_env,error=error)
  
  CALL allocate_matrix(sm,nrow=nrows,ncol=nrows,nblock_row=nblocks,&
       nblock_col=nblocks,first_row=first_row,last_row=last_row,&
       first_col=first_col,last_col=last_col,&
       matrix_name="aria",matrix_symmetry="none",&
       distribution_2d=atomic_distribution)

  DO iblock_row=1,SIZE(atomic_distribution%local_rows)
     DO iblock_col=1,SIZE(atomic_distribution%local_cols)
!FM        IF (atomic_distribution%local_cols(iblock_col)<=&
!FM             atomic_distribution%local_rows(iblock_row)) THEN
           NULLIFY(sm_block)
           IF (((MODULO(atomic_distribution%local_cols(iblock_col)+&
             atomic_distribution%local_rows(iblock_row),2)==0).AND.&
             (atomic_distribution%local_cols(iblock_col)<=&
             atomic_distribution%local_rows(iblock_row))).OR.&
              ((MODULO(atomic_distribution%local_cols(iblock_col)+&
             atomic_distribution%local_rows(iblock_row),2)/=0).AND.&
             (atomic_distribution%local_cols(iblock_col)>&
             atomic_distribution%local_rows(iblock_row)))) THEN
              CALL add_block_node(matrix=sm,&
                   block_row=atomic_distribution%local_rows(iblock_row),&
                   block_col=atomic_distribution%local_cols(iblock_col),&
                   BLOCK=sm_block)
!FM           ELSE
!FM              CALL add_block_node(matrix=sm,&
!FM                   block_row=atomic_distribution%local_cols(iblock_col),&
!FM                   block_col=atomic_distribution%local_rows(iblock_row),&
!FM                   BLOCK=sm_block)
!FM           END IF
           sm_block=0.0_wp
           DO irow=1,MIN(SIZE(sm_block,1),SIZE(sm_block,2))
              sm_block(irow,irow)=atomic_distribution%local_rows(iblock_row)*&
                   atomic_distribution%local_cols(iblock_col)
           END DO
        END IF
     END DO
  END DO

  CALL distribution_2d_release(atomic_distribution,error)

  DO icol=1,ncols
     DO irow=1,nrows
        CALL cp_fm_set_element(fm1,irow_global=irow,icol_global=icol,&
             alpha=REAL(icol+(irow-1)*ncols,wp),error=error)
     END DO
  END DO

!FM  sm_ncol_local=0
!FM  DO iblock=1,SIZE(atomic_distribution%local_cols)
!FM     sm_ncol_local=sm_ncol_local+last_col(atomic_distribution%local_cols(iblock))-&
!FM          first_col(atomic_distribution%local_cols(iblock))+1
!FM  END DO
!FM  ! alloc & init v_in
!FM  ALLOCATE(v_in (nchunk_local_in+1),stat=stat)
!FM  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM  DO ichunk=1,nchunk_local_in+1
!FM     ALLOCATE(v_in(ichunk)%array(chunk_size,sm_ncol_local), stat=stat)
!FM     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM  END DO
!FM
!FM  CALL fm2chunk(sm,fm1,v_in,shift=0,hole_pos=1,sm_transposed=.FALSE.,&
!FM       alpha=1.0_wp,error=error)
!FM  
!FM  WRITE(unit=unit_nr,fmt=*) "---- local chunks ----"
!FM  DO ichunk=1,SIZE(v_in)
!FM     WRITE(unit=unit_nr,fmt=*) "local chunk ",ichunk
!FM     CALL cp_2d_r_write(v_in(ichunk)%array,unit_nr=unit_nr,error=error)
!FM  END DO
!FM  WRITE(unit=unit_nr,fmt=*) "---- shuffle back ----"
!FM  CALL m_flush(unit_nr)
!FM
!FM  CALL chunk2fm(sm,fm2,v_in,shift=0,hole_pos=1,sm_transposed=.FALSE.,&
!FM       beta=0.0_wp,error=error)
!FM
!FM  CALL cp_fm_scale_and_add(alpha=1.0_wp,matrix_a=fm2,beta=-1.0_wp,&
!FM       matrix_b=fm1,error=error)
!FM  
!FM  CALL cp_fm_maxval(fm1,diff)
!FM  
!FM  WRITE(unit=unit_nr,fmt=*) "diff=",diff
!FM  CALL m_flush(unit_nr)
!FM  
!FM  IF (diff/=0.0_wp) THEN
!FM     CALL cp_fm_write(fm2,unit_nr=unit_nr, long_description=.TRUE., &
!FM          local=.TRUE., error=error)
!FM     CALL m_flush(unit_nr)
!FM     CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
!FM  END IF

  WRITE(unit=unit_nr,fmt=*) " * multiply! *"
  CALL m_flush(unit_nr)

  CALL cp_sm_fm_multiply_2d(sm=sm,fm_in=fm1,fm_out=fm2,transpose_sm=.FALSE.,&
       do_diagonal_blocks=.TRUE.,antisymmetric=.FALSE.,alpha=1.0_wp,beta=0.0_wp,ncol=ncols-3,error=error)
!FM  WRITE(unit=unit_nr,fmt=*) " * again *"
!FM  CALL m_flush(unit_nr)
!FM  CALL cp_sm_fm_multiply_2d(sm=sm,fm_in=fm1,fm_out=fm2,transpose_sm=.TRUE.,&
!FM       do_diagonal_blocks=.FALSE.,antisymmetric=.FALSE.,alpha=1.0_wp,beta=1.0_wp,error=error)
  CALL cp_sm_fm_multiply(sparse_matrix=sm,v_in=fm1,v_out=fm3,ncol=ncols-3,&
       para_env=fm1%matrix_struct%para_env,&
       alpha=1.0_wp,beta=0.0_wp,error=error)
  CALL cp_fm_scale_and_add(alpha=1.0_wp,matrix_a=fm3,beta=-1.0_wp,&
       matrix_b=fm2,error=error)
  
  CALL m_flush(unit_nr)
  CALL cp_fm_maxval(fm3,diff)
  
  WRITE(unit=unit_nr,fmt=*) "diff=",diff
  CALL m_flush(unit_nr)
  
  IF (diff>1.0e-15_wp) THEN
     WRITE(unit=unit_nr,fmt=*) "new_mult"
     CALL cp_fm_write(fm2,unit_nr=unit_nr, long_description=.TRUE., &
          local=.TRUE., error=error)
     CALL m_flush(unit_nr)
     WRITE(unit=unit_nr,fmt=*) "diff_mult"
     CALL cp_fm_write(fm3,unit_nr=unit_nr, long_description=.TRUE., &
          local=.TRUE., error=error)
     CALL m_flush(unit_nr)
     CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
  END IF

  STOP "OK"

  CALL cp_fm_release(fm1,error=error)
  CALL cp_fm_release(fm2,error=error)
  CALL deallocate_matrix(sm)
END SUBROUTINE cp_sm_fm_multiply_tests
!***************************************************************************

END MODULE cp_sm_fm_interactions

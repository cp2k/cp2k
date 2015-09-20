!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines for calculating RI-MP2 gradients
!> \par History
!>      10.2013 created [Mauro Del Ben]
! *****************************************************************************
MODULE mp2_ri_grad_util
! 
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
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
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_irecv,&
                                             mp_isend,&
                                             mp_sendrecv,&
                                             mp_sum,&
                                             mp_wait,&
                                             mp_waitall
  USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                             mp2_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_grad_util'

  PUBLIC :: complete_gamma

  TYPE index_map
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: map
  END TYPE

  CONTAINS

! *****************************************************************************
!> \brief complete the calculation of the Gamma matrices 
!> \param mp2_env ...
!> \param B_ia_Q ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ngroup ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_B_size ...
!> \param my_B_virtual_start ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param sub_proc_map ...
!> \param alpha_case ...
!> \author Mauro Del Ben 
! *****************************************************************************
  SUBROUTINE complete_gamma(mp2_env,B_ia_Q,dimen_RI,homo,virtual,para_env,para_env_sub,ngroup,&
                            my_group_L_size,my_group_L_start,my_group_L_end,&
                            my_B_size,my_B_virtual_start,&
                            ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                            sub_proc_map,alpha_case)

    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: B_ia_Q
    INTEGER                                  :: dimen_RI, homo, virtual
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER :: ngroup, my_group_L_size, my_group_L_start, my_group_L_end, &
      my_B_size, my_B_virtual_start
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
      sub_proc_map
    LOGICAL                                  :: alpha_case

    CHARACTER(LEN=*), PARAMETER :: routineN = 'complete_gamma', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_offset, col_size, dimen_ia, handle, i, i_global, &
      ia_global, iiB, iproc, itmp(2), j_global, jjB, kkB, my_ia_end, &
      my_ia_size, my_ia_start, my_P_end, my_P_size, my_P_start, ncol_local, &
      nrow_local, pos_group, pos_sub, proc_receive, proc_send, proc_shift, &
      rec_B_size, rec_B_virtual_end, rec_B_virtual_start, row, row_offset, &
      row_size, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_ia, ends_P, pos_info, &
                                                proc_map, sizes_ia, sizes_P, &
                                                starts_ia, starts_P
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ends, group_grid_2_mepos, &
                                                mepos_2_grid_group, sizes, &
                                                starts
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_2D, Gamma_2D
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C_rec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_Gamma, fm_Gamma_PQ, &
                                                fm_Gamma_PQ_temp, fm_ia_P, &
                                                fm_Y, PQ_half

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! reshape the data into a 2D array
    ! reorder the local data in such a way to help the next stage of matrix creation
    ! now the data inside the group are divided into a ia x K matrix
    dimen_ia=homo*virtual
    ALLOCATE(sizes_ia(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    sizes_ia=0
    ALLOCATE(starts_ia(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    starts_ia=0
    ALLOCATE(ends_ia(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    ends_ia=0

    DO iproc=0, para_env_sub%num_pe-1
      itmp=get_limit(dimen_ia,para_env_sub%num_pe,iproc)
      starts_ia(iproc)=itmp(1)
      ends_ia(iproc)=itmp(2)
      sizes_ia(iproc)=itmp(2)-itmp(1)+1
    END DO

    my_ia_start=starts_ia(para_env_sub%mepos)
    my_ia_end=ends_ia(para_env_sub%mepos)
    my_ia_size=sizes_ia(para_env_sub%mepos)

    ! reorder data (first (ia|K))
    ALLOCATE(BIb_C_2D(my_ia_size,my_group_L_size),STAT=stat)
    CPASSERT(stat==0)
    BIb_C_2D=0.0_dp

    DO iiB=1, homo
      DO jjB=1, my_B_size
        ia_global=(iiB-1)*virtual+my_B_virtual_start+jjB-1
        IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
          BIb_C_2D(ia_global-my_ia_start+1,1:my_group_L_size)=B_ia_Q(iiB,jjB,1:my_group_L_size)
        END IF
      END DO
    END DO

    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

      rec_B_size=sizes_B_virtual(proc_receive)
      rec_B_virtual_end=ends_B_virtual(proc_receive)
      rec_B_virtual_start=starts_B_virtual(proc_receive)

      ALLOCATE(BIb_C_rec(homo,rec_B_size,my_group_L_size),STAT=stat)
      CPASSERT(stat==0)
      BIb_C_rec=0.0_dp

      CALL  mp_sendrecv(B_ia_Q,proc_send,&
                        BIb_C_rec,proc_receive,&
                        para_env_sub%group)

      DO iiB=1, homo
        DO jjB=1, rec_B_size
          ia_global=(iiB-1)*virtual+rec_B_virtual_start+jjB-1
          IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
            BIb_C_2D(ia_global-my_ia_start+1,1:my_group_L_size)=BIb_C_rec(iiB,jjB,1:my_group_L_size)
          END IF
        END DO
      END DO

      DEALLOCATE(BIb_C_rec)
    END DO
    DEALLOCATE(B_ia_Q)

    ! reorder data (Gamma)
    ALLOCATE(Gamma_2D(my_ia_size,my_group_L_size),STAT=stat)
    CPASSERT(stat==0)
    Gamma_2D=0.0_dp

    DO iiB=1, homo
      DO jjB=1, my_B_size
        ia_global=(iiB-1)*virtual+my_B_virtual_start+jjB-1
        IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
          ! Closed-shell, open-shell alpha-alpha component
          IF (alpha_case) THEN 
             Gamma_2D(ia_global-my_ia_start+1,1:my_group_L_size)=&
             mp2_env%ri_grad%Gamma_P_ia(iiB,jjB,1:my_group_L_size)
          ELSE ! Open-shell beta-beta component
             Gamma_2D(ia_global-my_ia_start+1,1:my_group_L_size)=&
             mp2_env%ri_grad%Gamma_P_ia_beta(iiB,jjB,1:my_group_L_size)
          ENDIF
        END IF
      END DO
    END DO

    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

      rec_B_size=sizes_B_virtual(proc_receive)
      rec_B_virtual_end=ends_B_virtual(proc_receive)
      rec_B_virtual_start=starts_B_virtual(proc_receive)

      ALLOCATE(BIb_C_rec(homo,rec_B_size,my_group_L_size),STAT=stat)
      CPASSERT(stat==0)
      BIb_C_rec=0.0_dp

      IF (alpha_case) THEN 
         CALL  mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia,proc_send,&
                        BIb_C_rec,proc_receive,&
                        para_env_sub%group)
      ELSE
         CALL  mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia_beta,proc_send,&
                        BIb_C_rec,proc_receive,&
                        para_env_sub%group)
      ENDIF

      DO iiB=1, homo
        DO jjB=1, rec_B_size
          ia_global=(iiB-1)*virtual+rec_B_virtual_start+jjB-1
          IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
            Gamma_2D(ia_global-my_ia_start+1,1:my_group_L_size)=BIb_C_rec(iiB,jjB,1:my_group_L_size)
          END IF
        END DO
      END DO

      DEALLOCATE(BIb_C_rec)
    END DO
    IF (alpha_case) THEN 
       DEALLOCATE(mp2_env%ri_grad%Gamma_P_ia)
    ELSE
       DEALLOCATE(mp2_env%ri_grad%Gamma_P_ia_beta)
    ENDIF

    ! create the processor map and size arrays
    ALLOCATE(proc_map(-para_env%num_pe:2*para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    proc_map=0
    ALLOCATE(sizes(2,0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    sizes=0
    ALLOCATE(starts(2,0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    starts=0
    ALLOCATE(ends(2,0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    ends=0

    ALLOCATE(pos_info(0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    pos_info=0
    pos_info(para_env%mepos)=para_env_sub%mepos
    CALL mp_sum(pos_info,para_env%group)

    ALLOCATE(group_grid_2_mepos(0:para_env_sub%num_pe-1,0:ngroup-1),STAT=stat)
    CPASSERT(stat==0)
    group_grid_2_mepos=0
    ALLOCATE(mepos_2_grid_group(0:para_env%num_pe-1,2),STAT=stat)
    CPASSERT(stat==0)
    mepos_2_grid_group=0

    DO i=0,para_env%num_pe-1
      proc_map(i)=i
      proc_map(-i-1)=para_env%num_pe-i-1
      proc_map(para_env%num_pe+i)=i
      ! calculate postition of the group
      pos_group=i/para_env_sub%num_pe
      ! calculate postition in the subgroup
      pos_sub=pos_info(i)
      ! 1 -> rows, 2 -> cols
      ! fill sizes array
      sizes(1,i)=sizes_ia(pos_sub)
      sizes(2,i)=sizes_array(pos_group)
      ! fill starts array
      starts(1,i)=starts_ia(pos_sub)
      starts(2,i)=starts_array(pos_group)
      ! fill ends array
      ends(1,i)=ends_ia(pos_sub)
      ends(2,i)=ends_array(pos_group)
      ! fill the map from the grid of groups to process
      group_grid_2_mepos(pos_sub,pos_group)=i
      ! and the opposite, from the global pos to the grid pos
      mepos_2_grid_group(i,1)=pos_sub
      mepos_2_grid_group(i,2)=pos_group
    END DO

    ! WRITE(*,*) mepos_2_grid_group(:,1)
    ! WRITE(*,*) mepos_2_grid_group(:,2)

    ! create the blacs env
    NULLIFY(blacs_env)
    CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

    ! create the fm matrix Gamma
    CALL array2fm(Gamma_2D,dimen_ia,dimen_RI,para_env,blacs_env,proc_map,&
                  my_ia_start,my_ia_end,&
                  my_group_L_start,my_group_L_end,&
                  sizes,starts,ends,&
                  group_grid_2_mepos,para_env_sub%num_pe,ngroup,&
                  fm_Y)
    ! create the fm matrix B_ia_P
    CALL array2fm(BIb_C_2D,dimen_ia,dimen_RI,para_env,blacs_env,proc_map,&
                  my_ia_start,my_ia_end,&
                  my_group_L_start,my_group_L_end,&
                  sizes,starts,ends,&
                  group_grid_2_mepos,para_env_sub%num_pe,ngroup,&
                  fm_ia_P)

    ! create PQ_half
    ! get the ranges and info for the (P|Q)^(-1/2) matrix
    ALLOCATE(sizes_P(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    sizes_P=0
    ALLOCATE(starts_P(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    starts_P=0
    ALLOCATE(ends_P(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    ends_P=0
    DO iproc=0, para_env_sub%num_pe-1
      itmp=get_limit(dimen_RI,para_env_sub%num_pe,iproc)
      starts_P(iproc)=itmp(1)
      ends_P(iproc)=itmp(2)
      sizes_P(iproc)=itmp(2)-itmp(1)+1
    END DO
    my_P_start=starts_P(para_env_sub%mepos)
    my_P_end=ends_P(para_env_sub%mepos)
    my_P_size=sizes_P(para_env_sub%mepos) 
    sizes=0
    starts=0
    ends=0
    DO i=0,para_env%num_pe-1
      ! calculate postition of the group
      pos_group=i/para_env_sub%num_pe
      ! calculate postition in the subgroup
      pos_sub=pos_info(i)
      ! 1 -> rows, 2 -> cols
      ! fill sizes array
      sizes(1,i)=sizes_P(pos_sub)
      sizes(2,i)=sizes_array(pos_group)
      ! fill starts array
      starts(1,i)=starts_P(pos_sub)
      starts(2,i)=starts_array(pos_group)
      ! fill ends array
      ends(1,i)=ends_P(pos_sub)
      ends(2,i)=ends_array(pos_group)
    END DO
    ! since we will need (P|Q)^(-1/2) in the future, make a copy
    IF (alpha_case) THEN
       ALLOCATE(mp2_env%ri_grad%Gamma_PQ(my_P_size,my_group_L_size),STAT=stat)
       CPASSERT(stat==0)
       mp2_env%ri_grad%Gamma_PQ(:,:)=mp2_env%ri_grad%PQ_half
       CALL array2fm(mp2_env%ri_grad%Gamma_PQ,dimen_RI,dimen_RI,para_env,blacs_env,proc_map,&
            my_P_start,my_P_end,&
            my_group_L_start,my_group_L_end,&
            sizes,starts,ends,&
            group_grid_2_mepos,para_env_sub%num_pe,ngroup,&
            PQ_half)
    ELSE
       ALLOCATE(mp2_env%ri_grad%Gamma_PQ_beta(my_P_size,my_group_L_size),STAT=stat)
       CPASSERT(stat==0)
       mp2_env%ri_grad%Gamma_PQ_beta(:,:)=mp2_env%ri_grad%PQ_half
       CALL array2fm(mp2_env%ri_grad%Gamma_PQ_beta,dimen_RI,dimen_RI,para_env,blacs_env,proc_map,&
            my_P_start,my_P_end,&
            my_group_L_start,my_group_L_end,&
            sizes,starts,ends,&
            group_grid_2_mepos,para_env_sub%num_pe,ngroup,&
            PQ_half)
    ENDIF


    ! complete the gamma matrix Gamma_ia^P = sum_Q (Y_ia^Q * (Q|P)^(-1/2))
    NULLIFY(fm_Gamma)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_ia,&
                             ncol_global=dimen_RI,para_env=para_env)
    ! for now we don't force the blocks
    ! CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=num_rows,&
    !                          ncol_global=num_cols,para_env=para_env,&
    !                          nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)
    CALL cp_fm_create(fm_Gamma,fm_struct,name="fm_Gamma")
    CALL cp_fm_set_all(matrix=fm_Gamma,alpha=0.0_dp)
    CALL cp_fm_struct_release(fm_struct)
    ! perform the matrix multiplication
    CALL cp_gemm(transa="N",transb="T", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=1.0_dp,&
                    matrix_a=fm_Y, matrix_b=PQ_half, beta=0.0_dp,&
                    matrix_c=fm_Gamma)
    ! release the Y matrix
    CALL cp_fm_release(fm_Y)
 
    ! complete gamma small (fm_Gamma_PQ)
    NULLIFY(fm_Gamma_PQ_temp)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_RI,&
                             ncol_global=dimen_RI,para_env=para_env)
    ! create temp matrix
    CALL cp_fm_create(fm_Gamma_PQ_temp,fm_struct,name="fm_Gamma_PQ_temp")
    CALL cp_fm_set_all(matrix=fm_Gamma_PQ_temp,alpha=0.0_dp)
    CALL cp_gemm(transa="T",transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia, alpha=1.0_dp,&
                    matrix_a=fm_Gamma, matrix_b=fm_ia_P, beta=0.0_dp,&
                    matrix_c=fm_Gamma_PQ_temp)
    CALL cp_fm_release(fm_ia_P)
    ! create fm_Gamma_PQ matrix
    CALL cp_fm_create(fm_Gamma_PQ,fm_struct,name="fm_Gamma_PQ") 
    ! release the structure
    CALL cp_fm_struct_release(fm_struct)
    CALL cp_fm_set_all(matrix=fm_Gamma_PQ,alpha=0.0_dp)
    ! perfome matrix multiplication
    CALL cp_gemm(transa="N",transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp,&
                    matrix_a=fm_Gamma_PQ_temp, matrix_b=PQ_half, beta=0.0_dp,&
                    matrix_c=fm_Gamma_PQ)
    CALL cp_fm_release(fm_Gamma_PQ_temp)
    CALL cp_fm_release(PQ_half)

    IF(.FALSE.) THEN
      ALLOCATE(Gamma_2D(dimen_ia,dimen_RI))
      Gamma_2D=0.0_dp
      CALL cp_fm_get_info(matrix=fm_Gamma,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          Gamma_2D(i_global,j_global)=fm_Gamma%local_data(iiB,jjB)
        END DO
      END DO
      CALL mp_sum(Gamma_2D,para_env%group)
      IF(para_env%mepos==0) CALL write_array(Gamma_2D(1:dimen_ia,1:dimen_RI))
      DEALLOCATE(Gamma_2D)
      
      ALLOCATE(Gamma_2D(dimen_RI,dimen_RI))
      Gamma_2D=0.0_dp
      CALL cp_fm_get_info(matrix=fm_Gamma_PQ,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          Gamma_2D(i_global,j_global)=fm_Gamma_PQ%local_data(iiB,jjB)
        END DO
      END DO
      CALL mp_sum(Gamma_2D,para_env%group)
      IF(para_env%mepos==0) CALL write_array(Gamma_2D(1:dimen_RI,1:dimen_RI))
      DEALLOCATE(Gamma_2D)
    END IF

    ! restore original sizes information (dimen_ia x dimen_RI)
    ! DO i=0,para_env%num_pe-1
    !   ! calculate postition of the group
    !   pos_group=i/para_env_sub%num_pe
    !   ! calculate postition in the subgroup
    !   pos_sub=pos_info(i)
    !   ! fill the map from the grid of groups to process
    !   group_grid_2_mepos(pos_sub,pos_group)=i
    ! END DO
    CALL cp_fm_get_info(matrix=fm_Gamma,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local)
    sizes=0
    sizes(1,para_env%mepos)=nrow_local
    sizes(2,para_env%mepos)=ncol_local
    CALL mp_sum(sizes,para_env%group)
    
    ! now redistribute the data, in this case we go back
    ! to the original way the integrals were distributed
    CALL fm2array(Gamma_2D,dimen_ia,dimen_RI,para_env,proc_map,&
                  my_ia_size,my_ia_start,my_ia_end,&
                  my_group_L_size,my_group_L_start,my_group_L_end,&
                  sizes,group_grid_2_mepos,mepos_2_grid_group,&
                  para_env_sub%num_pe,ngroup,&
                  fm_Gamma)

    ! the same for the Gamma_PQ matrix
    CALL cp_fm_get_info(matrix=fm_Gamma_PQ,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local)
    sizes=0
    sizes(1,para_env%mepos)=nrow_local
    sizes(2,para_env%mepos)=ncol_local
    CALL mp_sum(sizes,para_env%group)
    IF (alpha_case) THEN
       CALL fm2array(mp2_env%ri_grad%Gamma_PQ,dimen_RI,dimen_RI,para_env,proc_map,&
                  my_P_size,my_P_start,my_P_end,&
                  my_group_L_size,my_group_L_start,my_group_L_end,&
                  sizes,group_grid_2_mepos,mepos_2_grid_group,&
                  para_env_sub%num_pe,ngroup,&
                  fm_Gamma_PQ)
    ELSE
       CALL fm2array(mp2_env%ri_grad%Gamma_PQ_beta,dimen_RI,dimen_RI,para_env,proc_map,&
                  my_P_size,my_P_start,my_P_end,&
                  my_group_L_size,my_group_L_start,my_group_L_end,&
                  sizes,group_grid_2_mepos,mepos_2_grid_group,&
                  para_env_sub%num_pe,ngroup,&
                  fm_Gamma_PQ)
    ENDIF

    IF(.FALSE.) THEN
      ALLOCATE(BIb_C_2D(dimen_ia,dimen_RI))
      BIb_C_2D=0.0_dp
      BIb_C_2D(my_ia_start:my_ia_end,my_group_L_start:my_group_L_end)=Gamma_2D
      CALL mp_sum(BIb_C_2D,para_env%group)
      IF(para_env%mepos==0) CALL write_array(BIb_C_2D(1:dimen_ia,1:dimen_RI))
      DEALLOCATE(BIb_C_2D)

      ALLOCATE(BIb_C_2D(dimen_RI,dimen_RI))
      BIb_C_2D=0.0_dp
      IF (alpha_case) THEN
         BIb_C_2D(my_P_start:my_P_end,my_group_L_start:my_group_L_end)= &
         mp2_env%ri_grad%Gamma_PQ
      ELSE
         BIb_C_2D(my_P_start:my_P_end,my_group_L_start:my_group_L_end)= &
         mp2_env%ri_grad%Gamma_PQ_beta
      ENDIF
      CALL mp_sum(BIb_C_2D,para_env%group)
      IF(para_env%mepos==0) THEN 
        WRITE(*,*) 'AAAAAAAAAAAAAAAAAA'
        WRITE(*,*) 
        CALL write_array(BIb_C_2D(1:dimen_RI,1:dimen_RI))
      END IF
      DEALLOCATE(BIb_C_2D)
    END IF

    ! create the Gamma_ia_P in DBCSR style
    CALL create_dbcsr_gamma(Gamma_2D,homo,virtual,dimen_ia,para_env_sub,sub_proc_map,&
                            my_ia_start,my_ia_end,my_group_L_size,&
                            sizes_ia,starts_ia,ends_ia,&
                            mp2_env,alpha_case)


    IF(.FALSE.) THEN
      ALLOCATE(Gamma_2D(dimen_ia,dimen_RI))
      ALLOCATE(BIb_C_rec(virtual,homo,dimen_RI))
      Gamma_2D=0.0_dp
      BIb_C_rec=0.0_dp
      DO kkB=1, my_group_L_size
        IF (alpha_case) THEN
           CALL cp_dbcsr_iterator_start(iter,mp2_env%ri_grad%G_P_ia(kkB)%matrix)
        ELSE
           CALL cp_dbcsr_iterator_start(iter,mp2_env%ri_grad%G_P_ia_beta(kkB)%matrix)
        ENDIF
        DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
           CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk,&
                row_size=row_size, col_size=col_size, &
                row_offset=row_offset, col_offset=col_offset)
           DO iiB=1,row_size
             i_global=row_offset+iiB-1
             DO jjB=1,col_size
               j_global=col_offset+jjB-1
               ia_global=(i_global-1)*virtual+j_global
               Gamma_2D(ia_global,my_group_L_start+kkB-1)=data_block(iiB,jjB)
               BIb_C_rec(j_global,i_global,my_group_L_start+kkB-1)=data_block(iiB,jjB)
             ENDDO
           ENDDO
        ENDDO
        CALL cp_dbcsr_iterator_stop(iter)
      END DO
      CALL mp_sum(Gamma_2D,para_env%group)
      CALL mp_sum(BIb_C_rec,para_env%group)
      ! IF(para_env%mepos==0) CALL write_array(Gamma_2D(1:dimen_ia,1:dimen_RI))
      IF(para_env%mepos==0) THEN
        DO iiB=1, dimen_RI
          WRITE(*,*) 'AAAAAAAAAA', iiB
          CALL write_array(BIb_C_rec(1:virtual,1:homo,iiB))
        END DO
      END IF
      DEALLOCATE(Gamma_2D)
      DEALLOCATE(BIb_C_rec)
    END IF

    DEALLOCATE(pos_info)
    DEALLOCATE(group_grid_2_mepos, mepos_2_grid_group)
    DEALLOCATE(sizes_ia, starts_ia, ends_ia)
    DEALLOCATE(sizes_P, starts_P, ends_P)
    DEALLOCATE(proc_map)
    DEALLOCATE(sizes,starts,ends)

    ! release fm_matrices
    ! CALL cp_fm_release(fm_Gamma)
    ! CALL cp_fm_release(fm_Gamma_PQ)

    ! release blacs_env
    CALL cp_blacs_env_release(blacs_env)

    CALL timestop(handle)

  END SUBROUTINE complete_gamma

! *****************************************************************************
!> \brief ...
!> \param mat2D ...
!> \param num_rows ...
!> \param num_cols ...
!> \param para_env ...
!> \param blacs_env ...
!> \param proc_map ...
!> \param my_start_row ...
!> \param my_end_row ...
!> \param my_start_col ...
!> \param my_end_col ...
!> \param sizes ...
!> \param starts ...
!> \param ends ...
!> \param group_grid_2_mepos ...
!> \param ngroup_row ...
!> \param ngroup_col ...
!> \param fm_mat ...
! *****************************************************************************
  SUBROUTINE array2fm(mat2D,num_rows,num_cols,para_env,blacs_env,proc_map,&
                      my_start_row,my_end_row,&
                      my_start_col,my_end_col,&
                      sizes,starts,ends,&
                      group_grid_2_mepos,ngroup_row,ngroup_col,&
                      fm_mat)

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mat2D
    INTEGER                                  :: num_rows, num_cols
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_map
    INTEGER                                  :: my_start_row, my_end_row, &
                                                my_start_col, my_end_col
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: sizes, starts, ends, &
                                                group_grid_2_mepos
    INTEGER                                  :: ngroup_row, ngroup_col
    TYPE(cp_fm_type), POINTER                :: fm_mat

    CHARACTER(LEN=*), PARAMETER :: routineN = 'array2fm', &
      routineP = moduleN//':'//routineN

    INTEGER :: dummy_proc, end_col_block, end_row_block, handle, handle2, &
      i_global, i_local, i_sub, iiB, iii, itmp(2), j_global, j_local, j_sub, &
      jjB, my_num_col_blocks, my_num_row_blocks, mypcol, myprow, ncol_block, &
      ncol_local, npcol, nprow, nrow_block, nrow_local, num_rec_cols, &
      number_of_rec, number_of_send, proc_receive, proc_send, proc_shift, &
      rec_col_end, rec_col_size, rec_col_start, rec_counter, rec_pcol, &
      rec_prow, rec_row_end, rec_row_size, rec_row_start, ref_send_pcol, &
      ref_send_prow, send_counter, send_pcol, send_prow, size_rec_buffer, &
      size_send_buffer, start_col_block, start_row_block
    INTEGER :: stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii_vet, index_col_rec, &
                                                map_rec_size, map_send_size, &
                                                req_send
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: blocks_ranges_col, &
      blocks_ranges_row, grid_2_mepos, grid_ref_2_send_pos, mepos_2_grid
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: part_col, part_row
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! create the full matrix, (num_rows x num_cols)
    NULLIFY(fm_mat)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=num_rows,&
                             ncol_global=num_cols,para_env=para_env)
    ! for now we don't force the blocks
    ! CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=num_rows,&
    !                          ncol_global=num_cols,para_env=para_env,&
    !                          nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)
    CALL cp_fm_create(fm_mat,fm_struct,name="fm_mat")
    CALL cp_fm_set_all(matrix=fm_mat,alpha=0.0_dp)
    CALL cp_fm_struct_release(fm_struct)

    ! start filling procedure
    ! fill the matrix
    CALL cp_fm_get_info(matrix=fm_mat,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block)
    myprow=fm_mat%matrix_struct%context%mepos(1)
    mypcol=fm_mat%matrix_struct%context%mepos(2)
    nprow =fm_mat%matrix_struct%context%num_pe(1)
    npcol =fm_mat%matrix_struct%context%num_pe(2)

    ! start the communication part
    ! 0) create array containing the processes position
    !    and supporting infos
    CALL timeset(routineN//"_info",handle2)
    ALLOCATE(grid_2_mepos(0:nprow-1,0:npcol-1),STAT=stat)
    CPASSERT(stat==0)
    grid_2_mepos=0
    ALLOCATE(mepos_2_grid(0:para_env%num_pe-1,2),STAT=stat)
    CPASSERT(stat==0)
    mepos_2_grid=0
    ! fill the info array
    grid_2_mepos(myprow,mypcol)=para_env%mepos
    mepos_2_grid(para_env%mepos,1)=myprow
    mepos_2_grid(para_env%mepos,2)=mypcol
    ! sum infos
    CALL mp_sum(grid_2_mepos,para_env%group)
    CALL mp_sum(mepos_2_grid,para_env%group)   

    ! 1) loop over my local data and define a map for the proc to send data
    ALLOCATE(map_send_size(0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    map_send_size=0
    DO jjB=my_start_col, my_end_col
      send_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                        fm_mat%matrix_struct%first_p_pos(2),npcol)
      DO iiB=my_start_row, my_end_row
        send_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                          fm_mat%matrix_struct%first_p_pos(1),nprow)
        proc_send=grid_2_mepos(send_prow,send_pcol)
        map_send_size(proc_send)=map_send_size(proc_send)+1
      END DO
    END DO

    ! 2) loop over my local data of fm_mat and define a map for the proc from which rec data
    ALLOCATE(map_rec_size(0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    map_rec_size=0
    part_row=FLOAT(num_rows)/FLOAT(ngroup_row)
    part_col=FLOAT(num_cols)/FLOAT(ngroup_col)


    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      ! check the group holding this element
      ! dirty way, if someone has a better idea ...
      rec_pcol=INT(FLOAT(j_global-1)/part_col)
      rec_pcol=MAX(0,rec_pcol)
      rec_pcol=MIN(rec_pcol,ngroup_col-1)
      DO
        itmp=get_limit(num_cols,ngroup_col,rec_pcol)
        IF(j_global>=itmp(1).AND.j_global<=itmp(2)) EXIT
        IF(j_global<itmp(1)) rec_pcol=rec_pcol-1
        IF(j_global>itmp(2)) rec_pcol=rec_pcol+1
      END DO

      DO iiB=1, nrow_local
        i_global=row_indices(iiB)
        ! check the process in the group holding this element
        rec_prow=INT(FLOAT(i_global-1)/part_row)
        rec_prow=MAX(0,rec_prow)
        rec_prow=MIN(rec_prow,ngroup_row-1)
        DO
          itmp=get_limit(num_rows,ngroup_row,rec_prow)
          IF(i_global>=itmp(1).AND.i_global<=itmp(2)) EXIT
          IF(i_global<itmp(1)) rec_prow=rec_prow-1
          IF(i_global>itmp(2)) rec_prow=rec_prow+1
        END DO

        proc_receive=group_grid_2_mepos(rec_prow,rec_pcol)

        map_rec_size(proc_receive)=map_rec_size(proc_receive)+1

      END DO  ! i_global
    END DO  ! j_global

    ! 3) check if the local data has to be stored in the new fm_mat 
    IF(map_rec_size(para_env%mepos)>0) THEN
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        IF(j_global>=my_start_col.AND.j_global<=my_end_col) THEN
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)
            IF(i_global>=my_start_row.AND.i_global<=my_end_row) THEN
              fm_mat%local_data(iiB,jjB)=mat2D(i_global-my_start_row+1,j_global-my_start_col+1)
            END IF
          END DO
        END IF
      END DO
    END IF
    CALL timestop(handle2)
 
    ! 4) reorder data in the send_buffer
    CALL timeset(routineN//"_buffer_send",handle2)
    ! count the number of messages to send
    number_of_send=0
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        number_of_send=number_of_send+1
      END IF
    END DO
    ! allocate the structure that will hold the messages to be sent
    ALLOCATE(buffer_send(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    ! grid_ref_2_send_pos is an array (map) that given a pair 
    ! (ref_send_prow,ref_send_pcol) returns
    ! the position in the buffer_send associated to that process
    ALLOCATE(grid_ref_2_send_pos(0:nprow-1,0:npcol-1),STAT=stat)
    CPASSERT(stat==0)
    grid_ref_2_send_pos=0
    ! finalize the allocation of buffer_send with the actual size
    ! of each message (actual size is size_send_buffer)
    send_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      size_send_buffer=map_send_size(proc_send)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        ! allocate the sending buffer (msg)
        ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
        CPASSERT(stat==0)
        buffer_send(send_counter)%msg=0.0_dp
        buffer_send(send_counter)%proc=proc_send
        ! get the pointer to prow, pcol of the process that has 
        ! to receive this message
        ref_send_prow=mepos_2_grid(proc_send,1)
        ref_send_pcol=mepos_2_grid(proc_send,2)
        ! save the rank of the process that has to receive this message
        grid_ref_2_send_pos(ref_send_prow,ref_send_pcol)=send_counter
      END IF
    END DO   
 
    ! loop over the locally held data and fill the buffer_send
    ! for doing that we need an array that keep track if the 
    ! sequential increase of the index for each message
    ALLOCATE(iii_vet(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    iii_vet=0
    DO iiB=my_start_row, my_end_row
      send_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                              fm_mat%matrix_struct%first_p_pos(1),nprow)
      DO jjB=my_start_col, my_end_col
        send_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                                fm_mat%matrix_struct%first_p_pos(2),npcol)
        ! we don't need to send to ourselves
        IF(grid_2_mepos(send_prow,send_pcol)==para_env%mepos) CYCLE

        send_counter=grid_ref_2_send_pos(send_prow,send_pcol)
        iii_vet(send_counter)=iii_vet(send_counter)+1
        iii=iii_vet(send_counter)
        buffer_send(send_counter)%msg(iii)=mat2D(iiB-my_start_row+1,jjB-my_start_col+1)
      END DO
    END DO

    DEALLOCATE(iii_vet)
    DEALLOCATE(grid_ref_2_send_pos)
    DEALLOCATE(mat2D)
    CALL timestop(handle2)

    ! 5) similarly to what done for the buffer_send 
    !    create the buffer for receive, post the message with irecv
    !    and send the messages with mp_isend
    CALL timeset(routineN//"_isendrecv",handle2)
    ! count the number of messages to be received
    number_of_rec=0
    DO proc_shift=1, para_env%num_pe-1
      proc_receive=proc_map(para_env%mepos-proc_shift)
      IF(map_rec_size(proc_receive)>0) THEN
        number_of_rec=number_of_rec+1
      END IF
    END DO
    
    ALLOCATE(buffer_rec(number_of_rec),STAT=stat)
    CPASSERT(stat==0)

    rec_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_receive=proc_map(para_env%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)
      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1
        ! prepare the buffer for receive
        ALLOCATE(buffer_rec(rec_counter)%msg(size_rec_buffer),STAT=stat)
        CPASSERT(stat==0)
        buffer_rec(rec_counter)%msg=0.0_dp
        buffer_rec(rec_counter)%proc=proc_receive
        ! post the message to be received
        CALL mp_irecv(buffer_rec(rec_counter)%msg,proc_receive,para_env%group,&
                      buffer_rec(rec_counter)%msg_req)
      END IF
    END DO

    ! send messages
    ALLOCATE(req_send(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    send_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        CALL mp_isend(buffer_send(send_counter)%msg,proc_send,para_env%group,&
                      buffer_send(send_counter)%msg_req)
        req_send(send_counter)=buffer_send(send_counter)%msg_req
      END IF
    END DO
    CALL timestop(handle2)

    ! 6) fill the fm_mat matrix with the messages received from the 
    !    other processes
    CALL timeset(routineN//"_fill",handle2)
    ! In order to perform this step efficiently first we have to know the  
    ! ranges of the blocks over which a given process hold its local data.
    ! Start with the rows ...
    my_num_row_blocks=1
    DO iiB=1, nrow_local-1
      IF(ABS(row_indices(iiB+1)-row_indices(iiB))==1) CYCLE
      my_num_row_blocks=my_num_row_blocks+1
    END DO
    ALLOCATE(blocks_ranges_row(2,my_num_row_blocks),STAT=stat)
    CPASSERT(stat==0)
    blocks_ranges_row=0
    blocks_ranges_row(1,1)=row_indices(1)
    iii=1
    DO iiB=1, nrow_local-1
      IF(ABS(row_indices(iiB+1)-row_indices(iiB))==1) CYCLE
      iii=iii+1
      blocks_ranges_row(2,iii-1)=row_indices(iiB)
      blocks_ranges_row(1,iii)=row_indices(iiB+1)
    END DO
    blocks_ranges_row(2,my_num_row_blocks)=row_indices(MAX(nrow_local,1))
    ! ... and columns
    my_num_col_blocks=1
    DO jjB=1, ncol_local-1
      IF(ABS(col_indices(jjB+1)-col_indices(jjB))==1) CYCLE
      my_num_col_blocks=my_num_col_blocks+1
    END DO
    ALLOCATE(blocks_ranges_col(2,my_num_col_blocks),STAT=stat)
    CPASSERT(stat==0)
    blocks_ranges_col=0
    blocks_ranges_col(1,1)=col_indices(1)
    iii=1
    DO jjB=1, ncol_local-1
      IF(ABS(col_indices(jjB+1)-col_indices(jjB))==1) CYCLE
      iii=iii+1
      blocks_ranges_col(2,iii-1)=col_indices(jjB)
      blocks_ranges_col(1,iii)=col_indices(jjB+1)
    END DO
    blocks_ranges_col(2,my_num_col_blocks)=col_indices(MAX(ncol_local,1))

    rec_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_receive=proc_map(para_env%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)

      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1

        rec_col_size=sizes(2,proc_receive)
        rec_col_start=starts(2,proc_receive)
        rec_col_end=ends(2,proc_receive)

        ! precompute the number of received columns and relative index
        num_rec_cols=0
        DO jjB=1, my_num_col_blocks
          start_col_block=MAX(blocks_ranges_col(1,jjB),rec_col_start)
          end_col_block=MIN(blocks_ranges_col(2,jjB),rec_col_end)
          DO j_sub=start_col_block, end_col_block
            num_rec_cols=num_rec_cols+1
          END DO
        END DO
        ALLOCATE(index_col_rec(num_rec_cols),STAT=stat)
        CPASSERT(stat==0)
        index_col_rec=0
        iii=0
        DO jjB=1, my_num_col_blocks
          start_col_block=MAX(blocks_ranges_col(1,jjB),rec_col_start)
          end_col_block=MIN(blocks_ranges_col(2,jjB),rec_col_end)
          DO j_sub=start_col_block, end_col_block
            iii=iii+1
            j_local=cp_fm_indxg2l(j_sub,ncol_block,dummy_proc,&
                                  fm_mat%matrix_struct%first_p_pos(2),npcol)
            index_col_rec(iii)=j_local
          END DO
        END DO

        rec_row_size=sizes(1,proc_receive)
        rec_row_start=starts(1,proc_receive)
        rec_row_end=ends(1,proc_receive)

        ! wait for the message 
        CALL mp_wait(buffer_rec(rec_counter)%msg_req)

        ! fill the local data in fm_mat
        iii=0
        DO iiB=1, my_num_row_blocks
          start_row_block=MAX(blocks_ranges_row(1,iiB),rec_row_start)
          end_row_block=MIN(blocks_ranges_row(2,iiB),rec_row_end)
          DO i_sub=start_row_block, end_row_block
            i_local=cp_fm_indxg2l(i_sub,nrow_block,dummy_proc,&
                                  fm_mat%matrix_struct%first_p_pos(1),nprow)
            DO jjB=1, num_rec_cols
              iii=iii+1
              j_local=index_col_rec(jjB)
              fm_mat%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
            END DO
          END DO
        END DO

        DEALLOCATE(buffer_rec(rec_counter)%msg)
        DEALLOCATE(index_col_rec)
      END IF
    END DO
    DEALLOCATE(buffer_rec)

    DEALLOCATE(blocks_ranges_row)
    DEALLOCATE(blocks_ranges_col)

    CALL timestop(handle2)

    ! 7) Finally wait for all messeges to be sent
    CALL timeset(routineN//"_waitall",handle2)
    CALL mp_waitall(req_send(:))
    DO send_counter=1, number_of_send
      DEALLOCATE(buffer_send(send_counter)%msg)
    END DO
    DEALLOCATE(buffer_send)
    CALL timestop(handle2)

    DEALLOCATE(map_send_size)
    DEALLOCATE(map_rec_size)
    DEALLOCATE(grid_2_mepos)
    DEALLOCATE(mepos_2_grid)

    CALL timestop(handle)

  END SUBROUTINE array2fm

! *****************************************************************************
!> \brief ...
!> \param mat2D ...
!> \param num_rows ...
!> \param num_cols ...
!> \param para_env ...
!> \param proc_map ...
!> \param my_rows ...
!> \param my_start_row ...
!> \param my_end_row ...
!> \param my_cols ...
!> \param my_start_col ...
!> \param my_end_col ...
!> \param sizes ...
!> \param group_grid_2_mepos ...
!> \param mepos_2_grid_group ...
!> \param ngroup_row ...
!> \param ngroup_col ...
!> \param fm_mat ...
! *****************************************************************************
  SUBROUTINE fm2array(mat2D,num_rows,num_cols,para_env,proc_map,&
                      my_rows,my_start_row,my_end_row,&
                      my_cols,my_start_col,my_end_col,&
                      sizes,group_grid_2_mepos,mepos_2_grid_group,&
                      ngroup_row,ngroup_col,&
                      fm_mat)

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mat2D
    INTEGER                                  :: num_rows, num_cols
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_map
    INTEGER                                  :: my_rows, my_start_row, &
                                                my_end_row, my_cols, &
                                                my_start_col, my_end_col
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: sizes, group_grid_2_mepos, &
                                                mepos_2_grid_group
    INTEGER                                  :: ngroup_row, ngroup_col
    TYPE(cp_fm_type), POINTER                :: fm_mat

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm2array', &
      routineP = moduleN//':'//routineN

    INTEGER :: dummy_proc, handle, handle2, i_global, iiB, iii, itmp(2), &
      j_global, jjB, mypcol, myprow, ncol_block, ncol_local, npcol, nprow, &
      nrow_block, nrow_local, num_rec_rows, number_of_rec, number_of_send, &
      proc_receive, proc_send, proc_shift, rec_col_size, rec_counter, &
      rec_pcol, rec_prow, rec_row_size, ref_send_pcol, ref_send_prow, &
      send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer, &
      stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii_vet, index_row_rec, &
                                                map_rec_size, map_send_size, &
                                                req_send
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: grid_2_mepos, &
                                                grid_ref_2_send_pos, &
                                                mepos_2_grid
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: part_col, part_row
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! allocate the array
    ALLOCATE(mat2D(my_rows,my_cols),STAT=stat)
    CPASSERT(stat==0)
    mat2D=0.0_dp

    ! start procedure
    ! get information of from the full matrix
    CALL cp_fm_get_info(matrix=fm_mat,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block)
    myprow=fm_mat%matrix_struct%context%mepos(1)
    mypcol=fm_mat%matrix_struct%context%mepos(2)
    nprow =fm_mat%matrix_struct%context%num_pe(1)
    npcol =fm_mat%matrix_struct%context%num_pe(2)
    
    ! start the communication part
    ! 0) create array containing the processes position
    !    and supporting infos
    CALL timeset(routineN//"_info",handle2)
    ALLOCATE(grid_2_mepos(0:nprow-1,0:npcol-1),STAT=stat)
    CPASSERT(stat==0)
    grid_2_mepos=0
    ALLOCATE(mepos_2_grid(0:para_env%num_pe-1,2),STAT=stat)
    CPASSERT(stat==0)
    mepos_2_grid=0
    ! ALLOCATE(mepos_2_grid_group(0:para_env%num_pe-1,2),STAT=stat)
    ! CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    ! mepos_2_grid_group=0
    ! fill the info array
    grid_2_mepos(myprow,mypcol)=para_env%mepos
    mepos_2_grid(para_env%mepos,1)=myprow
    mepos_2_grid(para_env%mepos,2)=mypcol
 
    ! sum infos
    CALL mp_sum(grid_2_mepos,para_env%group)
    CALL mp_sum(mepos_2_grid,para_env%group)

    ! 1) loop over my local data and define a map for the proc to send data
    ALLOCATE(map_send_size(0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    map_send_size=0
    part_row=FLOAT(num_rows)/FLOAT(ngroup_row)
    part_col=FLOAT(num_cols)/FLOAT(ngroup_col)   

    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      ! check the group holding this element
      ! dirty way, if someone has a better idea ...
      send_pcol=INT(FLOAT(j_global-1)/part_col)
      send_pcol=MAX(0,send_pcol)
      send_pcol=MIN(send_pcol,ngroup_col-1)
      DO
        itmp=get_limit(num_cols,ngroup_col,send_pcol)
        IF(j_global>=itmp(1).AND.j_global<=itmp(2)) EXIT
        IF(j_global<itmp(1)) send_pcol=send_pcol-1
        IF(j_global>itmp(2)) send_pcol=send_pcol+1
      END DO

      DO iiB=1, nrow_local
        i_global=row_indices(iiB)
        ! check the process in the group holding this element
        send_prow=INT(FLOAT(i_global-1)/part_row)
        send_prow=MAX(0,send_prow)
        send_prow=MIN(send_prow,ngroup_row-1)
        DO
          itmp=get_limit(num_rows,ngroup_row,send_prow)
          IF(i_global>=itmp(1).AND.i_global<=itmp(2)) EXIT
          IF(i_global<itmp(1)) send_prow=send_prow-1
          IF(i_global>itmp(2)) send_prow=send_prow+1
        END DO

        proc_send=group_grid_2_mepos(send_prow,send_pcol)
   
        map_send_size(proc_send)=map_send_size(proc_send)+1

      END DO  ! i_global
    END DO ! j_global

    ! WRITE(*,*) para_env%mepos,'AAAA',map_send_size

    ! 2) loop over my local data of the array and define a map for the proc from which rec data
    ALLOCATE(map_rec_size(0:para_env%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    map_rec_size=0
    DO jjB=my_start_col, my_end_col
      rec_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                             fm_mat%matrix_struct%first_p_pos(2),npcol)
      DO iiB=my_start_row, my_end_row
        rec_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                               fm_mat%matrix_struct%first_p_pos(1),nprow)
        proc_receive=grid_2_mepos(rec_prow,rec_pcol)
        map_rec_size(proc_receive)=map_rec_size(proc_receive)+1
      END DO
    END DO
    
    ! 3) check if the local data in the fm_mat has to be stored in the new array
    IF(map_rec_size(para_env%mepos)>0) THEN
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        IF(j_global>=my_start_col.AND.j_global<=my_end_col) THEN
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)
            IF(i_global>=my_start_row.AND.i_global<=my_end_row) THEN
              mat2D(i_global-my_start_row+1,j_global-my_start_col+1)=fm_mat%local_data(iiB,jjB)
            END IF
          END DO
        END IF
      END DO
    END IF
    CALL timestop(handle2)
    
    ! 4) reorder data in the send_buffer
    CALL timeset(routineN//"_buffer_send",handle2)
    ! count the number of messages to send
    number_of_send=0
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        number_of_send=number_of_send+1
      END IF
    END DO
    ! allocate the structure that will hold the messages to be sent
    ALLOCATE(buffer_send(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    ! grid_ref_2_send_pos is an array (map) that given a pair 
    ! (ref_send_prow,ref_send_pcol) returns
    ! the position in the buffer_send associated to that process

    ! ALLOCATE(grid_ref_2_send_pos(0:nprow-1,0:npcol-1),STAT=stat)
    ALLOCATE(grid_ref_2_send_pos(0:ngroup_row-1,0:ngroup_col-1),STAT=stat)
    CPASSERT(stat==0)
    grid_ref_2_send_pos=0

    ! finalize the allocation of buffer_send with the actual size
    ! of each message (actual size is size_send_buffer)
    send_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      size_send_buffer=map_send_size(proc_send)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        ! allocate the sending buffer (msg)
        ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
        CPASSERT(stat==0)
        buffer_send(send_counter)%msg=0.0_dp
        buffer_send(send_counter)%proc=proc_send
        ! get the pointer to prow, pcol of the process that has 
        ! to receive this message
        ref_send_prow=mepos_2_grid_group(proc_send,1)
        ref_send_pcol=mepos_2_grid_group(proc_send,2)
        ! save the rank of the process that has to receive this message
        grid_ref_2_send_pos(ref_send_prow,ref_send_pcol)=send_counter
      END IF
    END DO

    ! loop over the locally held data and fill the buffer_send
    ! for doing that we need an array that keep track if the 
    ! sequential increase of the index for each message
    ALLOCATE(iii_vet(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    iii_vet=0
    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      ! check the group holding this element
      ! dirty way, if someone has a better idea ...
      send_pcol=INT(FLOAT(j_global-1)/part_col)
      send_pcol=MAX(0,send_pcol)
      send_pcol=MIN(send_pcol,ngroup_col-1)
      DO
        itmp=get_limit(num_cols,ngroup_col,send_pcol)
        IF(j_global>=itmp(1).AND.j_global<=itmp(2)) EXIT
        IF(j_global<itmp(1)) send_pcol=send_pcol-1
        IF(j_global>itmp(2)) send_pcol=send_pcol+1
      END DO

      DO iiB=1, nrow_local
        i_global=row_indices(iiB)
        ! check the process in the group holding this element
        send_prow=INT(FLOAT(i_global-1)/part_row)
        send_prow=MAX(0,send_prow)
        send_prow=MIN(send_prow,ngroup_row-1)
        DO
          itmp=get_limit(num_rows,ngroup_row,send_prow)
          IF(i_global>=itmp(1).AND.i_global<=itmp(2)) EXIT
          IF(i_global<itmp(1)) send_prow=send_prow-1
          IF(i_global>itmp(2)) send_prow=send_prow+1
        END DO
        ! we don't need to send to ourselves
        IF(group_grid_2_mepos(send_prow,send_pcol)==para_env%mepos) CYCLE

        send_counter=grid_ref_2_send_pos(send_prow,send_pcol)
        iii_vet(send_counter)=iii_vet(send_counter)+1
        iii=iii_vet(send_counter)
        buffer_send(send_counter)%msg(iii)=fm_mat%local_data(iiB,jjB)
      END DO
    END DO

    DEALLOCATE(iii_vet)
    DEALLOCATE(grid_ref_2_send_pos)
    ! CALL cp_fm_release(fm_mat)
    CALL timestop(handle2)

    ! 5) similarly to what done for the buffer_send 
    !    create the buffer for receive, post the message with irecv
    !    and send the messages with mp_isend
    CALL timeset(routineN//"_isendrecv",handle2)
    ! count the number of messages to be received
    number_of_rec=0
    DO proc_shift=1, para_env%num_pe-1
      proc_receive=proc_map(para_env%mepos-proc_shift)
      IF(map_rec_size(proc_receive)>0) THEN
        number_of_rec=number_of_rec+1
      END IF
    END DO

    ALLOCATE(buffer_rec(number_of_rec),STAT=stat)
    CPASSERT(stat==0)
    
    rec_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_receive=proc_map(para_env%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)
      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1
        ! prepare the buffer for receive
        ALLOCATE(buffer_rec(rec_counter)%msg(size_rec_buffer),STAT=stat)
        CPASSERT(stat==0)
        buffer_rec(rec_counter)%msg=0.0_dp
        buffer_rec(rec_counter)%proc=proc_receive
        ! post the message to be received
        CALL mp_irecv(buffer_rec(rec_counter)%msg,proc_receive,para_env%group,&
                      buffer_rec(rec_counter)%msg_req)
      END IF
    END DO    

    ! send messages
    ALLOCATE(req_send(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    send_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        CALL mp_isend(buffer_send(send_counter)%msg,proc_send,para_env%group,&
                      buffer_send(send_counter)%msg_req)
        req_send(send_counter)=buffer_send(send_counter)%msg_req
      END IF
    END DO
    CALL timestop(handle2)

    ! 6) fill the fm_mat matrix with the messages received from the 
    !    other processes
    CALL timeset(routineN//"_fill",handle2)
    iiB=MAXVAL(sizes(1,:))
    ALLOCATE(index_row_rec(iiB),STAT=stat)
    CPASSERT(stat==0)
    index_row_rec=0
    rec_counter=0
    DO proc_shift=1, para_env%num_pe-1
      proc_receive=proc_map(para_env%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)

      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1

        rec_col_size=sizes(2,proc_receive)
        rec_row_size=sizes(1,proc_receive)
 
        ! precompute the indeces of the rows
        num_rec_rows=0
        DO iiB=1, rec_row_size
          i_global=cp_fm_indxl2g(iiB,nrow_block,mepos_2_grid(proc_receive,1),&
                                 fm_mat%matrix_struct%first_p_pos(1),nprow)
          IF(i_global>=my_start_row.AND.i_global<=my_end_row) THEN
            num_rec_rows=num_rec_rows+1
            index_row_rec(num_rec_rows)=i_global
          END IF
        END DO

        CALL mp_wait(buffer_rec(rec_counter)%msg_req)

        iii=0
        DO jjB=1, rec_col_size
          j_global=cp_fm_indxl2g(jjB,ncol_block,mepos_2_grid(proc_receive,2),&
                                 fm_mat%matrix_struct%first_p_pos(2),npcol)
          IF(j_global>=my_start_col.AND.j_global<=my_end_col) THEN 
            DO iiB=1, num_rec_rows
              i_global=index_row_rec(iiB)
              iii=iii+1
              mat2D(i_global-my_start_row+1,j_global-my_start_col+1)=buffer_rec(rec_counter)%msg(iii)
            END DO
          END IF
        END DO

        DEALLOCATE(buffer_rec(rec_counter)%msg)
      END IF
    END DO
    DEALLOCATE(buffer_rec)
    DEALLOCATE(index_row_rec)
    CALL cp_fm_release(fm_mat)
    CALL timestop(handle2)

    ! 7) Finally wait for all messeges to be sent
    CALL timeset(routineN//"_waitall",handle2)
    CALL mp_waitall(req_send(:))
    DO send_counter=1, number_of_send
      DEALLOCATE(buffer_send(send_counter)%msg)
    END DO
    DEALLOCATE(buffer_send)
    CALL timestop(handle2)

    CALL timestop(handle)

  END SUBROUTINE fm2array

! *****************************************************************************
!> \brief ...
!> \param Gamma_2D ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_ia ...
!> \param para_env_sub ...
!> \param sub_proc_map ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \param sizes_ia ...
!> \param starts_ia ...
!> \param ends_ia ...
!> \param mp2_env ...
!> \param alpha_case ...
! *****************************************************************************
  SUBROUTINE create_dbcsr_gamma(Gamma_2D,homo,virtual,dimen_ia,para_env_sub,sub_proc_map,&
                                my_ia_start,my_ia_end,my_group_L_size,&
                                sizes_ia,starts_ia,ends_ia,&
                                mp2_env,alpha_case)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Gamma_2D
    INTEGER                                  :: homo, virtual, dimen_ia
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sub_proc_map
    INTEGER                                  :: my_ia_start, my_ia_end, &
                                                my_group_L_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_ia, starts_ia, ends_ia
    TYPE(mp2_type), POINTER                  :: mp2_env
    LOGICAL                                  :: alpha_case

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dbcsr_gamma', &
      routineP = moduleN//':'//routineN

    INTEGER :: dummy_proc, handle, i_global, i_local, iaia, iiB, iii, &
      itmp(2), j_global, j_local, jjB, jjj, kkB, mypcol, myprow, ncol_block, &
      ncol_local, npcol, nprow, nrow_block, nrow_local, number_of_rec, &
      number_of_send, proc_receive, proc_send, proc_shift, rec_counter, &
      rec_iaia_end, rec_iaia_size, rec_iaia_start, rec_pcol, rec_prow, &
      ref_send_pcol, ref_send_prow, send_counter, send_pcol, send_prow, &
      size_rec_buffer, size_send_buffer, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii_vet, map_rec_size, &
                                                map_send_size, req_send
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: grid_2_mepos, &
                                                grid_ref_2_send_pos, &
                                                indeces_map_my, mepos_2_grid
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: part_ia
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_ia
    TYPE(index_map), ALLOCATABLE, &
      DIMENSION(:)                           :: indeces_rec
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send

! TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: G_P_ia

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! allocate G_P_ia (DBCSR)
    IF (alpha_case) THEN
       NULLIFY(mp2_env%ri_grad%G_P_ia)
       CALL cp_dbcsr_allocate_matrix_set(mp2_env%ri_grad%G_P_ia, &
            my_group_L_size)
    ELSE
       NULLIFY(mp2_env%ri_grad%G_P_ia_beta)
       CALL cp_dbcsr_allocate_matrix_set(mp2_env%ri_grad%G_P_ia_beta, &
            my_group_L_size)
    ENDIF

    ! create sub blacs env
    NULLIFY(blacs_env)
    CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_sub)

    ! create the fm_ia buffer matrix
    NULLIFY(fm_ia)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=homo,&
                             ncol_global=virtual,para_env=para_env_sub)
    CALL cp_fm_create(fm_ia,fm_struct,name="fm_ia")
    ! release structure 
    CALL cp_fm_struct_release(fm_struct)   
    ! release blacs_env
    CALL cp_blacs_env_release(blacs_env)

    ! get array information
    CALL cp_fm_get_info(matrix=fm_ia,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block)
    myprow=fm_ia%matrix_struct%context%mepos(1)
    mypcol=fm_ia%matrix_struct%context%mepos(2)
    nprow =fm_ia%matrix_struct%context%num_pe(1)
    npcol =fm_ia%matrix_struct%context%num_pe(2)

    ! 0) create array containing the processes position and supporting infos
    ALLOCATE(grid_2_mepos(0:nprow-1,0:npcol-1),STAT=stat)
    CPASSERT(stat==0)
    grid_2_mepos=0
    ALLOCATE(mepos_2_grid(0:para_env_sub%num_pe-1,2),STAT=stat)
    CPASSERT(stat==0)
    mepos_2_grid=0
    ! fill the info array
    grid_2_mepos(myprow,mypcol)=para_env_sub%mepos
    mepos_2_grid(para_env_sub%mepos,1)=myprow
    mepos_2_grid(para_env_sub%mepos,2)=mypcol
    ! sum infos
    CALL mp_sum(grid_2_mepos,para_env_sub%group)
    CALL mp_sum(mepos_2_grid,para_env_sub%group)

    ! loop over local index range and define the sending map
    ALLOCATE(map_send_size(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    map_send_size=0   
    DO iaia=my_ia_start, my_ia_end
      i_global=(iaia-1)/virtual+1
      j_global=MOD(iaia-1,virtual)+1
      send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
                              fm_ia%matrix_struct%first_p_pos(1),nprow)
      send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
                              fm_ia%matrix_struct%first_p_pos(2),npcol)
      proc_send=grid_2_mepos(send_prow,send_pcol)
      map_send_size(proc_send)=map_send_size(proc_send)+1
    END DO

    ! loop over local data of fm_ia and define the receiving map
    ALLOCATE(map_rec_size(0:para_env_sub%num_pe-1),STAT=stat)
    CPASSERT(stat==0)
    map_rec_size=0
    part_ia=FLOAT(dimen_ia)/FLOAT(para_env_sub%num_pe)

    DO iiB=1, nrow_local
      i_global=row_indices(iiB)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        iaia=(i_global-1)*virtual+j_global
        proc_receive=INT(FLOAT(iaia-1)/part_ia)
        proc_receive=MAX(0,proc_receive)
        proc_receive=MIN(proc_receive,para_env_sub%num_pe-1)
        DO 
          itmp=get_limit(dimen_ia,para_env_sub%num_pe,proc_receive)
          IF(iaia>=itmp(1).AND.iaia<=itmp(2)) EXIT
          IF(iaia<itmp(1)) proc_receive=proc_receive-1
          IF(iaia>itmp(2)) proc_receive=proc_receive+1
        END DO
        map_rec_size(proc_receive)=map_rec_size(proc_receive)+1
      END DO
    END DO

    ! allocate the buffer for sending data
    number_of_send=0
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      IF(map_send_size(proc_send)>0) THEN
        number_of_send=number_of_send+1
      END IF
    END DO
    ! allocate the structure that will hold the messages to be sent
    ALLOCATE(buffer_send(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    ! and the map from the grid of processess to the message position
    ALLOCATE(grid_ref_2_send_pos(0:nprow-1,0:npcol-1),STAT=stat)
    CPASSERT(stat==0)
    grid_ref_2_send_pos=0
    ! finally allocate each message
    send_counter=0
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      size_send_buffer=map_send_size(proc_send)
      IF(map_send_size(proc_send)>0) THEN
        send_counter=send_counter+1
        ! allocate the sending buffer (msg)
        ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
        CPASSERT(stat==0)
        buffer_send(send_counter)%proc=proc_send
        ! get the pointer to prow, pcol of the process that has 
        ! to receive this message
        ref_send_prow=mepos_2_grid(proc_send,1)
        ref_send_pcol=mepos_2_grid(proc_send,2)
        ! save the rank of the process that has to receive this message
        grid_ref_2_send_pos(ref_send_prow,ref_send_pcol)=send_counter
      END IF
    END DO
       

    ! allocate the buffer for receiving data
    number_of_rec=0
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
      IF(map_rec_size(proc_receive)>0) THEN
        number_of_rec=number_of_rec+1
      END IF
    END DO
    ! allocate the structure that will hold the messages to be received
    ! and relative indeces
    ALLOCATE(buffer_rec(number_of_rec),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(indeces_rec(number_of_rec),STAT=stat)
    CPASSERT(stat==0)
    ! finally allocate each message and fill the array of indeces
    rec_counter=0
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
      size_rec_buffer=map_rec_size(proc_receive)
      IF(map_rec_size(proc_receive)>0) THEN
        rec_counter=rec_counter+1
        ! prepare the buffer for receive
        ALLOCATE(buffer_rec(rec_counter)%msg(size_rec_buffer),STAT=stat)
        CPASSERT(stat==0)
        buffer_rec(rec_counter)%proc=proc_receive
        ! create the indeces array
        ALLOCATE(indeces_rec(rec_counter)%map(2,size_rec_buffer),STAT=stat)
        CPASSERT(stat==0)
        indeces_rec(rec_counter)%map=0
        rec_iaia_start=starts_ia(proc_receive)
        rec_iaia_end=ends_ia(proc_receive)
        rec_iaia_size=sizes_ia(proc_receive)
        iii=0
        DO iaia=rec_iaia_start, rec_iaia_end
          i_global=(iaia-1)/virtual+1
          j_global=MOD(iaia-1,virtual)+1
          rec_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
                                 fm_ia%matrix_struct%first_p_pos(1),nprow)
          rec_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
                                 fm_ia%matrix_struct%first_p_pos(2),npcol)
          IF(grid_2_mepos(rec_prow,rec_pcol)/=para_env_sub%mepos) CYCLE
          iii=iii+1
          i_local=cp_fm_indxg2l(i_global,nrow_block,dummy_proc,&
                                fm_ia%matrix_struct%first_p_pos(1),nprow)
          j_local=cp_fm_indxg2l(j_global,ncol_block,dummy_proc,&
                                fm_ia%matrix_struct%first_p_pos(2),npcol)
          indeces_rec(rec_counter)%map(1,iii)=i_local
          indeces_rec(rec_counter)%map(2,iii)=j_local
        END DO
      END IF
    END DO
    ! and create the index map for my local data
    IF(map_rec_size(para_env_sub%mepos)>0) THEN
      size_rec_buffer=map_rec_size(para_env_sub%mepos)
      ALLOCATE(indeces_map_my(2,size_rec_buffer),STAT=stat)
      CPASSERT(stat==0)
      indeces_map_my=0
      iii=0
      DO iaia=my_ia_start, my_ia_end
        i_global=(iaia-1)/virtual+1
        j_global=MOD(iaia-1,virtual)+1
        rec_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
                               fm_ia%matrix_struct%first_p_pos(1),nprow)
        rec_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
                               fm_ia%matrix_struct%first_p_pos(2),npcol)
        IF(grid_2_mepos(rec_prow,rec_pcol)/=para_env_sub%mepos) CYCLE
        iii=iii+1
        i_local=cp_fm_indxg2l(i_global,nrow_block,dummy_proc,&
                              fm_ia%matrix_struct%first_p_pos(1),nprow)
        j_local=cp_fm_indxg2l(j_global,ncol_block,dummy_proc,&
                              fm_ia%matrix_struct%first_p_pos(2),npcol)
        indeces_map_my(1,iii)=i_local
        indeces_map_my(2,iii)=j_local
      END DO     
    END IF

    ! auxiliary vector of indeces for the send buffer
    ALLOCATE(iii_vet(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    ! vector for the send requests
    ALLOCATE(req_send(number_of_send),STAT=stat)
    CPASSERT(stat==0)
    ! loop over auxiliary basis function and redistribute into a fm
    ! and then compy the fm into a dbcsr matrix
    DO kkB=1, my_group_L_size
      ! zero the matries of the buffers and post the messages to be received
      CALL cp_fm_set_all(matrix=fm_ia,alpha=0.0_dp)
      rec_counter=0
      DO proc_shift=1, para_env_sub%num_pe-1
        proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
        IF(map_rec_size(proc_receive)>0) THEN
          rec_counter=rec_counter+1
          buffer_rec(rec_counter)%msg=0.0_dp
          CALL mp_irecv(buffer_rec(rec_counter)%msg,proc_receive,para_env_sub%group,&
                        buffer_rec(rec_counter)%msg_req)
        END IF
      END DO
      ! fill the sending buffer and send the messages
      DO send_counter=1, number_of_send
        buffer_send(send_counter)%msg=0.0_dp
      END DO
      iii_vet=0
      jjj=0
      DO iaia=my_ia_start, my_ia_end
        i_global=(iaia-1)/virtual+1
        j_global=MOD(iaia-1,virtual)+1
        send_prow=cp_fm_indxg2p(i_global,nrow_block,dummy_proc,&
                                fm_ia%matrix_struct%first_p_pos(1),nprow)
        send_pcol=cp_fm_indxg2p(j_global,ncol_block,dummy_proc,&
                                fm_ia%matrix_struct%first_p_pos(2),npcol)
        proc_send=grid_2_mepos(send_prow,send_pcol)
        ! we don't need to send to ourselves
        IF(grid_2_mepos(send_prow,send_pcol)==para_env_sub%mepos) THEN
          ! filling fm_ia with local data
          jjj=jjj+1
          i_local=indeces_map_my(1,jjj)
          j_local=indeces_map_my(2,jjj)
          fm_ia%local_data(i_local,j_local)=Gamma_2D(iaia-my_ia_start+1,kkB)
        ELSE
          send_counter=grid_ref_2_send_pos(send_prow,send_pcol)
          iii_vet(send_counter)=iii_vet(send_counter)+1
          iii=iii_vet(send_counter)
          buffer_send(send_counter)%msg(iii)=Gamma_2D(iaia-my_ia_start+1,kkB)
        END IF
      END DO
      req_send=0  
      send_counter=0
      DO proc_shift=1, para_env_sub%num_pe-1
        proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
        IF(map_send_size(proc_send)>0) THEN
          send_counter=send_counter+1
          CALL mp_isend(buffer_send(send_counter)%msg,proc_send,para_env_sub%group,&
                        buffer_send(send_counter)%msg_req)
          req_send(send_counter)=buffer_send(send_counter)%msg_req
        END IF
      END DO
   
      ! receive the massages and fill the fm_ia
      rec_counter=0
      DO proc_shift=1, para_env_sub%num_pe-1
        proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)
        size_rec_buffer=map_rec_size(proc_receive)
        IF(map_rec_size(proc_receive)>0) THEN
          rec_counter=rec_counter+1
          ! wait for the message
          CALL mp_wait(buffer_rec(rec_counter)%msg_req)
          DO iii=1, size_rec_buffer
            i_local=indeces_rec(rec_counter)%map(1,iii)
            j_local=indeces_rec(rec_counter)%map(2,iii)
            fm_ia%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
          END DO
        END IF
      END DO

      ! wait all
      CALL mp_waitall(req_send(:))

      ! now create the DBCSR matrix and copy fm_ia into it
      IF (alpha_case) THEN
        ALLOCATE(mp2_env%ri_grad%G_P_ia(kkB)%matrix)
        CALL cp_dbcsr_init(mp2_env%ri_grad%G_P_ia(kkB)%matrix)
        CALL cp_dbcsr_m_by_n_from_template(mp2_env%ri_grad%G_P_ia(kkB)%matrix, &
             template=mp2_env%ri_grad%mo_coeff_o,&
             m=homo,n=virtual, sym=dbcsr_type_no_symmetry)
        CALL copy_fm_to_dbcsr(fm_ia,mp2_env%ri_grad%G_P_ia(kkB)%matrix,keep_sparsity=.FALSE.)
      ELSE 
        ALLOCATE(mp2_env%ri_grad%G_P_ia_beta(kkB)%matrix)
        CALL cp_dbcsr_init(mp2_env%ri_grad%G_P_ia_beta(kkB)%matrix)
        CALL cp_dbcsr_m_by_n_from_template(mp2_env%ri_grad%G_P_ia_beta(kkB)%matrix, &
             template=mp2_env%ri_grad%mo_coeff_o_beta,&
             m=homo,n=virtual, sym=dbcsr_type_no_symmetry)
        CALL copy_fm_to_dbcsr(fm_ia,mp2_env%ri_grad%G_P_ia_beta(kkB)%matrix,keep_sparsity=.FALSE.)
      ENDIF

    END DO

    ! deallocate stuff
    DEALLOCATE(Gamma_2D)
    DEALLOCATE(iii_vet)
    DEALLOCATE(req_send)
    IF(map_rec_size(para_env_sub%mepos)>0) THEN
      DEALLOCATE(indeces_map_my)
    END IF
    DO rec_counter=1, number_of_rec
      DEALLOCATE(indeces_rec(rec_counter)%map)
      DEALLOCATE(buffer_rec(rec_counter)%msg)
    END DO
    DEALLOCATE(indeces_rec)
    DEALLOCATE(buffer_rec)
    DO send_counter=1, number_of_send
      DEALLOCATE(buffer_send(send_counter)%msg)
    END DO
    DEALLOCATE(buffer_send)
    DEALLOCATE(map_send_size)
    DEALLOCATE(map_rec_size)
    DEALLOCATE(grid_2_mepos)
    DEALLOCATE(mepos_2_grid)

    ! release buffer matrix
    CALL cp_fm_release(fm_ia)

    CALL timestop(handle)

   END SUBROUTINE create_dbcsr_gamma

! *****************************************************************************
!> \brief ...
!> \param mat ...
! *****************************************************************************
   SUBROUTINE write_array(mat)
    REAL(KIND=dp), DIMENSION(:, :)           :: mat

    INTEGER                                  :: iii, jjj

     WRITE(*,*)
     DO iii=1, SIZE(mat,1)
       WRITE(*,*) iii
       DO jjj=1, SIZE(mat,2), 20
         WRITE(*,'(1000F10.5)') mat(iii,jjj:MIN(SIZE(mat,2),jjj+19))
       END DO
       WRITE(*,*)
     END DO
     WRITE(*,*)
  END SUBROUTINE

END MODULE mp2_ri_grad_util

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines for RPA with imaginary time
!> \par History
!>      10.2015 created [Jan Wilhelm]
! *****************************************************************************
MODULE rpa_im_time
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_copy, cp_dbcsr_filter, cp_dbcsr_get_num_blocks, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_set, cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_comm_split_direct,&
                                             mp_irecv,&
                                             mp_isend,&
                                             mp_max,&
                                             mp_sendrecv,&
                                             mp_sum,&
                                             mp_sync,&
                                             mp_waitall
  USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                             mp2_type
  USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time'

  PUBLIC :: fill_fm_with_2d_array, &
            compute_fm_mat_Q_omega
  CONTAINS

! *****************************************************************************
!> \brief compute the matrix Q(it) (intermediate) and Fourier transform it 
!>        directly to fm_mat_Q_omega(iw) (output)
!> \param fm_mat_Q_omega ...
!> \param fm_scaled_dm_occ ...
!> \param fm_scaled_dm_virt ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mat_Q ...
!> \param mat_M_munu_occ ...
!> \param mat_M_munu_virt ...
!> \param mat_M_munu_P_occ ...
!> \param mat_M_munu_P_virt ...
!> \param mat_B_munu ...
!> \param mat_munu ...
!> \param trace_2d ...
!> \param tau_wj_2d ...
!> \param tj ...
!> \param wj ...
!> \param first_cycle_im_time ...
!> \param do_minimax_quad ...
!> \param tau ...
!> \param tau_old ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param a_scaling ...
!> \param alpha ...
!> \param Eigenval ...
!> \param nmo ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param dimen_RI ...
!> \param max_size_RI_group ...
!> \param ngroup ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
! *****************************************************************************
  SUBROUTINE compute_fm_mat_Q_omega(fm_mat_Q_omega,fm_scaled_dm_occ,&
                                    fm_scaled_dm_virt,fm_mo_coeff_occ,fm_mo_coeff_virt,&
                                    fm_mat_Q,mat_M_munu_occ,mat_M_munu_virt,mat_M_munu_P_occ,&
                                    mat_M_munu_P_virt,mat_B_munu,mat_munu,trace_2d,tau_wj_2d,tj,wj,&
                                    first_cycle_im_time,do_minimax_quad,&
                                    tau,tau_old,e_fermi,eps_filter,a_scaling,alpha,&
                                    Eigenval,nmo,my_group_L_size,&
                                    my_group_L_start,dimen_RI,max_size_RI_group,ngroup,&
                                    num_integ_points, jquad,mp2_env,para_env,para_env_sub)

    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fm_mat_Q_omega
    TYPE(cp_fm_type), POINTER                :: fm_scaled_dm_occ, &
                                                fm_scaled_dm_virt, &
                                                fm_mo_coeff_occ, &
                                                fm_mo_coeff_virt, fm_mat_Q
    TYPE(cp_dbcsr_type), POINTER             :: mat_M_munu_occ, &
                                                mat_M_munu_virt
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_M_munu_P_occ, &
                                                mat_M_munu_P_virt, mat_B_munu
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: trace_2d, tau_wj_2d
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tj, wj
    LOGICAL                                  :: first_cycle_im_time, &
                                                do_minimax_quad
    REAL(KIND=dp)                            :: tau, tau_old, e_fermi, &
                                                eps_filter, a_scaling, alpha
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER :: nmo, my_group_L_size, my_group_L_start, dimen_RI, &
      max_size_RI_group, ngroup, num_integ_points, jquad
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_mat_Q_omega', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle4, handle5, &
                                                i_global, iiB, iquad, jjB, &
                                                LLL, ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: first_cycle_omega_loop
    REAL(KIND=dp)                            :: omega, omega_old, weight, &
                                                weight_old

    CALL timeset(routineN,handle)

    ! Multiply the occupied and the virtual MO coefficients with the factor exp((-e_i-e_F)*tau/2).
    ! Then, we simply get the sum over all occ states and virt. states by a simple matrix-matrix
    ! multiplication.

    ! get info of fm_mo_coeff_occ
    CALL cp_fm_get_info(matrix=fm_mo_coeff_occ,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices)


    CALL timeset(routineN//"_im_time_calc_scaled_dm",handle4)

    ! first, the occ
    IF(first_cycle_im_time) THEN

      DO jjB=1,nrow_local
        DO iiB=1,ncol_local
          i_global=col_indices(iiB)

          fm_mo_coeff_occ%local_data(jjB,iiB) = &
                       fm_mo_coeff_occ%local_data(jjB,iiB)*EXP(tau*0.5_dp*(Eigenval(i_global)-e_fermi))
        END DO
      END DO

    ELSE

      DO jjB=1,nrow_local
        DO iiB=1,ncol_local
          i_global=col_indices(iiB)

          fm_mo_coeff_occ%local_data(jjB,iiB) = &
                       fm_mo_coeff_occ%local_data(jjB,iiB)*EXP((tau-tau_old)*0.5_dp*(Eigenval(i_global)-e_fermi))
        END DO
      END DO

    END IF

    ! get info of fm_mo_coeff_virt
    CALL cp_fm_get_info(matrix=fm_mo_coeff_virt,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices)

    ! the same for virt
    IF(first_cycle_im_time) THEN

      DO jjB=1,nrow_local
        DO iiB=1,ncol_local
          i_global=col_indices(iiB)

          fm_mo_coeff_virt%local_data(jjB,iiB) = &
                       fm_mo_coeff_virt%local_data(jjB,iiB)*EXP(-tau*0.5_dp*(Eigenval(i_global)-e_fermi))
        END DO
      END DO

    ELSE

      DO jjB=1,nrow_local
        DO iiB=1,ncol_local
          i_global=col_indices(iiB)

          fm_mo_coeff_virt%local_data(jjB,iiB) = &
                       fm_mo_coeff_virt%local_data(jjB,iiB)*EXP(-(tau-tau_old)*0.5_dp*(Eigenval(i_global)-e_fermi))
        END DO
      END DO

    END IF

    CALL cp_gemm(transa="N",transb="T",m=nmo,n=nmo,k=nmo,alpha=1.0_dp,&
                    matrix_a=fm_mo_coeff_occ,matrix_b=fm_mo_coeff_occ,beta=0.0_dp,&
                    matrix_c=fm_scaled_dm_occ)

    CALL cp_gemm(transa="N",transb="T",m=nmo,n=nmo,k=nmo,alpha=1.0_dp,&
                    matrix_a=fm_mo_coeff_virt,matrix_b=fm_mo_coeff_virt,beta=0.0_dp,&
                    matrix_c=fm_scaled_dm_virt)

    CALL timestop(handle4)


    CALL timeset(routineN//"_im_time_repl_subgr",handle4)

    ! replicate fm_scaled_dm_occ to every subgroup into sparse matrix mat_M_munu_occ
    CALL replicate_mat_to_subgroup_im_time(mp2_env,para_env,para_env_sub,fm_scaled_dm_occ,nmo,&
                                           mat_munu,mat_M_munu_occ)

    CALL cp_dbcsr_filter(mat_M_munu_occ,eps_filter)

    ! the same for the virtual scaled density matrix
    CALL replicate_mat_to_subgroup_im_time(mp2_env,para_env,para_env_sub,fm_scaled_dm_virt,nmo,&
                                           mat_munu,mat_M_munu_virt)

    CALL cp_dbcsr_filter(mat_M_munu_virt,eps_filter)


    CALL timestop(handle4)

    CALL timeset(routineN//"_mult_3cERI_sc_DM",handle5)

    DO LLL=1,my_group_L_size

      ! transposed or not does not matter because both matrices being multiplied with each other
      ! are symmetric (other order as with occ, because for the trace we actually need one matrix 
      ! transposed (and transposing is easily done in the previous multiplication)
      CALL cp_dbcsr_multiply("T","N",1.0_dp,mat_M_munu_occ,mat_B_munu(LLL)%matrix, &
                             0.0_dp, mat_M_munu_P_occ(LLL)%matrix, filter_eps=eps_filter)

      CALL cp_dbcsr_filter(mat_M_munu_P_occ(LLL)%matrix,eps_filter)

      ! alpha=2.0_dp only for closed-shell (for imag. frequency it is a factor 4.0)
      CALL cp_dbcsr_multiply("T","N",alpha/2.0_dp,mat_B_munu(LLL)%matrix,mat_M_munu_virt, &
                             0.0_dp, mat_M_munu_P_virt(LLL)%matrix, filter_eps=eps_filter)

      CALL cp_dbcsr_filter(mat_M_munu_P_virt(LLL)%matrix,eps_filter)

    END DO

    ! in case we have different my_group_L_size for different groups, we generate the 
    ! the dbcsr matrix format for all mat_M_munu_P_virt matrices

    ! JW should check whether this is really necessary
    IF(my_group_L_size < max_size_RI_group) THEN

      DO LLL=my_group_L_size+1,max_size_RI_group

        CALL cp_dbcsr_copy(mat_M_munu_P_virt(LLL)%matrix, mat_M_munu_P_virt(LLL-1)%matrix)

        CALL cp_dbcsr_set(mat_M_munu_P_virt(LLL)%matrix,0.0_dp)

        CALL cp_dbcsr_filter(mat_M_munu_P_virt(LLL)%matrix,eps_filter)

      END DO

    END IF

    CALL timestop(handle5)

    trace_2d = 0.0_dp

    CALL dbcsr_3index_trace(trace_2d,mat_M_munu_P_occ,mat_M_munu_P_virt,para_env,para_env_sub,&
                            ngroup,my_group_L_size,my_group_L_start)


    ! fill fm_mat_Q 
    CALL fill_fm_with_2d_array(fm_mat_Q,trace_2d,para_env,my_group_L_start,my_group_L_size,&
                               dimen_RI,ngroup)


    CALL timeset(routineN//"_Fourier_transform",handle5)

    ! Fourier transform of Q(it) to Q(iw) (fm_mat_Q_omega)
    first_cycle_omega_loop = .TRUE.

    DO iquad=1, num_integ_points

      IF(do_minimax_quad) THEN
        omega=tj(iquad)
        weight = tau_wj_2d(iquad,jquad)
      ELSE
        omega=a_scaling/TAN(tj(iquad))
        weight = wj(jquad)
      END IF

      IF(first_cycle_omega_loop) THEN
        ! no multiplication with 2.0 as in Kresses paper (Kaltak, JCTC 10, 2498 (2014), Eq. 12)
        ! because this factor is already absorbed in the weight w_j  
        CALL cp_fm_scale(COS(omega*tau)*weight,fm_mat_Q)
      ELSE
        CALL cp_fm_scale(COS(omega*tau)/COS(omega_old*tau)*weight/weight_old,fm_mat_Q)
      END IF

      CALL cp_fm_scale_and_add(1.0_dp,fm_mat_Q_omega(iquad)%matrix,1.0_dp,fm_mat_Q)

      first_cycle_omega_loop = .FALSE.

      omega_old = omega
      weight_old = weight

    END DO

    CALL timestop(handle5)

    CALL timestop(handle)

  END SUBROUTINE 


! *****************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param fm_scaled_dm ...
!> \param nmo ...
!> \param mat_munu ...
!> \param mat_M_munu ...
! *****************************************************************************
  SUBROUTINE replicate_mat_to_subgroup_im_time(mp2_env,para_env,para_env_sub,fm_scaled_dm,&
                                               nmo,mat_munu,mat_M_munu)
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    TYPE(cp_fm_type), POINTER                :: fm_scaled_dm
    INTEGER                                  :: nmo
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(cp_dbcsr_type), POINTER             :: mat_M_munu

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'replicate_mat_to_subgroup_im_time', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_offset, col_size, handle, i_global, iiB, iproc, &
      itmp(2), j_global, jjB, max_row_col_local, my_mu_end, my_mu_size, &
      my_mu_start, ncol_local, ncol_rec, nrow_local, nrow_rec, proc_receive, &
      proc_receive_static, proc_send, proc_send_static, proc_shift, row, &
      row_offset, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_array, proc_map, &
                                                sizes_array, starts_array
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_col_row_info, &
                                                rec_col_row_info
    INTEGER, DIMENSION(:), POINTER           :: col_indices, col_indices_rec, &
                                                row_indices, row_indices_rec
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C, rec_C
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block, local_C, &
                                                local_C_internal
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_coeff

    CALL timeset(routineN,handle)

    ALLOCATE(sizes_array(0:para_env_sub%num_pe-1))
    ALLOCATE(starts_array(0:para_env_sub%num_pe-1))
    starts_array=0
    ALLOCATE(ends_array(0:para_env_sub%num_pe-1))
    ends_array=0

    DO iproc=0,para_env_sub%num_pe-1
       itmp=get_limit(nmo,para_env_sub%num_pe,iproc) 
       starts_array(iproc)=itmp(1)
       ends_array(iproc)=itmp(2)
       sizes_array(iproc)=itmp(2)-itmp(1)+1
    ENDDO
 
    my_mu_size=sizes_array(para_env_sub%mepos)
    my_mu_start=starts_array(para_env_sub%mepos)
    my_mu_end=ends_array(para_env_sub%mepos)

    ! local storage for the C matrix
    ALLOCATE(C(my_mu_size,nmo))
    C=0.0_dp

    ! proc_map, vector that replicate the processor numbers also
    ! for negative and positive number > num_pe
    ! needed to know which is the processor, to respect to another one,
    ! for a given shift
    ALLOCATE(proc_map(-para_env%num_pe:2*para_env%num_pe-1))
    DO iiB=0,para_env%num_pe-1
      proc_map(iiB)=iiB
      proc_map(-iiB-1)=para_env%num_pe-iiB-1
      proc_map(para_env%num_pe+iiB)=iiB
    END DO

    CALL cp_fm_get_info(matrix=fm_scaled_dm,&
                        matrix_struct=fm_struct_coeff,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        local_data=local_C_internal)

    ALLOCATE(local_C(nrow_local,ncol_local))
    local_C=local_C_internal(1:nrow_local,1:ncol_local)
    NULLIFY(local_C_internal) 

    max_row_col_local=MAX(nrow_local,ncol_local)
    CALL mp_max(max_row_col_local,para_env%group)

    ALLOCATE(local_col_row_info(0:max_row_col_local,2))
    local_col_row_info=0
    ! 0,1 nrows
    local_col_row_info(0,1)=nrow_local
    local_col_row_info(1:nrow_local,1)=row_indices(1:nrow_local)
    ! 0,2 ncols
    local_col_row_info(0,2)=ncol_local
    local_col_row_info(1:ncol_local,2)=col_indices(1:ncol_local)

    ALLOCATE(rec_col_row_info(0:max_row_col_local,2))

    ! accumulate data on C buffer starting from myself
    DO iiB=1, nrow_local
      i_global=row_indices(iiB)
      IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          C(i_global-my_mu_start+1,j_global)=local_C(iiB,jjB)
        END DO
      END IF
    END DO

    ! start ring communication for collecting the data from the other
    proc_send_static=proc_map(para_env%mepos+1)
    proc_receive_static=proc_map(para_env%mepos-1)
    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      proc_receive=proc_map(para_env%mepos-proc_shift)

      ! first exchange information on the local data
      rec_col_row_info=0
      CALL  mp_sendrecv(local_col_row_info,proc_send_static,rec_col_row_info,proc_receive_static,para_env%group)
      nrow_rec=rec_col_row_info(0,1)
      ncol_rec=rec_col_row_info(0,2)

      ALLOCATE(row_indices_rec(nrow_rec))
      row_indices_rec=rec_col_row_info(1:nrow_rec,1)

      ALLOCATE(col_indices_rec(ncol_rec))
      col_indices_rec=rec_col_row_info(1:ncol_rec,2)

      ALLOCATE(rec_C(nrow_rec,ncol_rec))
      rec_C=0.0_dp

      ! then send and receive the real data
      CALL  mp_sendrecv(local_C,proc_send_static,rec_C,proc_receive_static,para_env%group)

      ! accumulate the received data on C buffer 
      DO iiB=1, nrow_rec
        i_global=row_indices_rec(iiB)
        IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
          DO jjB=1, ncol_rec
            j_global=col_indices_rec(jjB)
            C(i_global-my_mu_start+1,j_global)=rec_C(iiB,jjB)
          END DO
        END IF
      END DO

      local_col_row_info(:,:)=rec_col_row_info
      DEALLOCATE(local_C)
      ALLOCATE(local_C(nrow_rec,ncol_rec))
      local_C=rec_C

      DEALLOCATE(col_indices_rec)
      DEALLOCATE(row_indices_rec)
      DEALLOCATE(rec_C)
    END DO

    DEALLOCATE(local_col_row_info)
    DEALLOCATE(rec_col_row_info)
    DEALLOCATE(proc_map)

    ! proc_map, for the sub_group
    ALLOCATE(proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
    DO iiB=0,para_env_sub%num_pe-1
      proc_map(iiB)=iiB
      proc_map(-iiB-1)=para_env_sub%num_pe-iiB-1
      proc_map(para_env_sub%num_pe+iiB)=iiB
    END DO


    ! all levels
    CALL build_mat_M_munu(mp2_env,para_env_sub,mat_M_munu,C,&
                          nmo,blk,row,col,row_size,col_size,row_offset,&
                          col_offset,i_global,j_global,my_mu_start,my_mu_end,&
                          mat_munu,iter,&
                          data_block,ends_array,proc_map, &
                          sizes_array,starts_array)


    DEALLOCATE(proc_map)
    DEALLOCATE(sizes_array)
    DEALLOCATE(starts_array)
    DEALLOCATE(ends_array)
    DEALLOCATE(local_C)

    CALL timestop(handle)

  END SUBROUTINE replicate_mat_to_subgroup_im_time

! *****************************************************************************
!> \brief fill dbcsr matrix mat_M_munu
!> \param mp2_env ...
!> \param para_env_sub ...
!> \param mat_M_munu ...
!> \param Cread ...
!> \param number_of_level ...
!> \param blk ...
!> \param row ...
!> \param col ...
!> \param row_size ...
!> \param col_size ...
!> \param row_offset ...
!> \param col_offset ...
!> \param i_global ...
!> \param j_global ...
!> \param my_mu_start ...
!> \param my_mu_end ...
!> \param mat_munu ...
!> \param iter ...
!> \param data_block ...
!> \param ends_array ...
!> \param proc_map ...
!> \param sizes_array ...
!> \param starts_array ...
!> \author Jan Wilhelm
! *****************************************************************************
  SUBROUTINE build_mat_M_munu(mp2_env,para_env_sub,mat_M_munu,Cread,&
                                    number_of_level,blk,row,col,row_size,col_size,row_offset,&
                                    col_offset,i_global,j_global,my_mu_start,my_mu_end,&
                                    mat_munu,iter,&
                                    data_block,ends_array,proc_map, &
                                    sizes_array,starts_array)
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    TYPE(cp_dbcsr_type), POINTER             :: mat_M_munu
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Cread
    INTEGER :: number_of_level, blk, row, col, row_size, col_size, &
      row_offset, col_offset, i_global, j_global, my_mu_start, my_mu_end
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(cp_dbcsr_iterator)                  :: iter
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_array, proc_map, &
                                                sizes_array, starts_array

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_mat_M_munu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, proc_receive, &
                                                proc_send, proc_shift, &
                                                rec_mu_end, rec_mu_size, &
                                                rec_mu_start
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rec_C

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_copy(mat_M_munu,mat_munu%matrix)
    CALL cp_dbcsr_set(mat_M_munu,0.0_dp)


    ! accumulate data on mat_M_munu starting from myself
    CALL cp_dbcsr_iterator_start(iter, mat_M_munu)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO i=1,row_size
         i_global=row_offset+i-1
         IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
           DO j=1,col_size
             j_global=col_offset+j-1
             data_block(i,j)=Cread(i_global-my_mu_start+1,col_offset+j-1)
           ENDDO
         END IF
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    ! start ring communication in the subgroup for collecting the data from the other
    ! proc (occupied)
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=proc_map(para_env_sub%mepos-proc_shift)

      rec_mu_start=starts_array(proc_receive)
      rec_mu_end=ends_array(proc_receive)
      rec_mu_size=sizes_array(proc_receive)

      ALLOCATE(rec_C(rec_mu_size,number_of_level))
      rec_C=0.0_dp

      ! then send and receive the real data
      CALL  mp_sendrecv(Cread,proc_send,rec_C,proc_receive,para_env_sub%group)

      ! accumulate data on mat_M_munu the data received from proc_rec
      CALL cp_dbcsr_iterator_start(iter, mat_M_munu)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
              row_size=row_size, col_size=col_size, &
              row_offset=row_offset, col_offset=col_offset)
         DO i=1,row_size
           i_global=row_offset+i-1
           IF(i_global>=rec_mu_start.AND.i_global<=rec_mu_end) THEN
             DO j=1,col_size
               j_global=col_offset+j-1
               data_block(i,j)=rec_C(i_global-rec_mu_start+1,col_offset+j-1)
             ENDDO
           END IF
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      DEALLOCATE(rec_C)

    END DO
    CALL cp_dbcsr_filter(mat_M_munu,mp2_env%mp2_gpw%eps_filter)

    DEALLOCATE(Cread)

    CALL timestop(handle)

  END SUBROUTINE build_mat_M_munu  


! *****************************************************************************
!> \brief dbcsr_3index_trace
!> \param trace_2d ...
!> \param mat_occ ...
!> \param mat_virt ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ngroup ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
! *****************************************************************************
    SUBROUTINE dbcsr_3index_trace(trace_2d,mat_occ,mat_virt,para_env,para_env_sub,& 
                                  ngroup,my_group_L_size,my_group_L_start)

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: trace_2d
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_occ, mat_virt
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: ngroup, my_group_L_size, &
                                                my_group_L_start

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_3index_trace', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, col_occ, col_size, col_virt, comm_exch, exch_color, &
      group_iter, handle, handle2, i_block, i_block_occ, i_block_virt, KKK, &
      LLL, my_group_K_start, n_blocks_local_occ, n_blocks_local_virt, &
      n_entries_local_occ, n_entries_local_virt, offset, offset_occ, &
      offset_virt, proc_rec_static, proc_send_static, row, row_occ, row_size, &
      row_virt, size_occ, size_virt
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_rec, proc_send
    INTEGER, DIMENSION(:, :), POINTER        :: buffer_info_occ, &
                                                buffer_info_virt, req_array
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer_data_occ, &
                                                buffer_data_virt
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    REAL(KIND=dp), EXTERNAL                  :: DDOT
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_para_env_type), POINTER          :: para_env_exch
    TYPE(integ_mat_buffer_type), ALLOCATABLE :: buffer_send
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec

    CALL timeset(routineN,handle)

    NULLIFY(data_block, para_env_exch)

    n_blocks_local_occ = 0
    n_blocks_local_virt = 0

    ! get number of local blocks of mat_M_occ and mat_M_virt
    DO LLL = 1,my_group_L_size

      n_blocks_local_occ = n_blocks_local_occ + cp_dbcsr_get_num_blocks(mat_occ(LLL)%matrix)
      n_blocks_local_virt = n_blocks_local_virt + cp_dbcsr_get_num_blocks(mat_virt(LLL)%matrix)

    END DO

    ! get number of entries in mat_M_occ and mat_M_virt
    n_entries_local_occ = 0
    n_entries_local_virt = 0

    DO LLL = 1,my_group_L_size

      ! n_entries_local_occ
      CALL cp_dbcsr_iterator_start(iter, mat_occ(LLL)%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))

        CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                          row_size=row_size, col_size=col_size)

        n_entries_local_occ = n_entries_local_occ + row_size*col_size

      END DO
      CALL cp_dbcsr_iterator_stop(iter)

      ! n_entries_local_virt
      CALL cp_dbcsr_iterator_start(iter, mat_virt(LLL)%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))

        CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                          row_size=row_size, col_size=col_size)

        n_entries_local_virt = n_entries_local_virt + row_size*col_size

      END DO
      CALL cp_dbcsr_iterator_stop(iter)

    END DO

    ALLOCATE(buffer_info_occ(n_blocks_local_occ,5))

    ALLOCATE(buffer_data_occ(n_entries_local_occ))

    ! 1) fill buffer_data_occ with blocks from mat_occ(LLL)%matrix
    ! 2) fill block_info_occ with LLL-counter (1), row/col in original dbcsr matrix (2/3), 
    !    block_size = row_size*col_size (4) and offset in buffer (5)

    offset = 0
    i_block = 1

    DO LLL = 1,my_group_L_size
      CALL cp_dbcsr_iterator_start(iter, mat_occ(LLL)%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
  
        CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                          row_size=row_size, col_size=col_size)
  
        ! write 1D array of data_block into buffer_data_occ
        buffer_data_occ(offset+1:offset+row_size*col_size) = RESHAPE(data_block, (/row_size*col_size/))
  
        buffer_info_occ(i_block,1) = LLL
        buffer_info_occ(i_block,2) = row
        buffer_info_occ(i_block,3) = col
        buffer_info_occ(i_block,4) = row_size*col_size
        buffer_info_occ(i_block,5) = offset
  
        offset  = offset + row_size*col_size
        i_block = i_block + 1
 
      END DO ! dbcsr iterator
      CALL cp_dbcsr_iterator_stop(iter)

      CALL cp_dbcsr_set(mat_occ(LLL)%matrix,0.0_dp)
      ! deallocate the blocks in mat_occ
      CALL cp_dbcsr_filter(mat_occ(LLL)%matrix,1.0_dp)

    END DO ! LLL


    ! the same for the virtuals
    ALLOCATE(buffer_info_virt(n_blocks_local_virt,6))
    ALLOCATE(buffer_data_virt(n_entries_local_virt))

    offset = 0
    i_block = 1

    DO LLL = 1,my_group_L_size
      CALL cp_dbcsr_iterator_start(iter, mat_virt(LLL)%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))

        CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                          row_size=row_size, col_size=col_size)

        ! write 1D array of data_block into buffer_data_virt
        buffer_data_virt(offset+1:offset+row_size*col_size) = RESHAPE(data_block, (/row_size*col_size/))

        buffer_info_virt(i_block,1) = LLL
        buffer_info_virt(i_block,2) = row
        buffer_info_virt(i_block,3) = col
        buffer_info_virt(i_block,4) = row_size*col_size
        buffer_info_virt(i_block,5) = offset
        buffer_info_virt(i_block,6) = my_group_L_start

        offset  = offset + row_size*col_size
        i_block = i_block + 1

      END DO ! dbcsr iterator
      CALL cp_dbcsr_iterator_stop(iter)

      CALL cp_dbcsr_set(mat_virt(LLL)%matrix,0.0_dp)
      ! deallocate the blocks in mat_virt
      CALL cp_dbcsr_filter(mat_virt(LLL)%matrix,1.0_dp)

    END DO ! LLL

    CALL timeset(routineN//"_local_DDOT",handle2)

    ! check locally which data is there to multiply
    ! JW write a subroutine for that!!!
    DO i_block_virt = 1,n_blocks_local_virt

      KKK              = buffer_info_virt(i_block_virt,1)
      row_virt         = buffer_info_virt(i_block_virt,2)
      col_virt         = buffer_info_virt(i_block_virt,3)
      size_virt        = buffer_info_virt(i_block_virt,4)
      offset_virt      = buffer_info_virt(i_block_virt,5)
      my_group_K_start = buffer_info_virt(i_block_virt,6)


      DO i_block_occ = 1,n_blocks_local_occ

        row_occ = buffer_info_occ(i_block_occ,2)
        IF(row_occ /= row_virt) CYCLE

        col_occ = buffer_info_occ(i_block_occ,3)
        IF(col_occ /= col_virt) CYCLE

        LLL        = buffer_info_occ(i_block_occ,1)
        size_occ   = buffer_info_occ(i_block_occ,4)
        offset_occ = buffer_info_occ(i_block_occ,5)

        CPASSERT(size_occ==size_virt)

        trace_2d(LLL,KKK+my_group_K_start-1) = trace_2d(LLL,KKK+my_group_K_start-1) + &
                                               DDOT(size_occ,buffer_data_occ(offset_occ+1:offset_occ+size_occ),1,&
                                                     buffer_data_virt(offset_virt+1:offset_virt+size_virt),1)

      END DO ! i_block_occ

    END DO ! i_block_virt

    CALL timestop(handle2)

    ! Get process, from which we get and to whom we send the blocks
    exch_color = MODULO(para_env%mepos,para_env%num_pe/ngroup)
    CALL mp_comm_split_direct(para_env%group,comm_exch,exch_color)
    NULLIFY(para_env_exch)
    CALL cp_para_env_create(para_env_exch,comm_exch)

    proc_rec_static = para_env_exch%mepos + 1
    IF(proc_rec_static >= para_env_exch%num_pe) proc_rec_static = proc_rec_static - para_env_exch%num_pe 

    proc_send_static = para_env_exch%mepos - 1
    IF(proc_send_static < 0) proc_send_static = proc_send_static + para_env_exch%num_pe

    ALLOCATE(proc_rec(ngroup-1))
    ALLOCATE(proc_send(ngroup-1))

    DO group_iter=1,ngroup-1

      proc_rec(group_iter) = para_env_exch%mepos + group_iter
      IF(proc_rec(group_iter) >= para_env_exch%num_pe) THEN
        proc_rec(group_iter) = proc_rec(group_iter) - para_env_exch%num_pe
      END IF

      proc_send(group_iter) = para_env_exch%mepos - group_iter
      IF(proc_send(group_iter) < 0) THEN
        proc_send(group_iter) = proc_send(group_iter) + para_env_exch%num_pe
      END IF

    END DO

    CALL timeset(routineN//"_comm_1",handle2)

    ALLOCATE(req_array(ngroup-1,2))
    ALLOCATE(buffer_rec(ngroup-1))
    ALLOCATE(buffer_send)

    DO group_iter=1,ngroup-1
      ALLOCATE(buffer_rec(group_iter)%sizes(2))
    END DO
    ALLOCATE(buffer_send%sizes(2))

    ! get and send number of local virtual blocks in advance with  non-blocking communication
    buffer_send%sizes(1) = n_blocks_local_virt
    buffer_send%sizes(2) = n_entries_local_virt

    DO group_iter=1,ngroup-1
      CALL mp_irecv(buffer_rec(group_iter)%sizes,proc_rec(group_iter),para_env_exch%group,req_array(group_iter,2),tag=group_iter)
      CALL mp_isend(buffer_send%sizes,proc_send(group_iter),para_env_exch%group,req_array(group_iter,1),tag=group_iter)
    END DO
    CALL mp_waitall(req_array(:,:))

    CALL timestop(handle2)


    DO group_iter=1,ngroup-1

      ALLOCATE(buffer_rec(group_iter)%indx(buffer_rec(group_iter)%sizes(1),6))
      ALLOCATE(buffer_rec(group_iter)%msg(buffer_rec(group_iter)%sizes(2)))

      ! get infos
      CALL timeset(routineN//"_comm_2",handle2)
      CALL mp_irecv(buffer_rec(group_iter)%indx,proc_rec(group_iter),para_env_exch%group,req_array(group_iter,2),tag=group_iter)
      CALL mp_isend(buffer_info_virt,proc_send(group_iter),para_env_exch%group,req_array(group_iter,1),tag=group_iter)
      CALL timestop(handle2)

      ! wait (makes it blocking send/recv)
      CALL timeset(routineN//"_wait",handle2)
      CALL mp_waitall(req_array(group_iter,:))
      CALL timestop(handle2)

      ! get actual data
      CALL timeset(routineN//"_comm_3",handle2)
      CALL mp_irecv(buffer_rec(group_iter)%msg,proc_rec(group_iter),para_env_exch%group,req_array(group_iter,2),tag=group_iter)
      CALL mp_isend(buffer_data_virt,proc_send(group_iter),para_env_exch%group,req_array(group_iter,1),tag=group_iter)
      CALL timestop(handle2)

      ! wait (makes it blocking send/recv)
      CALL timeset(routineN//"_wait",handle2)
      CALL mp_waitall(req_array(group_iter,:))
      CALL timestop(handle2)

      CALL timeset(routineN//"_DDOT",handle2)

      ! perform local multiplications
      DO i_block_virt = 1,buffer_rec(group_iter)%sizes(1)
  
        KKK              = buffer_rec(group_iter)%indx(i_block_virt,1)
        row_virt         = buffer_rec(group_iter)%indx(i_block_virt,2)
        col_virt         = buffer_rec(group_iter)%indx(i_block_virt,3)
        size_virt        = buffer_rec(group_iter)%indx(i_block_virt,4)
        offset_virt      = buffer_rec(group_iter)%indx(i_block_virt,5)
        my_group_K_start = buffer_rec(group_iter)%indx(i_block_virt,6)

        DO i_block_occ = 1,n_blocks_local_occ
 
          row_occ = buffer_info_occ(i_block_occ,2)
          IF(row_occ /= row_virt) CYCLE
  
          col_occ = buffer_info_occ(i_block_occ,3)
          IF(col_occ /= col_virt) CYCLE
  
          LLL              = buffer_info_occ(i_block_occ,1)
          size_occ         = buffer_info_occ(i_block_occ,4)
          offset_occ       = buffer_info_occ(i_block_occ,5)
 
          CPASSERT(size_occ==size_virt)
 
          trace_2d(LLL,KKK+my_group_K_start-1) = trace_2d(LLL,KKK+my_group_K_start-1) + &
                                                 DDOT(size_occ,buffer_data_occ(offset_occ+1:offset_occ+size_occ),1,&
                                                      buffer_rec(group_iter)%msg(offset_virt+1:offset_virt+size_virt),1)
        END DO ! i_block_occ
  
      END DO ! i_block_virt

      DEALLOCATE(buffer_rec(group_iter)%indx)

      DEALLOCATE(buffer_rec(group_iter)%msg)

      CALL timestop(handle2)

    END DO ! communication in groups

    CALL timeset(routineN//"_comm_trace",handle2)

    ! para_env_sub is para env of the underlying dbcsr matrix
    CALL mp_sum(trace_2d(:,:),para_env_sub%group)

    CALL timestop(handle2)

    CALL timeset(routineN//"_sync",handle2)

    CALL mp_sync(para_env_exch%group)

    CALL timestop(handle2)

    CALL cp_para_env_release(para_env_exch)

    DEALLOCATE(buffer_info_occ,buffer_info_virt)
    DEALLOCATE(buffer_data_occ,buffer_data_virt)
    DEALLOCATE(proc_send,proc_rec)
    DEALLOCATE(req_array)
    DO group_iter=1,ngroup-1
      DEALLOCATE(buffer_rec(group_iter)%sizes)
    END DO
    DEALLOCATE(buffer_send%sizes)
    DEALLOCATE(buffer_rec,buffer_send)

    CALL timestop(handle)

    END SUBROUTINE


! *****************************************************************************
!> \brief fill_fm_with_2d_array
!> \param fm_mat_Q ...
!> \param trace_2d ...
!> \param para_env ...
!> \param my_group_L_start ...
!> \param my_group_L_size ...
!> \param dimen_RI ...
!> \param ngroup ...
! *****************************************************************************
    SUBROUTINE fill_fm_with_2d_array(fm_mat_Q,trace_2d,para_env,my_group_L_start,my_group_L_size,&
                                     dimen_RI,ngroup)

    TYPE(cp_fm_type), POINTER                :: fm_mat_Q
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: trace_2d
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: my_group_L_start, &
                                                my_group_L_size, dimen_RI, &
                                                ngroup

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_with_2d_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: igroup, my_group_K_size, &
                                                my_group_K_start, proc_rec_M, &
                                                proc_send_M
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: my_group_K_rec, &
                                                my_group_K_send
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: trace_2d_rec, trace_2d_send

! set the local data

    CALL cp_fm_set_submatrix(fm=fm_mat_Q,new_values=trace_2d,start_row=my_group_L_start,start_col=1,&
                             n_rows=my_group_L_size, n_cols=dimen_RI)      


    ! Get process, from which we get and to whom we send trace2d

    proc_send_M = para_env%mepos + para_env%num_pe/ngroup

    IF(proc_send_M >= para_env%num_pe) THEN
      proc_send_M = proc_send_M - para_env%num_pe
    END IF

    proc_rec_M = para_env%mepos - para_env%num_pe/ngroup

    IF(proc_rec_M < 0) THEN
      proc_rec_M = proc_rec_M + para_env%num_pe
    END IF

    ALLOCATE(my_group_K_send(2))
    ALLOCATE(my_group_K_rec(2))
    my_group_K_send(1) = my_group_L_start
    my_group_K_send(2) = my_group_L_size

    ALLOCATE(trace_2d_send(my_group_L_size,dimen_RI))
    
    trace_2d_send = trace_2d

    DO igroup=1,ngroup

      ! send and receive my_group_K_start and my_group_K_size
      CALL mp_sendrecv(my_group_K_send,proc_send_M,my_group_K_rec,proc_rec_M,para_env%group)

      my_group_K_start = my_group_K_rec(1)
      my_group_K_size  = my_group_K_rec(2)
      
      ALLOCATE(trace_2d_rec(my_group_K_size,dimen_RI))

      CALL mp_sendrecv(trace_2d_send,proc_send_M,trace_2d_rec,proc_rec_M,para_env%group)

      ! send and receive trace_2d
      CALL cp_fm_set_submatrix(fm=fm_mat_Q,new_values=trace_2d_rec,start_row=my_group_K_start,start_col=1,&
                               n_rows=my_group_K_size, n_cols=dimen_RI)    

      DEALLOCATE(trace_2d_send)

      ALLOCATE(trace_2d_send(my_group_K_size,dimen_RI))

      trace_2d_send = trace_2d_rec

      DEALLOCATE(trace_2d_rec)

      my_group_K_send(1) = my_group_K_start
      my_group_K_send(2) = my_group_K_size

    END DO

    DEALLOCATE(trace_2d_send)

    DEALLOCATE(my_group_K_send,my_group_K_rec)

    END SUBROUTINE fill_fm_with_2d_array

END MODULE rpa_im_time 
